
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotionsTimeDev.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

// Include general usage libraries
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// Include local modules
#include "../../config.hpp"
#include "../../containers/mpi_config.hpp"

#include "td_hdf5_exporter.hpp"
#include "../../hydrostatics.hpp"
#include "../../math/generalized_alpha.hpp"
#include "../../math/sparse/sparse_containers.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../mesh/rigid_body_mesh.hpp"
#include "formulation_kernel_backend_t.hpp"
#include "input_t.hpp"


/**
 * @brief Per-body external force functor for the time-domain Generalized-Alpha integrator.
 *
 * GeneralizedAlpha<T> requires T to be a pointer to a callable struct with the
 * 6-argument operator() signature below.  This struct stores a raw pointer to
 * the hydro_forces vector owned by TimeSolver so that the force values can be
 * updated in-place between time steps without re-assigning the functor.
 *
 * A linear ramp-up factor is applied when ramp_time > 0:
 *   scale = min( t / ramp_time, 1 )
 * This smoothly brings the external forces from zero to their full value over
 * the first ramp_time seconds, reducing impulsive transients at t = 0.
 * Set ramp_time = 0 (default) to disable the ramp entirely.
 */
struct TimeDomainFextStruct
{
    cusfloat*   forces      = nullptr;  ///< Raw ptr into TimeSolver::_hydro_forces[ib]
    int         dofs_np     = 0;        ///< Number of DOFs (= 6)
    cusfloat    ramp_time   = 0.0;      ///< Ramp-up duration [s]; 0 = disabled

    void operator()(
        cusfloat    t,
        cusfloat  /* dt  */,
        cusfloat* /* pos */,
        cusfloat* /* vel */,
        cusfloat* /* acc */,
        cusfloat* fext
    )
    {
        // Linear ramp factor: rises from 0 to 1 over [0, ramp_time].
        // When ramp_time <= 0 the factor is always 1 (no ramp).
        const cusfloat scale = ( ramp_time > static_cast<cusfloat>( 0.0 ) )
                               ? ( ( t < ramp_time ) ? ( t / ramp_time )
                                                     : static_cast<cusfloat>( 1.0 ) )
                               : static_cast<cusfloat>( 1.0 );

        for ( int id = 0; id < dofs_np; id++ )
            fext[id] = scale * forces[id];
    }
};


/**
 * @brief Main time-domain BEM solver.
 *
 * Iterates over time steps from t=0 to t=SimTime, solving the
 * radiation-diffraction problem via a source formulation with Duhamel
 * convolution for the time-domain wave Green's function.
 *
 * For free-floating bodies the structural dynamics are integrated using
 * the Generalized-Alpha method.
 *
 * Template parameters:
 * @tparam N     Spatial Gauss points per edge (= NUM_GP).
 * @tparam NGPT  Time Gauss points for the Duhamel convolution.
 */
template<std::size_t N, int NGPT>
class TimeSolver
{
private:
    // External objects (not owned)
    InputT*     _input      = nullptr;
    MpiConfig*  _mpi_config = nullptr;

    // Objects created by this solver (owned)
    MeshGroup*                          _mesh_gp        = nullptr;
    RigidBodyMesh**                     _rb_meshes      = nullptr;  ///< One RigidBodyMesh per body (owned)
    Hydrostatics**                      _hydrostatics   = nullptr;
    FormulationKernelBackendT<N, NGPT>* _kernel         = nullptr;

    // Per-body 6-DOF dynamics integrators (nullptr if is_fix == true)
    using FextFunctor = TimeDomainFextStruct*;
    std::vector<GeneralizedAlpha<FextFunctor>*>  _gen_alpha;
    std::vector<TimeDomainFextStruct>            _fext_structs;  ///< One functor struct per body

    // Source intensity history: _sigma_hist[k] is the solution k steps ago
    std::vector<std::vector<cusfloat>>  _sigma_hist;

    // Per-body CoG tracking (6-DOF: surge, sway, heave, roll, pitch, yaw)
    std::vector<std::array<cusfloat, 6>>    _body_pos;
    std::vector<std::array<cusfloat, 6>>    _body_vel;
    std::vector<std::array<cusfloat, 6>>    _body_acc;

    // Per-body hydrodynamic force vector (updated each step; referenced by _fext_structs)
    std::vector<std::vector<cusfloat>>      _hydro_forces;

    // Hydrostatic stiffness per body (6×6, row-major)
    std::vector<std::array<cusfloat, 36>>   _hydrostiff;

    // Structural mass matrix per body (6×6, diagonal)
    std::vector<std::array<cusfloat, 36>>   _struct_mass;

    // CSRMatrix objects whose arrays are referenced by the MKL sparse handles
    // inside _gen_alpha.  They must outlive the GeneralizedAlpha objects.
    std::vector<CSRMatrix*>     _csr_mass;
    std::vector<CSRMatrix*>     _csr_stiff;
    std::vector<CSRMatrix*>     _csr_damp;

    // Reusable scratch buffer for body kinematic + wave BC (avoids per-step heap allocation).
    // Sized to np each step via assign(); capacity is retained between steps.
    std::vector<cusfloat>                   _body_vel_bc;
    // Time derivative of _body_vel_bc (body acceleration + wave acceleration BC)
    std::vector<cusfloat>                   _body_acc_bc;

    // Helper methods
    void _initialize_mesh_group( );
    void _initialize_hydrostatics( );
    void _initialize_ic_positions( );
    void _apply_initial_displacement( );
    void _compute_hydrostatic_initial_forces( );
    void _initialize_structural_dynamics( );
    void _initialize_kernel( );

    /**
     * @brief Destroy and recreate MeshGroup + BEM kernel from the current
     *        state of _rb_meshes.  Called after every mesh position update.
     */
    void _rebuild_mesh_group( );

    /**
     * @brief Compute hydrodynamic pressure forces on a body from the current source intensities.
     *
     * @param body_id   Index of the body to compute forces for.
     * @param forces    Output 6-DOF force vector [Fx, Fy, Fz, Mx, My, Mz].
     */
    void _compute_hydro_forces( int body_id, cusfloat* forces );

    /**
     * @brief Accumulate gravitational force (-m*g in heave) into an existing
     *        6-DOF force vector.
     *
     * Since gravity acts at the CoG there are no gravitational moments about
     * the CoG, so only forces[2] (heave) is modified.
     *
     * @param body_id   Index of the body.
     * @param forces    6-DOF force vector to which the gravity contribution is added.
     */
    void _compute_gravitational_forces( int body_id, cusfloat* forces );

    /**
     * @brief Update mesh panel positions based on current CoG states.
     *
     * For the first approach the panels are kept at their initial positions.
     */
    void _update_mesh_positions( );

    /**
     * @brief Compute the body kinematic and incident wave boundary conditions.
     *
     * For each panel, assembles:
     *   bc[j] = n · (v_lin + omega × r)   (body motion, DIFFRAC panels only)
     *         - Re{ (grad phi_wave) · n * exp(-i*w*t) }  (incident wave diffraction)
     *
     * Also computes bc2 = d(bc)/dt:
     *   bc2[j] = n · (a_lin + alpha × r)  (body acceleration, DIFFRAC panels only)
     *          - Re{ (grad d(phi_wave)/dt) · n * exp(-i*w*t) }  (wave acceleration BC)
     *
     * @param t    Current simulation time [s].
     * @param bc   Output vector sized n_panels; values zeroed before accumulation.
     * @param bc2  Output vector sized n_panels for the time derivative of bc.
     */
    void _compute_body_vel_bc( cusfloat t, std::vector<cusfloat>& bc, std::vector<cusfloat>& bc2 );

    // ---------------------------------------------------------------
    // ParaView VTU / PVD output
    // ---------------------------------------------------------------

    /// Directory where VTU files and the PVD index are written.
    std::string _paraview_dir;

    /// Ordered list of (time [s], relative VTU filename) pairs for the PVD file.
    std::vector<std::pair<double, std::string>> _pvd_entries;

    /**
     * @brief Create the ParaView output directory and clear the PVD entry list.
     *        No-op when out_pressure is false.
     */
    void _init_paraview_output( );

    /**
     * @brief Write the PVD collection file after the simulation is complete.
     *        No-op when out_pressure is false.
     */
    void _finalize_paraview_output( );

    /**
     * @brief Write output for the current time step.
     *
     * @param t     Current simulation time.
     * @param step  Current step index.
     */
    void _output_step( cusfloat t, int step );

    // ---------------------------------------------------------------
    // HDF5 time-series output
    // ---------------------------------------------------------------

    /// Owns the HDF5 exporter object; nullptr until _init_hdf5_output() is called.
    TimeDomainHDF5Exporter* _hdf5_exporter = nullptr;

    /**
     * @brief Open the HDF5 file and initialise all extendable datasets.
     *        No-op when out_pressure is false or HDF5 is not compiled in.
     */
    void _init_hdf5_output( );

    /**
     * @brief Flush and close the HDF5 file after the simulation is complete.
     *        No-op when out_pressure is false or HDF5 is not compiled in.
     */
    void _finalize_hdf5_output( );

public:
    // Constructor and destructor
    TimeSolver( InputT* input, MpiConfig* mpi_config );
    ~TimeSolver( );

    /**
     * @brief Execute the full time-domain simulation.
     */
    void run( );
};

// Include template definitions
#include "time_solver.txx"
