
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
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
#include <unordered_map>
#include <vector>

// Include local modules
#include "../../config.hpp"
#include "../../containers/circular_buffer.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/source_node.hpp"
#include "../../interfaces/gwtfcns_interface_t.hpp"
#include "../../math/integration.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../waves/regular_wave_fo.hpp"
#include "input_t.hpp"

// MPI
#include "mpi.h"


// ---------------------------------------------------------------------------
// Stable per-face σ history key.
//
// The Duhamel convolution must reach back into σ values stored at past time
// steps.  The naive layout (a per-BEM-panel-index vector) is broken by the
// wetted-surface remesher: panel index `i` at step T is not the same physical
// face as index `i` at step T-1 once a panel emerges from or submerges into
// the water — geometry from the current step gets multiplied by σ from a
// different face.
//
// Instead we key σ history by the underlying physical face, encoded as a
// single 64-bit integer:
//
//     key = (int64_t(parent_body_id) << 32) | uint32_t(parent_face_id)
//
// Sub-panels of an FS-intersecting parent face share the parent's
// (body, face) id, so they aggregate naturally into a single per-face σ.
// ---------------------------------------------------------------------------
using SigmaFaceKey = std::int64_t;

inline SigmaFaceKey make_face_key( int body_id, int face_id )
{
    return ( static_cast<std::int64_t>( body_id ) << 32 )
         |   static_cast<std::uint32_t>( face_id );
}

/// One σ snapshot per physical face, indexed by SigmaFaceKey.
using SigmaFaceMap  = std::unordered_map<SigmaFaceKey, cusfloat>;

/// Circular history of σ snapshots — `_sigma_hist[k]` is the snapshot k steps ago.
using SigmaHistType = CircularBuffer<SigmaFaceMap>;


/**
 * @brief Time-domain BEM formulation kernel backend.
 *
 * Manages the steady BEM system matrix (real-valued, MPI-parallel via ScaLAPACK)
 * and the Duhamel convolution over the source intensity history.
 *
 * @tparam N    Spatial Gauss points per edge (= NUM_GP = 2).
 * @tparam NGPT Number of time Gauss points for the Duhamel time integration.
 */
template<std::size_t N, int NGPT>
struct FormulationKernelBackendT
{
private:
    // ScaLAPACK distributed solver (real, MPI-parallel)
    SclReal*    _solver             = nullptr;

    // Clean steady sysmat (never passed to Solve; used to restore working copy)
    cusfloat*   _sysmat_steady      = nullptr;

    // Working sysmat local block (gets overwritten with LU factors by pgesv)
    cusfloat*   _sysmat             = nullptr;

    // RHS and solution vectors (size = n_panels, full vector on every process)
    cusfloat*   _rhs                = nullptr;  ///< Normal-velocity BC RHS  (sigma solve)
    cusfloat*   _rhs_dt             = nullptr;  ///< Time-derivative BC RHS  (sigma_dt solve)
    cusfloat*   _sigma              = nullptr;
    cusfloat*   _sigma_dt           = nullptr;

    // Decomposed RHS contributions (debug / ParaView export)
    cusfloat*   _rhs_body_kin       = nullptr;  ///< Body kinematic BC part only
    cusfloat*   _rhs_duhamel        = nullptr;  ///< Duhamel convolution part only
    cusfloat*   _rhs_wave           = nullptr;  ///< Incident-wave diffraction BC part only

    // Total velocity-potential derivatives (radiation + incident wave; sum of the two
    // split arrays below).  Updated at the end of compute_potential_derivatives().
    cusfloat*   _phi_dt         = nullptr;  ///< d(phi_total)/dt  at each panel
    cusfloat*   _phi_dx         = nullptr;  ///< d(phi_total)/dx  at each panel
    cusfloat*   _phi_dy         = nullptr;  ///< d(phi_total)/dy  at each panel
    cusfloat*   _phi_dz         = nullptr;  ///< d(phi_total)/dz  at each panel

    // Split velocity-potential derivatives: radiation contribution
    // (Duhamel convolution memory kernel + steady Rankine, both σ-driven)
    cusfloat*   _phi_dt_rad     = nullptr;  ///< d(phi_rad)/dt  at each panel (static + memory)
    cusfloat*   _phi_dx_rad     = nullptr;  ///< d(phi_rad)/dx  at each panel (static + memory)
    cusfloat*   _phi_dy_rad     = nullptr;  ///< d(phi_rad)/dy  at each panel (static + memory)
    cusfloat*   _phi_dz_rad     = nullptr;  ///< d(phi_rad)/dz  at each panel (static + memory)

    // Sub-split of the radiation contribution:
    //   _static = steady-Rankine (σ̇·G₀) part assembled in compute_potential_derivatives
    //   _memory = Duhamel free-surface convolution part assembled in build_rhs
    // The sum equals _phi_dt_rad / _phi_dx_rad / _phi_dy_rad / _phi_dz_rad.
    cusfloat*   _phi_dt_rad_static = nullptr;  ///< d(phi_rad_static)/dt at each panel
    cusfloat*   _phi_dx_rad_static = nullptr;  ///< d(phi_rad_static)/dx at each panel
    cusfloat*   _phi_dy_rad_static = nullptr;  ///< d(phi_rad_static)/dy at each panel
    cusfloat*   _phi_dz_rad_static = nullptr;  ///< d(phi_rad_static)/dz at each panel
    cusfloat*   _phi_dt_rad_memory = nullptr;  ///< d(phi_rad_memory)/dt at each panel
    cusfloat*   _phi_dx_rad_memory = nullptr;  ///< d(phi_rad_memory)/dx at each panel
    cusfloat*   _phi_dy_rad_memory = nullptr;  ///< d(phi_rad_memory)/dy at each panel
    cusfloat*   _phi_dz_rad_memory = nullptr;  ///< d(phi_rad_memory)/dz at each panel

    // Split velocity-potential derivatives: incident wave contribution
    // (analytical first-order wave potential; independent of body motion)
    cusfloat*   _phi_dt_wave    = nullptr;  ///< d(phi_wave)/dt  at each panel
    cusfloat*   _phi_dx_wave    = nullptr;  ///< d(phi_wave)/dx  at each panel
    cusfloat*   _phi_dy_wave    = nullptr;  ///< d(phi_wave)/dy  at each panel
    cusfloat*   _phi_dz_wave    = nullptr;  ///< d(phi_wave)/dz  at each panel

    // Time-domain Green's function integrator
    GWTFcnsInterfaceT<N*N>  _gwtfcns_interf;

    // Pointers to external objects (not owned)
    InputT*              _input      = nullptr;
    MeshGroup*           _mesh_gp    = nullptr;
    MpiConfig*           _mpi_config = nullptr;
    const RegularWaveFO* _inc_wave   = nullptr;   ///< incident wave (owned by TimeSolver)

    // Panel count
    int         _n_panels   = 0;

    // Current simulation time (cached from build_rhs; used by compute_potential_derivatives)
    cusfloat    _t_current  = static_cast<cusfloat>( 0.0 );

    // Preallocated scratch buffers for compute_potential_derivatives() — steady contribution.
    // Allocated once in initialize() and cleared at the start of each call to avoid
    // repeated heap allocation/deallocation every time step.
    cusfloat*   _acc_phi_dt          = nullptr;  ///< MPI partial-sum accumulator: d(phi_steady)/dt
    cusfloat*   _acc_phi_dx          = nullptr;  ///< MPI partial-sum accumulator: d(phi_steady)/dx
    cusfloat*   _acc_phi_dy          = nullptr;  ///< MPI partial-sum accumulator: d(phi_steady)/dy
    cusfloat*   _acc_phi_dz          = nullptr;  ///< MPI partial-sum accumulator: d(phi_steady)/dz
    cusfloat*   _steady_field_points = nullptr;  ///< Panel centre coordinates cache (size 3*n_panels)

    // Build the steady BEM matrix and LU-factorize it
    void _build_steady_matrix( );

    /**
     * @brief Accumulate Duhamel contributions using the trapezoidal (history-step) rule.
     *
     * For each sub-interval k, sigma is treated as constant and equal to
     * sigma_hist[k][src] (the value at the older end of the interval).
     * The spatial–temporal Green's function is integrated with the existing
     * Gauss–Legendre rule in time (NGPT points per sub-interval).
     *
     * Updates: _rhs, _rhs_dt, _rhs_duhamel, _phi_dt_rad/_dx/_dy/_dz_rad.
     * All output arrays must be pre-zeroed by the caller (build_rhs).
     * No MPI reduction is performed here; that is left to build_rhs.
     */
    void _accumulate_duhamel_trapezoidal(
                                            const SigmaHistType&    sigma_hist,
                                            cusfloat                t_current,
                                            cusfloat                dt
                                        );

    /**
     * @brief Accumulate Duhamel contributions using GK quadrature with
     *        piecewise-linear sigma interpolation across each sub-interval.
     *
     * Within each history sub-interval [t-(k+1)*dt, t-k*dt], sigma(tau) is
     * linearly interpolated between sigma_hist[k] (older end) and
     * sigma_hist[k-1] (newer end, for k≥1).  For k=0 the newest interval
     * has no available sigma at t_current, so sigma_hist[0] is used at both
     * endpoints (constant extrapolation, identical to the trapezoidal rule
     * for that single interval).
     *
     * The sigma weighting is folded directly into the time quadrature via
     * quadrature_panel_time_t_sigma, so the outputs already include sigma.
     *
     * Updates: _rhs, _rhs_dt, _rhs_duhamel, _phi_dt_rad/_dx/_dy/_dz_rad.
     * All output arrays must be pre-zeroed by the caller (build_rhs).
     * No MPI reduction is performed here; that is left to build_rhs.
     */
    void _accumulate_duhamel_gk_sigma(
                                            const SigmaHistType&    sigma_hist,
                                            cusfloat                t_current,
                                            cusfloat                dt
                                      );

public:
    // Constructor and destructor
    FormulationKernelBackendT( ) = default;
    ~FormulationKernelBackendT( );

    /**
     * @brief Initialize the backend: build steady matrix and factorize.
     *
     * @param mesh_gp    Pointer to the mesh group containing panels.
     * @param input      Pointer to the time-domain input data.
     * @param mpi_config Pointer to MPI configuration.
     */
    void    initialize( MeshGroup* mesh_gp, InputT* input, MpiConfig* mpi_config,
                        const RegularWaveFO* inc_wave );

    /**
     * @brief Build the RHS vector for time step t_current using the Duhamel convolution.
     *
     * Partial contributions are distributed over local ScaLAPACK columns
     * and summed via MPI_Allreduce.
     *
     * @param t_current     Current simulation time.
     * @param sigma_hist    Circular history of source intensities; sigma_hist[k] is the solution
     *                      at time t_current - (k+1)*dt.
     * @param dt            Time step size.
     * @param body_vel_bc   Optional pre-computed normal body velocity at each collocation
     *                      point (size = n_panels). Pass nullptr for stationary bodies.
     */
    void    build_rhs(
                        cusfloat                t_current,
                        const SigmaHistType&    sigma_hist,
                        cusfloat                dt,
                        const cusfloat*         body_vel_bc  = nullptr,
                        const cusfloat*         body_acc_bc  = nullptr,
                        const cusfloat*         body_kin_bc  = nullptr,  ///< body-motion BC only (stored in _rhs_body_kin for debug)
                        const cusfloat*         wave_bc      = nullptr   ///< wave diffraction BC only (stored in _rhs_wave for debug)
                     );

    /**
     * @brief Solve the linear system A*sigma = rhs using ScaLAPACK (pgesv),
     *        then broadcast the solution to all MPI processes.
     */
    void    solve( );

    /**
     * @brief Back-substitute the current σ̇ right-hand side (_rhs_dt) using the
     *        LU factors left by the most recent solve().
     *
     * Solves A·σ̇ = _rhs_dt reusing the stored LU (no re-factorization), then
     * broadcasts σ̇.  Leaves σ untouched.  The caller must have run build_rhs()
     * (to assemble _rhs_dt and clear the radiation arrays) and, earlier on the
     * same geometry, solve() (to factorize).  Intended for sweeps that reuse one
     * factorization across many acceleration RHSs — e.g. the 6 unit-acceleration
     * solves of the infinite-frequency added mass.
     */
    void    backsub_sigma_dt( );

    /**
     * @brief Add the incident wave and steady Rankine contributions to the
     *        velocity-potential derivative fields, completing the rad/wave split.
     *
     * After build_rhs() the radiation (Duhamel) part of dphi/dt, dphi/dx,
     * dphi/dy, dphi/dz is stored in _phi_dt_rad/_dx_rad/_dy_rad/_dz_rad.
     * This method:
     *   (a) fills _phi_dt_wave/_dx_wave/_dy_wave/_dz_wave from the analytical
     *       first-order incident wave potential;
     *   (b) adds the steady Rankine (σ × Green kernel) contribution to the
     *       radiation arrays (steady Rankine also depends on body motion);
     *   (c) sets the total arrays _phi_dt/_dx/_dy/_dz = rad + wave.
     *
     * Must be called after solve() on the same time step.
     */
    void    compute_potential_derivatives( );

    /**
     * @brief Access the current source intensity solution vector.
     */
    cusfloat*   get_sigma( )    { return this->_sigma;    }

    /**
     * @brief Access the current time-derivative of the source intensity solution.
     */
    cusfloat*   get_sigma_dt( )     { return this->_sigma_dt;     }

    /**
     * @brief Access the Duhamel integral for d(phi)/dt at each collocation point.
     *        Returns the TOTAL derivative (radiation + incident wave).
     */
    cusfloat*   get_phi_dt( )   { return this->_phi_dt;   }

    /**
     * @brief Access the Duhamel integral for d(phi)/dx at each collocation point.
     *        Returns the TOTAL derivative (radiation + incident wave).
     */
    cusfloat*   get_phi_dx( )   { return this->_phi_dx;   }

    /**
     * @brief Access the Duhamel integral for d(phi)/dy at each collocation point.
     *        Returns the TOTAL derivative (radiation + incident wave).
     */
    cusfloat*   get_phi_dy( )   { return this->_phi_dy;   }

    /**
     * @brief Access the Duhamel integral for d(phi)/dz at each collocation point.
     *        Returns the TOTAL derivative (radiation + incident wave).
     */
    cusfloat*   get_phi_dz( )   { return this->_phi_dz;   }

    // --- Radiation-only derivatives (Duhamel convolution + steady Rankine) ---

    /** @brief d(phi_rad)/dt at each panel centre (radiation contribution only). */
    cusfloat*   get_phi_dt_rad( )  { return this->_phi_dt_rad;  }
    /** @brief d(phi_rad)/dx at each panel centre (radiation contribution only). */
    cusfloat*   get_phi_dx_rad( )  { return this->_phi_dx_rad;  }
    /** @brief d(phi_rad)/dy at each panel centre (radiation contribution only). */
    cusfloat*   get_phi_dy_rad( )  { return this->_phi_dy_rad;  }
    /** @brief d(phi_rad)/dz at each panel centre (radiation contribution only). */
    cusfloat*   get_phi_dz_rad( )  { return this->_phi_dz_rad;  }

    // --- Radiation sub-split: steady-Rankine (static) vs. Duhamel (memory) ---

    cusfloat*   get_phi_dt_rad_static( ) { return this->_phi_dt_rad_static; }
    cusfloat*   get_phi_dx_rad_static( ) { return this->_phi_dx_rad_static; }
    cusfloat*   get_phi_dy_rad_static( ) { return this->_phi_dy_rad_static; }
    cusfloat*   get_phi_dz_rad_static( ) { return this->_phi_dz_rad_static; }
    cusfloat*   get_phi_dt_rad_memory( ) { return this->_phi_dt_rad_memory; }
    cusfloat*   get_phi_dx_rad_memory( ) { return this->_phi_dx_rad_memory; }
    cusfloat*   get_phi_dy_rad_memory( ) { return this->_phi_dy_rad_memory; }
    cusfloat*   get_phi_dz_rad_memory( ) { return this->_phi_dz_rad_memory; }

    // --- Incident-wave-only derivatives (analytical first-order wave potential) ---

    /** @brief d(phi_wave)/dt at each panel centre (incident wave contribution only). */
    cusfloat*   get_phi_dt_wave( ) { return this->_phi_dt_wave; }
    /** @brief d(phi_wave)/dx at each panel centre (incident wave contribution only). */
    cusfloat*   get_phi_dx_wave( ) { return this->_phi_dx_wave; }
    /** @brief d(phi_wave)/dy at each panel centre (incident wave contribution only). */
    cusfloat*   get_phi_dy_wave( ) { return this->_phi_dy_wave; }
    /** @brief d(phi_wave)/dz at each panel centre (incident wave contribution only). */
    cusfloat*   get_phi_dz_wave( ) { return this->_phi_dz_wave; }

    /**
     * @brief Return the number of panels.
     */
    int get_n_panels( ) const { return this->_n_panels; }

    // --- Decomposed RHS contributions (populated by build_rhs; useful for debugging) ---

    /** @brief Body kinematic BC contribution to the RHS (body motion only). */
    cusfloat*   get_rhs_body_kin( )  { return this->_rhs_body_kin; }
    /** @brief Duhamel convolution contribution to the RHS (memory-kernel term). */
    cusfloat*   get_rhs_duhamel( )   { return this->_rhs_duhamel;  }
    /** @brief Incident wave diffraction BC contribution to the RHS. */
    cusfloat*   get_rhs_wave( )      { return this->_rhs_wave;     }

    /**
     * @brief Export the RHS and source-intensity (sigma) fields to a VTK
     *        UnstructuredGrid (.vtu) file for ParaView visualisation.
     *
     * Writes the panel mesh as an unstructured grid (triangles / quads) with
     * two cell-data scalar arrays: "rhs" (_rhs) and "sigma" (_sigma).
     * Only MPI rank 0 writes; all other ranks return immediately.
     *
     * @param filename  Output file path (should end in .vtu).
     */
    void export_vtu( const std::string& filename ) const;
};

// Include template definitions
#include "formulation_kernel_backend_t.txx"
