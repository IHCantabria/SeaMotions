
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

// Include general usage libraries
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Include local modules
#include "time_solver.hpp"
#include "../../inout/vtu.hpp"
#include "../../tools.hpp"
#include "../../math/math_constants.hpp"
#include "../../waves/wave_dispersion_base_fo.hpp"
#include "../../waves/waves_common.hpp"


/*****************************************************************************
 * Constructor / Destructor
 *****************************************************************************/

template<std::size_t N, int NGPT>
TimeSolver<N, NGPT>::TimeSolver( InputT* input, MpiConfig* mpi_config )
    : _input( input )
    , _mpi_config( mpi_config )
{
    this->_initialize_mesh_group( );
    this->_initialize_hydrostatics( );
    this->_initialize_structural_dynamics( );
    this->_initialize_kernel( );
}


template<std::size_t N, int NGPT>
TimeSolver<N, NGPT>::~TimeSolver( )
{
    // Clean up kernel
    if ( this->_kernel != nullptr )
    {
        delete this->_kernel;
        this->_kernel = nullptr;
    }

    // Clean up hydrostatics
    if ( this->_hydrostatics != nullptr )
    {
        for ( int i=0; i<this->_input->bodies_np; i++ )
        {
            if ( this->_hydrostatics[i] != nullptr )
            {
                delete this->_hydrostatics[i];
                this->_hydrostatics[i] = nullptr;
            }
        }
        delete [] this->_hydrostatics;
        this->_hydrostatics = nullptr;
    }

    // Clean up structural dynamics integrators
    for ( auto* ga : this->_gen_alpha )
    {
        delete ga;
    }
    this->_gen_alpha.clear( );

    // Clean up CSRMatrix objects (must be freed after _gen_alpha since MKL
    // sparse handles inside GA hold raw pointers into their arrays)
    for ( auto* m : this->_csr_mass  ) { delete m; }
    for ( auto* m : this->_csr_stiff ) { delete m; }
    for ( auto* m : this->_csr_damp  ) { delete m; }
    this->_csr_mass.clear( );
    this->_csr_stiff.clear( );
    this->_csr_damp.clear( );

    // Clean up mesh group
    if ( this->_mesh_gp != nullptr )
    {
        delete this->_mesh_gp;
        this->_mesh_gp = nullptr;
    }

    // Clean up rigid-body meshes (owned by TimeSolver)
    if ( this->_rb_meshes != nullptr )
    {
        for ( int i=0; i<this->_input->bodies_np; i++ )
        {
            if ( this->_rb_meshes[i] != nullptr )
            {
                delete this->_rb_meshes[i];
                this->_rb_meshes[i] = nullptr;
            }
        }
        delete [] this->_rb_meshes;
        this->_rb_meshes = nullptr;
    }

    if ( this->_hdf5_exporter != nullptr )
    {
        delete this->_hdf5_exporter;
        this->_hdf5_exporter = nullptr;
    }
}


/*****************************************************************************
 * Initialization methods
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_initialize_mesh_group( )
{
    std::cout << " -> Building rigid-body mesh group..." << std::endl;

    const int bodies_np = this->_input->bodies_np;

    // ------------------------------------------------------------------
    // For each body create a RigidBodyMesh from the pre-loaded combined
    // (total) mesh so we inherit all sub-meshes (hull + lids).  The
    // RigidBodyMesh::move() / check_underwater_panels() pair then keeps
    // the BEM panel set consistent with the current body kinematics.
    // ------------------------------------------------------------------
    this->_rb_meshes = new RigidBodyMesh*[bodies_np];

    for ( int i=0; i<bodies_np; i++ )
    {
        BodyDef* body = this->_input->bodies[i];

        // Choose source mesh: prefer combined mesh_total (hull + lids),
        // fall back to plain hull mesh.
        const Mesh* src = ( body->mesh_total != nullptr )
                          ? body->mesh_total
                          : body->mesh;

        // Nominal draft = depth of deepest mesh node below the waterplane
        // (z = 0).  Used for bookkeeping only; check_underwater_panels()
        // determines which panels are submerged from panel geometry alone.
        const cusfloat draft = ( src->z_min < 0.0 ) ? -src->z_min : 0.0;

        this->_rb_meshes[i] = new RigidBodyMesh(
                                                    *src,
                                                    body->cog,
                                                    draft
                                                );

        // Classify panels as underwater / FS-intersecting / above-water
        // and create free-surface-refined sub-panels for those that cross
        // the waterplane.
        this->_rb_meshes[i]->check_underwater_panels( );

        std::cout << "    Body " << i
                  << ": " << this->_rb_meshes[i]->uw_elems_np
                  << " UW panels + "
                  << this->_rb_meshes[i]->fs_panels_np
                  << " FS-refined panels" << std::endl;
    }

    // Build the combined MeshGroup.  Because MeshGroup now calls the
    // virtual get_elems_np() / get_panel(), the RigidBodyMesh overrides
    // are used automatically to expose only the underwater panel set.
    Mesh** rb_as_mesh = new Mesh*[bodies_np];
    for ( int i=0; i<bodies_np; i++ )
    {
        rb_as_mesh[i] = this->_rb_meshes[i];
    }

    this->_mesh_gp = new MeshGroup(
                                        rb_as_mesh,
                                        bodies_np,
                                        false       // is_wl_points: not needed for time domain
                                    );

    delete [] rb_as_mesh;

    // Define mirror panels for image sources (finite water depth)
    this->_mesh_gp->define_mirror_panels( );

    // Build source nodes 1-to-1 with the filtered underwater panel set
    this->_mesh_gp->define_source_nodes( 0 );
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_initialize_hydrostatics( )
{
    std::cout << " -> Calculating hydrostatics..." << std::endl;

    const int bodies_np = this->_input->bodies_np;

    this->_hydrostatics = new Hydrostatics*[bodies_np];
    for ( int i=0; i<bodies_np; i++ )
    {
        this->_hydrostatics[i] = new Hydrostatics(
                                                        this->_input->bodies[i]->mesh,
                                                        this->_input->water_density,
                                                        this->_input->grav_acc,
                                                        this->_input->bodies[i]->mass,
                                                        this->_input->bodies[i]->cog,
                                                        this->_input->bodies[i]->rad_inertia,
                                                        this->_mpi_config
                                                   );
    }

    // Cache the hydrostatic stiffness matrix per body
    this->_hydrostiff.resize( bodies_np );
    for ( int i=0; i<bodies_np; i++ )
    {
        for ( int k=0; k<36; k++ )
        {
            this->_hydrostiff[i][k] = this->_hydrostatics[i]->hydstiffmat[k];
        }
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_initialize_structural_dynamics( )
{
    std::cout << " -> Initializing structural dynamics..." << std::endl;

    const int bodies_np     = this->_input->bodies_np;
    const int dofs_np       = this->_input->dofs_np;   // = 6
    const cusfloat dt       = this->_input->dt;
    const cusfloat rho_inf  = static_cast<cusfloat>( 0.8 );    // spectral radius for GA

    this->_body_pos.resize( bodies_np );
    this->_body_vel.resize( bodies_np );
    this->_body_acc.resize( bodies_np );
    this->_struct_mass.resize( bodies_np );
    this->_gen_alpha.resize( bodies_np, nullptr );

    // Initialise per-body hydro force vectors and functor structs.
    // The structs hold raw pointers into _hydro_forces so that forces can
    // be updated in-place each step without touching the GA object itself.
    this->_hydro_forces.assign( bodies_np, std::vector<cusfloat>( dofs_np, 0.0 ) );
    this->_fext_structs.resize( bodies_np );
    for ( int ib=0; ib<bodies_np; ib++ )
    {
        this->_fext_structs[ib].forces    = this->_hydro_forces[ib].data( );
        this->_fext_structs[ib].dofs_np   = dofs_np;
        this->_fext_structs[ib].ramp_time = this->_input->ramp_time;
    }

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        // Zero initial conditions
        this->_body_pos[ib].fill( 0.0 );
        this->_body_vel[ib].fill( 0.0 );
        this->_body_acc[ib].fill( 0.0 );

        if ( this->_input->bodies[ib]->is_fix )
        {
            // Fixed body: no structural integrator
            this->_gen_alpha[ib] = nullptr;
            continue;
        }

        // Build structural mass matrix (diagonal, 6 DOF)
        // Diagonal = [mass, mass, mass, Ixx, Iyy, Izz]
        cusfloat mass_dense[36];
        clear_vector( 36, mass_dense );
        const cusfloat mb   = this->_input->bodies[ib]->mass;
        const cusfloat* ria = this->_input->bodies[ib]->inertia;     // 6-component inertia
        mass_dense[ 0] = mb;                // surge
        mass_dense[ 7] = mb;                // sway
        mass_dense[14] = mb;                // heave
        mass_dense[21] = ria[0];            // roll  (Ixx)
        mass_dense[28] = ria[3];            // pitch (Iyy)
        mass_dense[35] = ria[5];            // yaw   (Izz)

        CSRMatrix* mass_mat = new CSRMatrix( dofs_np, mass_dense );

        // Hydrostatic stiffness matrix
        cusfloat stiff_dense[36];
        for ( int k=0; k<36; k++ )
        {
            stiff_dense[k] = this->_hydrostiff[ib][k];
        }
        CSRMatrix* stiff_mat = new CSRMatrix( dofs_np, stiff_dense );

        // No mechanical damping for the first approach
        cusfloat damp_dense[36];
        clear_vector( 36, damp_dense );
        CSRMatrix* damp_mat = new CSRMatrix( dofs_np, damp_dense );

        // Per-DOF kinematic restrictions read from the body definition file
        // (false = free, true = fixed per DOF: surge, sway, heave, roll, pitch, yaw).
        // Convert bool[6] to int[6] as required by GeneralizedAlpha.
        int restrictions[6];
        for ( int k=0; k<6; k++ ) restrictions[k] = this->_input->bodies[ib]->dof_restrictions[k] ? 1 : 0;

        cusfloat* y0_pos = this->_body_pos[ib].data( );
        cusfloat* y0_vel = this->_body_vel[ib].data( );
        cusfloat* y0_acc = this->_body_acc[ib].data( );

        // Pass a pointer to the functor struct (6-arg operator() required by GA)
        this->_gen_alpha[ib] = new GeneralizedAlpha<FextFunctor>(
                                                                        &this->_fext_structs[ib],
                                                                        mass_mat,
                                                                        stiff_mat,
                                                                        damp_mat,
                                                                        dt,
                                                                        static_cast<cusfloat>( 0.0 ),
                                                                        rho_inf,
                                                                        y0_pos,
                                                                        y0_vel,
                                                                        y0_acc,
                                                                        restrictions
                                                                    );

        // MKL sparse matrix handles inside GeneralizedAlpha hold raw pointers
        // into the CSRMatrix arrays.  Keep the CSRMatrix objects alive for the
        // lifetime of the solver; they are deleted in ~TimeSolver().
        this->_csr_mass.push_back( mass_mat );
        this->_csr_stiff.push_back( stiff_mat );
        this->_csr_damp.push_back( damp_mat );
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_initialize_kernel( )
{
    std::cout << " -> Building steady BEM system matrix..." << std::endl;

    this->_kernel = new FormulationKernelBackendT<N, NGPT>( );
    this->_kernel->initialize( this->_mesh_gp, this->_input, this->_mpi_config );
}


/*****************************************************************************
 * Mesh group rebuild (after rigid-body motion)
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_rebuild_mesh_group( )
{
    const int bodies_np = this->_input->bodies_np;

    // ------------------------------------------------------------------
    // The panel topology may have changed (new FS-refined panels, panels
    // emerging from / submerging into the water).  Destroy the existing
    // kernel and MeshGroup before reconstructing them.
    // ------------------------------------------------------------------
    if ( this->_kernel != nullptr )
    {
        delete this->_kernel;
        this->_kernel = nullptr;
    }
    if ( this->_mesh_gp != nullptr )
    {
        delete this->_mesh_gp;
        this->_mesh_gp = nullptr;
    }

    // Rebuild MeshGroup from the current RigidBodyMesh state.
    // Virtual dispatch in MeshGroup will pick up the filtered underwater
    // panel set from each RigidBodyMesh.
    Mesh** rb_as_mesh = new Mesh*[bodies_np];
    for ( int i=0; i<bodies_np; i++ )
    {
        rb_as_mesh[i] = this->_rb_meshes[i];
    }

    this->_mesh_gp = new MeshGroup(
                                        rb_as_mesh,
                                        bodies_np,
                                        false
                                    );
    delete [] rb_as_mesh;

    this->_mesh_gp->define_mirror_panels( );

    // Rebuild source nodes for the new panel geometry
    this->_mesh_gp->define_source_nodes( 0 );

    // Rebuild the BEM steady matrix for the new panel geometry.
    this->_kernel = new FormulationKernelBackendT<N, NGPT>( );
    this->_kernel->initialize( this->_mesh_gp, this->_input, this->_mpi_config );
}


/*****************************************************************************
 * Hydrodynamic force computation
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_hydro_forces( int body_id, cusfloat* forces )
{
    // Compute time-domain radiation + diffraction pressure from sigma.
    //
    // The pressure on a panel is related to the time derivative of the potential:
    //   p = -rho * d(phi)/dt
    // For the source formulation in the time domain, phi(t) = integral sigma(tau)*G(t-tau) dtau.
    //
    // At the current time step, the instantaneous contribution to the pressure
    // from the current sigma is approximated as:
    //   delta_phi = sigma[j] * G(t=0) * area_j  (per panel)
    //
    // For a first implementation we compute forces by integrating the panel
    // source intensity as a surface pressure proxy:
    //   F_dof = rho * sum_panels{ sigma[i] * area_i * n_i[dof] }
    //
    // This corresponds to the added-mass / radiation-force contribution at the
    // current instant without the convolution history (which is handled via the
    // RHS Duhamel integral in the BEM solve).

    const int dofs_np           = this->_input->dofs_np;  // 6
    const cusfloat rho          = this->_input->water_density;
    const cusfloat g            = this->_input->grav_acc;

    // Zero the forces vector
    clear_vector( dofs_np, forces );

    // Panels owned by this body
    const int panel_start       = this->_mesh_gp->panels_cnp[body_id];
    const int panel_end         = this->_mesh_gp->panels_cnp[body_id + 1];

    // Velocity-potential derivatives at panel centres (computed by kernel)
    const cusfloat* phi_dt  = this->_kernel->get_phi_dt( );
    const cusfloat* phi_dx  = this->_kernel->get_phi_dx( );
    const cusfloat* phi_dy  = this->_kernel->get_phi_dy( );
    const cusfloat* phi_dz  = this->_kernel->get_phi_dz( );

    // Centre of gravity of this body (for moment arm)
    const cusfloat* cog = this->_input->bodies[body_id]->cog;

    for ( int ip=panel_start; ip<panel_end; ip++ )
    {
        PanelGeom* panel = this->_mesh_gp->panels[ip];

        // Only submerged panels contribute to hydrodynamic pressure
        if ( panel->center[2] >= static_cast<cusfloat>( 0.0 ) )
            continue;

        // Bernoulli pressure (linearised Bernoulli + hydrostatic term):
        //   p = -rho * ( dphi/dt + 0.5 * |grad phi|^2 + g * z )
        const cusfloat press = -rho * (   phi_dt[ip]
                                        + static_cast<cusfloat>( 0.5 ) * (   pow2s( phi_dx[ip] )
                                                                            + pow2s( phi_dy[ip] )
                                                                            + pow2s( phi_dz[ip] ) )
                                        + g * panel->center[2] );

        // Force: F_i = p_i * area_i * n_i
        const cusfloat dF0 = press * panel->area * panel->normal_vec[0];
        const cusfloat dF1 = press * panel->area * panel->normal_vec[1];
        const cusfloat dF2 = press * panel->area * panel->normal_vec[2];

        forces[0] += dF0;
        forces[1] += dF1;
        forces[2] += dF2;

        // Moments: M = r × F   where r = panel_centre - CoG
        const cusfloat rx = panel->center[0] - cog[0];
        const cusfloat ry = panel->center[1] - cog[1];
        const cusfloat rz = panel->center[2] - cog[2];

        forces[3] += ry * dF2 - rz * dF1;  // roll
        forces[4] += rz * dF0 - rx * dF2;  // pitch
        forces[5] += rx * dF1 - ry * dF0;  // yaw
    }
}


/*****************************************************************************
 * Mesh update
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_body_vel_bc( 
                                                    cusfloat t, 
                                                    std::vector<cusfloat>& bc,
                                                    std::vector<cusfloat>& bc2
                                                )
{
    const int np = this->_kernel->get_n_panels( );
    std::cout << "  [DBG _vel_bc] entry np=" << np << "\n" << std::flush;
    bc.assign( np, static_cast<cusfloat>( 0.0 ) );
    bc2.assign( np, static_cast<cusfloat>( 0.0 ) );

    const cusfloat w  = this->_input->ang_freq;
    const cusfloat aw = this->_input->wave_amp;
    const cusfloat h  = this->_input->water_depth;
    const cusfloat g  = this->_input->grav_acc;

    // Convert heading from degrees to radians
    const cusfloat mu = this->_input->wave_heading * PI / static_cast<cusfloat>( 180.0 );

    const bool has_wave = ( aw > static_cast<cusfloat>( 0.0 ) && w > static_cast<cusfloat>( 0.0 ) );

    // Precompute wave number (only when waves are present)
    cusfloat k = static_cast<cusfloat>( 0.0 );
    if ( has_wave )
    {
        k = w2k( w, h, g );
    }

    for ( int j=0; j<np; j++ )
    {
        PanelGeom*      panel = this->_mesh_gp->panels[j];
        const cusfloat* n     = panel->normal_vec;  // 6-DOF: [0..2] physical normal, [3..5] = r×n

        if ( panel->type != PanelTypeE::DIFFRAC ) { continue; }

        // ------------------------------------------------------------------
        // 1. Body kinematic BC:  bc[j] += sum_{id=0}^{5} vel[id] * n[id]
        //    normal_vec[0..2] projects linear velocity; [3..5] = r×n
        //    projects rotational velocity.
        // ------------------------------------------------------------------
        const int ib = this->_mesh_gp->get_body_id( j );
        if ( ib >= 0 && !this->_input->bodies[ib]->is_fix )
        {
            const cusfloat* vel = this->_body_vel[ib].data( );
            const cusfloat* acc = this->_body_acc[ib].data( );
            for ( int id=0; id<6; id++ )
            {
                bc[j]   += vel[id] * n[id];
                bc2[j]  += acc[id] * n[id];
            }
        }

        // ------------------------------------------------------------------
        // 2. Incident wave diffraction BC (negative of wave normal velocity):
        //    bc[j] -= Re{ (phi_dx*n[0] + phi_dy*n[1] + phi_dz*n[2])
        //                  * exp(-i*w*t) }
        //    where Re{ z * exp(-i*w*t) } = Re{z}*cos(wt) + Im{z}*sin(wt)
        // ------------------------------------------------------------------
        if ( has_wave )
        {
            const cusfloat x        = panel->center[0];
            const cusfloat y        = panel->center[1];
            const cusfloat z        = panel->center[2];

            const cusfloat vx       = wave_potential_fo_time_dx( aw, w, k, h, g, x, y, z, mu, t );
            const cusfloat vy       = wave_potential_fo_time_dy( aw, w, k, h, g, x, y, z, mu, t );
            const cusfloat vz       = wave_potential_fo_time_dz( aw, w, k, h, g, x, y, z, mu, t );

            const cusfloat vx_dt    = wave_potential_fo_time_dtdx( aw, w, k, h, g, x, y, z, mu, t );
            const cusfloat vy_dt    = wave_potential_fo_time_dtdy( aw, w, k, h, g, x, y, z, mu, t );
            const cusfloat vz_dt    = wave_potential_fo_time_dtdz( aw, w, k, h, g, x, y, z, mu, t );

            // Negative of incident wave normal velocity (diffraction condition)
            bc[j] -= ( vx * n[0] + vy * n[1] + vz * n[2] );
            bc2[j] -= ( vx_dt * n[0] + vy_dt * n[1] + vz_dt * n[2] );
        }
    }
}

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_update_mesh_positions( )
{
    std::cout << "  [DBG _update_mesh] entry\n" << std::flush;
    const int bodies_np = this->_input->bodies_np;
    bool any_moved = false;

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix ) { continue; }

        // _body_pos[ib] = { surge(x), sway(y), heave(z),
        //                   roll(rx),  pitch(ry), yaw(rz) }
        const cusfloat* pos = this->_body_pos[ib].data( );

        // Apply rigid-body motion: rotate about the CoG then translate.
        // RigidBodyMesh::move() always starts from the original (backup)
        // node positions, so accumulated floating-point errors are avoided.
        std::cout << "  [DBG _update_mesh] body " << ib << " -> move()\n" << std::flush;
        this->_rb_meshes[ib]->move(
                                        pos[0],  // surge  (dx)
                                        pos[1],  // sway   (dy)
                                        pos[2],  // heave  (dz)
                                        pos[3],  // roll   (drx)
                                        pos[4],  // pitch  (dry)
                                        pos[5]   // yaw    (drz)
                                  );

        // Reclassify panels: underwater / FS-intersecting / above and
        // regenerate the free-surface-refined panel set.
        std::cout << "  [DBG _update_mesh] body " << ib << " -> check_underwater_panels()\n" << std::flush;
        this->_rb_meshes[ib]->check_underwater_panels( );

        any_moved = true;
    }

    // Rebuild the MeshGroup and the BEM steady matrix to reflect the new
    // panel geometry across all bodies that have moved.
    if ( any_moved )
    {
        std::cout << "  [DBG _update_mesh] -> _rebuild_mesh_group()\n" << std::flush;
        this->_rebuild_mesh_group( );
    }
    std::cout << "  [DBG _update_mesh] done\n" << std::flush;
}


/*****************************************************************************
 * ParaView output helpers
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_init_paraview_output( )
{
    if ( !this->_input->out_pressure ) { return; }

    namespace fs = std::filesystem;

    const fs::path results_dir  = fs::path( this->_input->folder_path )
                                  / fs::path( RESULTS_FOLDER_NAME );
    const fs::path paraview_dir = results_dir / fs::path( RESULTS_PARAVIEW_FOLDER_NAME );

    if ( !fs::exists( results_dir ) )
    {
        fs::create_directory( results_dir );
    }
    if ( !fs::exists( paraview_dir ) )
    {
        fs::create_directory( paraview_dir );
    }

    this->_paraview_dir = paraview_dir.string( );
    this->_pvd_entries.clear( );

    std::cout << " -> ParaView pressure output: " << this->_paraview_dir << std::endl;
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_finalize_paraview_output( )
{
    if ( !this->_input->out_pressure || this->_pvd_entries.empty( ) ) { return; }

    namespace fs = std::filesystem;

    const std::string pvd_path = ( fs::path( this->_paraview_dir ) / "pressure.pvd" ).string( );
    if ( write_pvd( pvd_path, this->_pvd_entries ) )
    {
        std::cout << " -> ParaView PVD written: " << pvd_path << std::endl;
    }
}


/*****************************************************************************
 * HDF5 time-series output helpers
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_init_hdf5_output( )
{
    if ( !this->_input->out_pressure ) { return; }

    const int np = this->_kernel->get_n_panels( );
    const int nb = this->_input->bodies_np;

    this->_hdf5_exporter = new TimeDomainHDF5Exporter( );
    this->_hdf5_exporter->initialize( this->_input->folder_path, np, nb );
}

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_finalize_hdf5_output( )
{
    if ( this->_hdf5_exporter != nullptr )
    {
        this->_hdf5_exporter->close( );
        std::cout << " -> HDF5 time-series closed." << std::endl;
    }
}


/*****************************************************************************
 * Output step
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_output_step( cusfloat t, int step )
{
    // ---------------------------------------------------------------
    // Console progress (every 10 steps)
    // ---------------------------------------------------------------
    const int bodies_np = this->_input->bodies_np;
    const int dofs_np   = this->_input->dofs_np;

    if ( step % 10 == 0 )
    {
        std::cout << "  t = " << t;
        for ( int ib=0; ib<bodies_np; ib++ )
        {
            std::cout << "  body[" << ib << "] pos:";
            for ( int id=0; id<dofs_np; id++ )
            {
                std::cout << " " << this->_body_pos[ib][id];
            }
        }
        std::cout << std::endl;
    }

    // ---------------------------------------------------------------
    // Pressure field: compute and route to VTU and/or HDF5
    // ---------------------------------------------------------------
    if ( !this->_input->out_pressure ) { return; }

    const int np = this->_kernel->get_n_panels( );
    if ( np <= 0 ) { return; }

    const cusfloat  rho = this->_input->water_density;
    const cusfloat  g   = this->_input->grav_acc;

    const cusfloat* phi_dt = this->_kernel->get_phi_dt( );
    const cusfloat* phi_dx = this->_kernel->get_phi_dx( );
    const cusfloat* phi_dy = this->_kernel->get_phi_dy( );
    const cusfloat* phi_dz = this->_kernel->get_phi_dz( );

    // ---- Build geometry arrays (for VTU) and decomposed pressure components ----
    std::vector<cusfloat>   nodes_x, nodes_y, nodes_z;
    std::vector<int32_t>    connectivity, offsets;
    std::vector<uint8_t>    types;

    const std::size_t snp = static_cast<std::size_t>( np );
    std::vector<cusfloat>   phi_dt_comp     ( snp );
    std::vector<cusfloat>   kinetic_comp    ( snp );
    std::vector<cusfloat>   hydrostatic_comp( snp );
    std::vector<cusfloat>   pressure        ( snp );

    nodes_x.reserve( snp * 4 );
    nodes_y.reserve( snp * 4 );
    nodes_z.reserve( snp * 4 );
    connectivity.reserve( snp * 4 );
    offsets.reserve( snp );
    types.reserve( snp );

    int32_t node_count = 0;

    for ( int ip = 0; ip < np; ip++ )
    {
        PanelGeom* panel = this->_mesh_gp->panels[ip];
        const int  nn    = panel->num_nodes;

        for ( int k = 0; k < nn; k++ )
        {
            nodes_x.push_back( panel->x[k] );
            nodes_y.push_back( panel->y[k] );
            nodes_z.push_back( panel->z[k] );
            connectivity.push_back( node_count + k );
        }
        node_count += static_cast<int32_t>( nn );
        offsets.push_back( node_count );

        // VTK cell type: 5 = VTK_TRIANGLE, 9 = VTK_QUAD, 7 = VTK_POLYGON
        uint8_t vtype = ( nn == 3 ) ? 5 : ( nn == 4 ) ? 9 : 7;
        types.push_back( vtype );

        // Decomposed Bernoulli pressure components
        const std::size_t sip = static_cast<std::size_t>( ip );
        phi_dt_comp[sip]      = -rho * phi_dt[ip];
        kinetic_comp[sip]     = -rho * static_cast<cusfloat>( 0.5 )
                                     * (   phi_dx[ip] * phi_dx[ip]
                                         + phi_dy[ip] * phi_dy[ip]
                                         + phi_dz[ip] * phi_dz[ip] );
        hydrostatic_comp[sip] = -rho * g * panel->center[2];
        pressure[sip]         = phi_dt_comp[sip] + kinetic_comp[sip] + hydrostatic_comp[sip];
    }

    // ---- ParaView VTU: pressure_XXXXXX.vtu ----
    std::ostringstream ss;
    ss << "pressure_" << std::setfill( '0' ) << std::setw( 6 ) << step << ".vtu";
    const std::string vtu_name = ss.str( );

    namespace fs = std::filesystem;
    const std::string vtu_path = ( fs::path( this->_paraview_dir ) / vtu_name ).string( );

    write_vtu_panel_pressure(
                                vtu_path,
                                static_cast<std::size_t>( node_count ),
                                nodes_x.data( ),
                                nodes_y.data( ),
                                nodes_z.data( ),
                                snp,
                                connectivity.data( ),
                                offsets.data( ),
                                types.data( ),
                                pressure.data( ),
                                phi_dt_comp.data( ),
                                kinetic_comp.data( ),
                                hydrostatic_comp.data( )
                            );

    // Record entry for PVD (relative path from the paraview dir's parent = results dir)
    const std::string pvd_rel = ( fs::path( RESULTS_PARAVIEW_FOLDER_NAME ) / vtu_name ).string( );
    this->_pvd_entries.emplace_back( static_cast<double>( t ), pvd_rel );

    // ---- HDF5 time-series append ----
    if ( this->_hdf5_exporter != nullptr && this->_hdf5_exporter->is_open( ) )
    {
        // Integrate panel pressure components into per-body 6-DOF forces.
        const int n_bodies = this->_input->bodies_np;

        std::vector<std::array<cusfloat, 6>> body_force_total      ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_phi_dt     ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_kinetic    ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_hydrostatic( static_cast<std::size_t>( n_bodies ) );

        for ( int ib = 0; ib < n_bodies; ib++ )
        {
            body_force_total[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_phi_dt[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_kinetic[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_hydrostatic[ib].fill( static_cast<cusfloat>( 0.0 ) );

            const int       panel_start = this->_mesh_gp->panels_cnp[ib];
            const int       panel_end   = this->_mesh_gp->panels_cnp[ib + 1];
            const cusfloat* cog         = this->_input->bodies[ib]->cog;

            for ( int ip = panel_start; ip < panel_end; ip++ )
            {
                PanelGeom* panel = this->_mesh_gp->panels[ip];

                // Only submerged panels contribute
                if ( panel->center[2] >= static_cast<cusfloat>( 0.0 ) ) { continue; }

                const std::size_t sip  = static_cast<std::size_t>( ip );
                const cusfloat    area = panel->area;

                const cusfloat rx = panel->center[0] - cog[0];
                const cusfloat ry = panel->center[1] - cog[1];
                const cusfloat rz = panel->center[2] - cog[2];

                // Accumulate F = p * area * n  and  M = r × F  for one component
                auto accumulate_force = [&]( std::array<cusfloat, 6>& f, cusfloat p )
                {
                    const cusfloat dF0 = p * area * panel->normal_vec[0];
                    const cusfloat dF1 = p * area * panel->normal_vec[1];
                    const cusfloat dF2 = p * area * panel->normal_vec[2];
                    f[0] += dF0;
                    f[1] += dF1;
                    f[2] += dF2;
                    f[3] += ry * dF2 - rz * dF1;
                    f[4] += rz * dF0 - rx * dF2;
                    f[5] += rx * dF1 - ry * dF0;
                };

                accumulate_force( body_force_phi_dt[ib],       phi_dt_comp[sip]      );
                accumulate_force( body_force_kinetic[ib],      kinetic_comp[sip]     );
                accumulate_force( body_force_hydrostatic[ib],  hydrostatic_comp[sip] );
                accumulate_force( body_force_total[ib],        pressure[sip]         );
            }
        }

        this->_hdf5_exporter->append_step(
            t,
            pressure.data( ),
            phi_dt_comp.data( ),
            kinetic_comp.data( ),
            hydrostatic_comp.data( ),
            np,
            this->_body_pos,
            this->_body_vel,
            this->_body_acc,
            body_force_total,
            body_force_phi_dt,
            body_force_kinetic,
            body_force_hydrostatic
        );
    }
}


/*****************************************************************************
 * Main time loop
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::run( )
{
    std::cout << std::endl;
    std::cout << " =====================================================" << std::endl;
    std::cout << " ==          Time-domain solver running             ==" << std::endl;
    std::cout << " =====================================================" << std::endl;

    const cusfloat  sim_time    = this->_input->sim_time;
    const cusfloat  dt          = this->_input->dt;
    const int       bodies_np   = this->_input->bodies_np;
    const int       dofs_np     = this->_input->dofs_np;
    const int       n_steps     = static_cast<int>( sim_time / dt );

    std::cout << " -> Total steps: " << n_steps << "  (dt=" << dt << " s, T=" << sim_time << " s)" << std::endl;

    // Reserve sigma history with a maximum memory footprint of ~1000 steps
    // (for a proper solver this should be managed more carefully)
    const int max_hist = n_steps;
    this->_sigma_hist.reserve( max_hist );

    // Print a progress line roughly every 1 % of the total steps
    const int print_interval = std::max( 1, n_steps / 100 );

    // Prepare ParaView output directory (no-op when out_pressure is false)
    this->_init_paraview_output( );

    // Prepare HDF5 time-series file (no-op when out_pressure is false or HDF5 absent)
    this->_init_hdf5_output( );

    for ( int step=0; step<n_steps; step++ )
    {
        const cusfloat t = static_cast<cusfloat>( step ) * dt;

        // -----------------------------------------------------------------
        // 1. Update mesh positions (rigid-body kinematics)
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (1) _update_mesh_positions\n" << std::flush;
        this->_update_mesh_positions( );

        // -----------------------------------------------------------------
        // 2. Build RHS using Duhamel convolution over sigma history
        //    + body kinematic BC + incident wave diffraction BC
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (2a) _compute_body_vel_bc\n" << std::flush;
        this->_compute_body_vel_bc( 
                                        t, 
                                        this->_body_vel_bc, 
                                        this->_body_acc_bc 
                                    );

        std::cout << "[DBG] step=" << step << "  (2b) build_rhs\n" << std::flush;
        this->_kernel->build_rhs( 
                                        t, 
                                        this->_sigma_hist, 
                                        dt, 
                                        this->_body_vel_bc.data( ),
                                        this->_body_acc_bc.data( )
                                );

        // -----------------------------------------------------------------
        // 3. Solve for current source intensities sigma(t)
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (3) solve\n" << std::flush;
        this->_kernel->solve( );

        // Store current sigma in history (newest at front)
        std::cout << "[DBG] step=" << step << "  (3b) sigma_hist.insert\n" << std::flush;
        const int   np          = this->_kernel->get_n_panels( );
        const cusfloat* sigma   = this->_kernel->get_sigma( );
        this->_sigma_hist.insert(
                                    this->_sigma_hist.begin( ),
                                    std::vector<cusfloat>( sigma, sigma + np )
                                );

        // -----------------------------------------------------------------
        // 4. Compute potential derivatives dt, dx, dy and dz. They will be  
        //    used to calculate the pressure over the panels and then the 
        //    hydrodynamic forces.
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (4) compute_potential_derivatives\n" << std::flush;
        this->_kernel->compute_potential_derivatives( );

        // -----------------------------------------------------------------
        // 5. Compute hydrodynamic forces from sigma into the member vectors
        //    (_hydro_forces[ib] is referenced directly by _fext_structs[ib])
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (5) compute_hydro_forces\n" << std::flush;
        for ( int ib=0; ib<bodies_np; ib++ )
        {
            this->_compute_hydro_forces( ib, this->_hydro_forces[ib].data( ) );
        }

        // -----------------------------------------------------------------
        // 6. Advance structural dynamics (GeneralizedAlpha)
        //    The force functor (_fext_structs[ib]) already points to
        //    _hydro_forces[ib], so GA will see the updated values without
        //    needing fext_fcn to be re-assigned.
        // -----------------------------------------------------------------
        for ( int ib=0; ib<bodies_np; ib++ )
        {
            if ( this->_gen_alpha[ib] == nullptr ) { continue; }

            auto& ga = *this->_gen_alpha[ib];

            // Advance one time step
            std::cout << "[DBG] step=" << step << "  (6) ga.step body=" << ib << "\n" << std::flush;
            ga.step( );

            // Read back updated displacement / velocity / acceleration
            for ( int id=0; id<dofs_np; id++ )
            {
                this->_body_pos[ib][id] = ga.y_pos[id];
                this->_body_vel[ib][id] = ga.y_vel[id];
                this->_body_acc[ib][id] = ga.y_acc[id];
            }
        }

        // -----------------------------------------------------------------
        // 7. Write output
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (7) _output_step\n" << std::flush;
        this->_output_step( t, step );

        // -----------------------------------------------------------------
        // 8. Console progress report
        // -----------------------------------------------------------------
        if ( step == 0 || ( step + 1 ) % print_interval == 0 || step == n_steps - 1 )
        {
            const int pct = static_cast<int>( 100.0 * ( step + 1 ) / n_steps );

            std::cout << std::fixed << std::setprecision( 4 );
            std::cout << "   Step " << std::setw( 6 ) << ( step + 1 ) << "/" << n_steps
                      << "  [" << std::setw( 3 ) << pct << "%]"
                      << "  t = " << std::setw( 10 ) << t << " s";

            for ( int ib = 0; ib < bodies_np; ib++ )
            {
                std::cout << "  |  Body " << ib << ": [";
                for ( int id = 0; id < dofs_np; id++ )
                {
                    if ( id > 0 ) { std::cout << ", "; }
                    std::cout << std::setw( 12 ) << this->_body_pos[ib][id];
                }
                std::cout << " ]";
            }
            std::cout << "\n";
        }
    }

    std::cout << " -> Time-domain simulation complete." << std::endl;

    // Write the PVD collection file (no-op when out_pressure is false)
    this->_finalize_paraview_output( );

    // Flush and close the HDF5 time-series file (no-op when out_pressure is false or HDF5 absent)
    this->_finalize_hdf5_output( );
}
