
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
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
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
    // Build the incident wave once: the single description of the first-order
    // wave (amplitude, dispersion, heading and soft-start ramp) shared by the
    // diffraction BC and the kernel's Froude-Krylov potential.
    this->_inc_wave = RegularWaveFO(
                                        this->_input->wave_amp,
                                        this->_input->ang_freq,
                                        this->_input->wave_heading,
                                        this->_input->water_depth,
                                        this->_input->grav_acc,
                                        this->_input->ramp_time
                                   );

    this->_initialize_mesh_group( );
    this->_initialize_hydrostatics( );
    this->_initialize_ic_positions( );
    this->_apply_initial_displacement( );
    this->_initialize_kernel( );
    this->_compute_hydrostatic_initial_forces( );
    this->_initialize_structural_dynamics( );
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
void TimeSolver<N, NGPT>::_initialize_ic_positions( )
{
    std::cout << " -> Loading initial conditions..." << std::endl;

    const int bodies_np = this->_input->bodies_np;
    const int dofs_np   = this->_input->dofs_np;  // = 6

    this->_body_pos.resize( bodies_np );
    this->_body_vel.resize( bodies_np );
    this->_body_acc.resize( bodies_np );
    this->_struct_mass.resize( bodies_np );
    this->_gen_alpha.resize( bodies_np, nullptr );

    // Initialise per-body hydro force vectors and functor structs.
    // The structs hold a raw pointer into _hydro_forces so the complete load can
    // be updated in-place each step without touching the GA object itself.  The
    // incident-wave excitation is soft-started at its source (amplitude ramp),
    // so _hydro_forces already carries the fully-ramped total load and the
    // functor applies it without any further scaling.
    this->_hydro_forces.assign( bodies_np, std::vector<cusfloat>( dofs_np, 0.0 ) );
    this->_fext_structs.resize( bodies_np );
    for ( int ib=0; ib<bodies_np; ib++ )
    {
        this->_fext_structs[ib].forces  = this->_hydro_forces[ib].data( );
        this->_fext_structs[ib].dofs_np = dofs_np;
    }

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        const BodyDef* bd = this->_input->bodies[ib];
        for ( int k=0; k<dofs_np; k++ )
        {
            this->_body_pos[ib][k] = bd->ic_pos[k];
            this->_body_vel[ib][k] = bd->ic_vel[k];
            this->_body_acc[ib][k] = bd->ic_acc[k];
        }
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_apply_initial_displacement( )
{
    // Check whether any free body has a non-zero initial position.
    const int bodies_np = this->_input->bodies_np;
    const int dofs_np   = this->_input->dofs_np;

    bool any_displaced = false;
    for ( int ib=0; ib<bodies_np && !any_displaced; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix ) { continue; }
        for ( int k=0; k<dofs_np; k++ )
        {
            if ( this->_body_pos[ib][k] != static_cast<cusfloat>( 0.0 ) )
            {
                any_displaced = true;
                break;
            }
        }
    }

    if ( !any_displaced ) { return; }

    std::cout << " -> Applying initial displacement to mesh panels..." << std::endl;

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix ) { continue; }

        const cusfloat* pos = this->_body_pos[ib].data( );
        this->_rb_meshes[ib]->move(
                                        pos[0],  // surge  (dx)
                                        pos[1],  // sway   (dy)
                                        pos[2],  // heave  (dz)
                                        pos[3],  // roll   (drx)
                                        pos[4],  // pitch  (dry)
                                        pos[5]   // yaw    (drz)
                                  );
        this->_rb_meshes[ib]->check_underwater_panels( );
    }

    // Rebuild the MeshGroup to reflect the displaced panel geometry.
    // The BEM kernel will be constructed immediately after this call by
    // _initialize_kernel(), so only the group is rebuilt here.
    if ( this->_mesh_gp != nullptr )
    {
        delete this->_mesh_gp;
        this->_mesh_gp = nullptr;
    }

    Mesh** rb_as_mesh = new Mesh*[bodies_np];
    for ( int i=0; i<bodies_np; i++ ) { rb_as_mesh[i] = this->_rb_meshes[i]; }

    this->_mesh_gp = new MeshGroup( rb_as_mesh, bodies_np, false );
    delete [] rb_as_mesh;

    this->_mesh_gp->define_mirror_panels( );
    this->_mesh_gp->define_source_nodes( 0 );
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_hydrostatic_initial_forces( )
{
    // Populate _hydro_forces with the static load on the (possibly displaced) mesh
    // before GeneralizedAlpha::initialize() calls the force functor.
    //
    // _compute_hydro_forces() uses the linearised Bernoulli pressure:
    //   p = -rho * ( phi_dt  +  0.5*|grad(phi)|^2  +  g*z )
    // At this point all four phi arrays (_phi_dt, _phi_dx, _phi_dy, _phi_dz) are
    // zero-initialised by mkl_calloc inside _initialize_kernel(), so the dynamic
    // (radiation/diffraction) and nonlinear Bernoulli terms both vanish.
    // Only the hydrostatic term survives: p = -rho*g*z  (rho*g*|z| upward).
    // _compute_gravitational_forces() then adds the body weight so that for a body
    // floating in static equilibrium the net initial force — and therefore the
    // GA-corrected initial acceleration — is zero.

    std::cout << " -> Computing initial hydrostatic forces..." << std::endl;

    const int bodies_np = this->_input->bodies_np;

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix ) { continue; }

        this->_compute_hydro_forces( ib, this->_hydro_forces[ib].data( ) );
        this->_compute_gravitational_forces( ib, this->_hydro_forces[ib].data( ) );
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_infinite_added_mass( bool verbose )
{
    const int      bodies_np = this->_input->bodies_np;
    const int      np        = this->_kernel->get_n_panels( );
    const cusfloat rho       = this->_input->water_density;
    const cusfloat dt        = this->_input->dt;

    this->_added_mass_inf.assign( static_cast<std::size_t>( bodies_np ),
                                  std::array<cusfloat, 36>{ } );

    // Empty Duhamel history -> the unit-acceleration solves carry no memory
    // contribution, isolating the instantaneous (added-mass) reaction.
    SigmaHistType empty_hist;

    // Reusable BC buffers.  vel_bc stays zero so σ (and the velocity/gradient
    // parts of the static pressure) vanish, leaving only the σ̇-driven reaction.
    std::vector<cusfloat> vel_bc( np, static_cast<cusfloat>( 0.0 ) );
    std::vector<cusfloat> acc_bc( np, static_cast<cusfloat>( 0.0 ) );

    if ( verbose )
    {
        std::cout << "\n -> Computing infinite-frequency added mass (A_inf)..."
                  << std::endl;
    }

    for ( int ib = 0; ib < bodies_np; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix ) { continue; }

        const cusfloat*           cog = this->_input->bodies[ib]->cog;
        std::array<cusfloat, 36>& A   = this->_added_mass_inf[ static_cast<std::size_t>( ib ) ];

        for ( int jdof = 0; jdof < 6; jdof++ )
        {
            // Unit rigid-body acceleration in DOF jdof of body ib.  The
            // normal-velocity projection of a unit generalized acceleration is
            // exactly normal_vec[jdof] (6-DOF normal: [0..2] = n, [3..5] = r×n).
            std::fill( acc_bc.begin( ), acc_bc.end( ), static_cast<cusfloat>( 0.0 ) );
            for ( int j = 0; j < np; j++ )
            {
                PanelGeom* panel = this->_mesh_gp->panels[j];
                if ( panel->type != PanelTypeE::DIFFRAC ) { continue; }
                if ( this->_mesh_gp->get_body_id( j ) != ib ) { continue; }
                acc_bc[j] = panel->normal_vec[jdof];
            }

            // Solve A·σ̇ = acc_bc and assemble the static pressure field
            // (no waves, no memory).  These calls clobber the kernel's working
            // state, which is fine: the time loop rebuilds it from scratch on
            // its first step.
            //
            // build_rhs() is called for EVERY DOF — it is cheap (no
            // factorization) and, crucially, it clears the radiation/static
            // potential arrays that compute_potential_derivatives() accumulates
            // into; skipping it would let A∞ columns pile up cumulatively.
            // The system matrix is identical for all 6 DOFs (same geometry), so
            // only the FIRST DOF factorizes (solve()); the remaining five reuse
            // those LU factors via a back-substitution of the freshly assembled
            // _rhs_dt — one O(np³) factorization instead of six.  (σ stays zero
            // throughout since vel_bc = 0, so only σ̇ matters.)
            this->_kernel->build_rhs( static_cast<cusfloat>( 0.0 ),
                                      empty_hist,
                                      dt,
                                      vel_bc.data( ),
                                      acc_bc.data( ),
                                      nullptr,
                                      nullptr );
            if ( jdof == 0 ) { this->_kernel->solve( );            }
            else             { this->_kernel->backsub_sigma_dt( ); }
            this->_kernel->compute_potential_derivatives( );

            const cusfloat* phi_s = this->_kernel->get_phi_dt_rad_static( );

            // Integrate the static radiation pressure over body ib's own
            // submerged panels to get the 6-DOF reaction force.
            cusfloat force6[6];
            clear_vector( 6, force6 );
            for ( int j = 0; j < np; j++ )
            {
                PanelGeom* panel = this->_mesh_gp->panels[j];
                if ( panel->type != PanelTypeE::DIFFRAC ) { continue; }
                if ( this->_mesh_gp->get_body_id( j ) != ib ) { continue; }
                if ( panel->center[2] >= static_cast<cusfloat>( 0.0 ) ) { continue; }
                const cusfloat p_static = -rho * phi_s[j];
                this->_add_panel_pressure_force( panel, cog, p_static, force6 );
            }

            // Added-mass relation F_i = −A∞_ij · a_j with a_j = 1
            //   => A∞[i][jdof] = −F_i.
            for ( int i = 0; i < 6; i++ )
            {
                A[ static_cast<std::size_t>( i * 6 + jdof ) ] = -force6[i];
            }
        }

        // ---- Positive-definiteness / sign sanity check -------------------
        // diag_ok is always evaluated so the warning below fires even on the
        // quiet (periodic-update) path; the full 6x6 dump is verbose-only.
        bool diag_ok = true;
        for ( int i = 0; i < 6; i++ )
        {
            if ( A[ static_cast<std::size_t>( i*6 + i ) ] <= static_cast<cusfloat>( 0.0 ) )
                diag_ok = false;
        }

        if ( verbose )
        {
            std::cout << "\n   Body " << ib << "  A_inf (6x6) [row-major]:\n";
            for ( int i = 0; i < 6; i++ )
            {
                std::cout << "     ";
                for ( int j = 0; j < 6; j++ )
                {
                    std::cout << std::scientific << std::setprecision( 4 )
                              << std::setw( 13 )
                              << static_cast<double>( A[ static_cast<std::size_t>( i*6 + j ) ] )
                              << " ";
                }
                std::cout << "\n";
            }

            cusfloat max_asym = static_cast<cusfloat>( 0.0 );
            for ( int i = 0; i < 6; i++ )
                for ( int j = i + 1; j < 6; j++ )
                {
                    const cusfloat a = std::abs( A[ static_cast<std::size_t>( i*6 + j ) ]
                                               - A[ static_cast<std::size_t>( j*6 + i ) ] );
                    if ( a > max_asym ) max_asym = a;
                }

            std::cout << "   heave-heave A_inf[2][2] = "
                      << std::scientific << std::setprecision( 6 )
                      << static_cast<double>( A[ 2*6 + 2 ] )
                      << "   (must be > 0 for a physical added mass)\n";
            std::cout << "   diagonal all positive : "
                      << ( diag_ok ? "YES" : "NO  <-- sign problem!" ) << "\n";
            std::cout << "   max |A_ij - A_ji|     : "
                      << std::scientific << std::setprecision( 4 )
                      << static_cast<double>( max_asym ) << "\n";
        }

        if ( !diag_ok )
        {
            std::cout << "   [WARN] A_inf (body " << ib << ") has a non-positive "
                         "diagonal entry: the static radiation term is mis-signed and\n"
                         "          would act as NEGATIVE added mass (shorter period / "
                         "blow-up).\n";
        }
        std::cout << std::flush;
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_initialize_structural_dynamics( )
{
    std::cout << " -> Initializing structural dynamics..." << std::endl;

    // Compute A∞ on the equilibrium geometry (verbose: prints the matrix and
    // the positive-definiteness / symmetry sanity check).
    this->_compute_infinite_added_mass( true );

    const int bodies_np = this->_input->bodies_np;

    // Spectral radius at infinite frequency for the generalized-alpha scheme.
    // 1.0 = no algorithmic damping (undamped trapezoidal Newmark), which lets
    // the high-frequency mode injected by the staggered BEM<->structure coupling
    // (lagged BCs + per-step remesh) persist and slowly grow — it surfaces first
    // in the acceleration channel (the phi_dt_rad_static diagnostic, ~2x the
    // heave frequency).  A value < 1 adds high-frequency dissipation that
    // annihilates that spurious mode while leaving the physical heave resonance
    // essentially untouched.
    this->_ga_rho_inf = static_cast<cusfloat>( 0.85 );

    // Per-body backing CSR matrices, indexed by body id (nullptr for fixed
    // bodies).  Indexed storage — rather than push_back — lets a single body's
    // mass matrix be rebuilt in place when A∞ is refreshed mid-run.
    this->_csr_mass .assign( bodies_np, nullptr );
    this->_csr_stiff.assign( bodies_np, nullptr );
    this->_csr_damp .assign( bodies_np, nullptr );

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix )
        {
            // Fixed body: no structural integrator
            this->_gen_alpha[ib] = nullptr;
            continue;
        }

        // Build the integrator from the equilibrium-state initial conditions
        // (_body_pos/vel/acc already hold the input ICs) at t0 = 0.
        this->_build_body_integrator( ib, static_cast<cusfloat>( 0.0 ) );
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_build_body_integrator( int ib, cusfloat t0 )
{
    // (Re)build body ib's GeneralizedAlpha integrator with the current A∞ and
    // the current kinematics (_body_pos/vel/acc) as initial conditions.  Used
    // both at start-up and when A∞ is refreshed during the run; on a rebuild the
    // body state is carried over continuously and only the inertia operator (and
    // the integrator's factorized LHS) change.
    const int dofs_np = this->_input->dofs_np;   // = 6
    const cusfloat dt = this->_input->dt;

    // Structural mass matrix (diagonal, 6 DOF): [mass, mass, mass, Ixx, Iyy, Izz].
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

    // Cummins form: augment the structural mass with the infinite-frequency
    // added mass, (M + A∞) on the LHS.  The instantaneous (σ̇-driven) radiation
    // reaction is thereby represented implicitly here instead of as an explicit
    // RHS force — which removes the explicit added-mass instability (A∞ > M for
    // this box).  A∞ is symmetrized, (A∞ + A∞ᵀ)/2, since it is theoretically
    // symmetric and the raw operator carries ~18% asymmetry from the i==j
    // solid-angle treatment.  The matching static pressure is dropped from the
    // RHS in _compute_hydro_forces (see the `press` assembly there).
    {
        const std::array<cusfloat, 36>& Ainf =
            this->_added_mass_inf[ static_cast<std::size_t>( ib ) ];
        for ( int r = 0; r < 6; r++ )
            for ( int c = 0; c < 6; c++ )
                mass_dense[r*6 + c] += static_cast<cusfloat>( 0.5 )
                    * ( Ainf[r*6 + c] + Ainf[c*6 + r] );
    }

    // Hydrostatic stiffness and mechanical damping are both zero in the GA:
    // hydrostatic restoring is computed every step via pressure integration and
    // passed as F_ext (so K in the GA must be zero to avoid double-counting),
    // and there is no mechanical damping in this model.
    cusfloat stiff_dense[36];   clear_vector( 36, stiff_dense );
    cusfloat damp_dense[36];    clear_vector( 36, damp_dense );

    // Per-DOF kinematic restrictions (false = free, true = fixed:
    // surge, sway, heave, roll, pitch, yaw).  bool[6] -> int[6] for the GA.
    int restrictions[6];
    for ( int k=0; k<6; k++ )
        restrictions[k] = this->_input->bodies[ib]->dof_restrictions[k] ? 1 : 0;

    // Tear down any existing integrator BEFORE freeing its backing CSR: the GA's
    // MKL sparse handles alias the CSRMatrix arrays, so the order matters.
    if ( this->_gen_alpha[ib] != nullptr ) { delete this->_gen_alpha[ib]; this->_gen_alpha[ib] = nullptr; }
    if ( this->_csr_mass [ib] != nullptr ) { delete this->_csr_mass [ib]; this->_csr_mass [ib] = nullptr; }
    if ( this->_csr_stiff[ib] != nullptr ) { delete this->_csr_stiff[ib]; this->_csr_stiff[ib] = nullptr; }
    if ( this->_csr_damp [ib] != nullptr ) { delete this->_csr_damp [ib]; this->_csr_damp [ib] = nullptr; }

    this->_csr_mass [ib] = new CSRMatrix( dofs_np, mass_dense );
    this->_csr_stiff[ib] = new CSRMatrix( dofs_np, stiff_dense );
    this->_csr_damp [ib] = new CSRMatrix( dofs_np, damp_dense );

    // Initial conditions = the body's current kinematics (the input ICs at
    // start-up, or the live state on a mid-run rebuild).
    this->_gen_alpha[ib] = new GeneralizedAlpha<FextFunctor>(
                                                                &this->_fext_structs[ib],
                                                                this->_csr_mass [ib],
                                                                this->_csr_stiff[ib],
                                                                this->_csr_damp [ib],
                                                                dt,
                                                                t0,
                                                                this->_ga_rho_inf,
                                                                this->_body_pos[ib].data( ),
                                                                this->_body_vel[ib].data( ),
                                                                this->_body_acc[ib].data( ),
                                                                restrictions
                                                            );
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_update_added_mass( cusfloat t_now )
{
    // Recompute A∞ on the current (moved) wetted geometry and fold it back into
    // each free body's structural mass.  Called from the time loop every
    // ``added_mass_update_niter`` steps, after _update_mesh_positions() has
    // rebuilt the kernel for the current pose.  _compute_infinite_added_mass
    // clobbers only the kernel's transient working state, which the same step's
    // build_rhs()/solve() re-establish, so it is safe here.
    const int bodies_np = this->_input->bodies_np;

    this->_compute_infinite_added_mass( false );   // quiet recompute

    for ( int ib=0; ib<bodies_np; ib++ )
    {
        if ( this->_input->bodies[ib]->is_fix ) { continue; }
        this->_build_body_integrator( ib, t_now );   // rebuild mass + GA, carry state
    }
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_initialize_kernel( )
{
    std::cout << " -> Building steady BEM system matrix..." << std::endl;

    this->_kernel = new FormulationKernelBackendT<N, NGPT>( );
    this->_kernel->initialize( this->_mesh_gp, this->_input, this->_mpi_config, &this->_inc_wave );
    this->_resize_panel_pressure_cache( );
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_resize_panel_pressure_cache( )
{
    const std::size_t np = static_cast<std::size_t>( this->_kernel->get_n_panels( ) );
    this->_panel_phi_dt_comp            .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_phi_dt_rad_comp        .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_phi_dt_rad_static_comp .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_phi_dt_rad_memory_comp .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_phi_dt_wave_comp       .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_kinetic_comp           .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_hydrostatic_comp       .assign( np, static_cast<cusfloat>( 0.0 ) );
    this->_panel_pressure               .assign( np, static_cast<cusfloat>( 0.0 ) );
}


/*****************************************************************************
 * Mesh group rebuild (after rigid-body motion)
 ****************************************************************************/

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
    this->_kernel->initialize( this->_mesh_gp, this->_input, this->_mpi_config, &this->_inc_wave );
    this->_resize_panel_pressure_cache( );
}


/*****************************************************************************
 * Bernoulli pressure decomposition helper
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_panel_pressure_components(
    cusfloat  rho,
    cusfloat  g,
    cusfloat  phi_dt_val,
    cusfloat  phi_dx_val,
    cusfloat  phi_dy_val,
    cusfloat  phi_dz_val,
    cusfloat  z,
    cusfloat& p_phi_dt,
    cusfloat& p_kinetic,
    cusfloat& p_hydrostatic )
{
    p_phi_dt      = -rho * phi_dt_val;
    p_kinetic     = -rho * static_cast<cusfloat>( 0.5 )
                        * (   phi_dx_val * phi_dx_val
                            + phi_dy_val * phi_dy_val
                            + phi_dz_val * phi_dz_val );
    p_hydrostatic = -rho * g * z;
}


/*****************************************************************************
 * Panel pressure → 6-DOF force helper
 *****************************************************************************/


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_add_panel_pressure_force(
    PanelGeom*      panel,
    const cusfloat* cog,
    cusfloat        press,
    cusfloat*       forces )
{
    // Force: F_i = -p * area * n_i
    // (negative sign reverts the inward-pointing normal convention inherited
    //  from the frequency-domain formulation)
    const cusfloat dF0 = -press * panel->area * panel->normal_vec[0];
    const cusfloat dF1 = -press * panel->area * panel->normal_vec[1];
    const cusfloat dF2 = -press * panel->area * panel->normal_vec[2];

    forces[0] += dF0;
    forces[1] += dF1;
    forces[2] += dF2;

    // Moment: M = r × F   where r = panel_centre - CoG
    const cusfloat rx = panel->center[0] - cog[0];
    const cusfloat ry = panel->center[1] - cog[1];
    const cusfloat rz = panel->center[2] - cog[2];

    forces[3] += ry * dF2 - rz * dF1;  // roll
    forces[4] += rz * dF0 - rx * dF2;  // pitch
    forces[5] += rx * dF1 - ry * dF0;  // yaw
}


/*****************************************************************************
 * Hydrodynamic force computation
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_hydro_forces( int body_id, cusfloat* forces )
{
    // Compute the hydrodynamic pressure force on each panel using the full
    // linearised Bernoulli equation plus the hydrostatic term:
    //
    //   p = -rho * ( d(phi)/dt  +  0.5 * |grad(phi)|^2  +  g * z )
    //
    // where:
    //   d(phi)/dt  — time derivative of the velocity potential (radiation +
    //                diffraction); computed by the Duhamel convolution and
    //                stored in _phi_dt by build_rhs() / compute_potential_derivatives().
    //   0.5*|grad(phi)|^2 — nonlinear Bernoulli (kinetic) term; assembled from
    //                _phi_dx, _phi_dy, _phi_dz.  Currently disabled (multiplied
    //                by 0) pending validation of the gradient fields.
    //   g * z      — hydrostatic pressure (rho*g*|z| upward for z < 0).
    //
    // The 6-DOF force and moment contributions from panel ip are:
    //   dF_k = -p * area * n_k           (k = 0,1,2 : Fx, Fy, Fz)
    //   dM_k = -p * area * (r x n)_k     (k = 3,4,5 : Mx, My, Mz)
    // where r = panel_centre - CoG.
    // Normal signs are reversed to undo the inward-pointing convention
    // inherited from the frequency-domain formulation.

    const int dofs_np           = this->_input->dofs_np;  // 6
    const cusfloat rho          = this->_input->water_density;
    const cusfloat g            = this->_input->grav_acc;

    // Zero the base forces vector
    clear_vector( dofs_np, forces );

    // Panels owned by this body
    const int panel_start       = this->_mesh_gp->panels_cnp[body_id];
    const int panel_end         = this->_mesh_gp->panels_cnp[body_id + 1];

    // Velocity-potential derivatives at panel centres (computed by kernel)
    const cusfloat* phi_dt      = this->_kernel->get_phi_dt( );       // total (rad + wave)
    const cusfloat* phi_dx      = this->_kernel->get_phi_dx( );
    const cusfloat* phi_dy      = this->_kernel->get_phi_dy( );
    const cusfloat* phi_dz      = this->_kernel->get_phi_dz( );
    const cusfloat* phi_dt_rad        = this->_kernel->get_phi_dt_rad( );        // radiation only
    const cusfloat* phi_dt_rad_static = this->_kernel->get_phi_dt_rad_static( ); // steady σ̇·G₀ part
    const cusfloat* phi_dt_rad_memory = this->_kernel->get_phi_dt_rad_memory( ); // Duhamel convolution
    const cusfloat* phi_dt_wave       = this->_kernel->get_phi_dt_wave( );       // incident wave only

    // Centre of gravity of this body (for moment arm)
    const cusfloat* cog = this->_input->bodies[body_id]->cog;

    for ( int ip=panel_start; ip<panel_end; ip++ )
    {
        PanelGeom* panel = this->_mesh_gp->panels[ip];
        const std::size_t sip = static_cast<std::size_t>( ip );

        if ( panel->center[2] >= static_cast<cusfloat>( 0.0 ) )
        {
            // Above waterplane: zero the pressure cache slots and skip force contribution
            this->_panel_phi_dt_comp[sip]            = static_cast<cusfloat>( 0.0 );
            this->_panel_phi_dt_rad_comp[sip]        = static_cast<cusfloat>( 0.0 );
            this->_panel_phi_dt_rad_static_comp[sip] = static_cast<cusfloat>( 0.0 );
            this->_panel_phi_dt_rad_memory_comp[sip] = static_cast<cusfloat>( 0.0 );
            this->_panel_phi_dt_wave_comp[sip]       = static_cast<cusfloat>( 0.0 );
            this->_panel_kinetic_comp[sip]           = static_cast<cusfloat>( 0.0 );
            this->_panel_hydrostatic_comp[sip]       = static_cast<cusfloat>( 0.0 );
            this->_panel_pressure[sip]               = static_cast<cusfloat>( 0.0 );
            continue;
        }

        // Decompose Bernoulli pressure into its three components and cache
        cusfloat p_phi_dt, p_kinetic, p_hydrostatic;
        _compute_panel_pressure_components(
            rho, g,
            phi_dt[ip], phi_dx[ip], phi_dy[ip], phi_dz[ip],
            panel->center[2],
            p_phi_dt, p_kinetic, p_hydrostatic );

        // Radiation and incident-wave splits of the phi_dt pressure component
        const cusfloat p_phi_dt_rad        = -rho * phi_dt_rad[ip];
        const cusfloat p_phi_dt_rad_static = -rho * phi_dt_rad_static[ip];
        const cusfloat p_phi_dt_rad_memory = -rho * phi_dt_rad_memory[ip];
        const cusfloat p_phi_dt_wave       = -rho * phi_dt_wave[ip];

        this->_panel_phi_dt_comp[sip]            = p_phi_dt;
        this->_panel_phi_dt_rad_comp[sip]        = p_phi_dt_rad;
        this->_panel_phi_dt_rad_static_comp[sip] = p_phi_dt_rad_static;
        this->_panel_phi_dt_rad_memory_comp[sip] = p_phi_dt_rad_memory;
        this->_panel_phi_dt_wave_comp[sip]       = p_phi_dt_wave;
        this->_panel_kinetic_comp[sip]           = p_kinetic;
        this->_panel_hydrostatic_comp[sip]       = p_hydrostatic;
        this->_panel_pressure[sip]               = p_phi_dt + p_kinetic + p_hydrostatic;

        // Total hydrodynamic pressure on the RHS.  The static (σ̇·G₀) radiation
        // reaction is EXCLUDED — it is carried implicitly by the (M + A∞) mass
        // matrix on the LHS (Cummins form); including it here would double-count
        // the added mass and reintroduce the explicit-added-mass instability.
        // What remains is the radiation MEMORY + the diffraction + the incident
        // (Froude-Krylov) wave pressure + the hydrostatic restoring pressure.
        // The incident wave is already soft-started at the source (the amplitude
        // ramp in _compute_body_vel_bc and the kernel), so the wave excitation
        // is applied here exactly ONCE — there is no separate, post-hoc ramped
        // wave force vector.
        //   press = p_phi_dt + p_hydrostatic
        //         = −ρ·( phi_dt_rad_memory + phi_dt_diffraction + phi_dt_wave )
        //           + p_hydrostatic
        const cusfloat press = p_phi_dt + p_hydrostatic;

        // Accumulate the total hydrodynamic force contribution
        _add_panel_pressure_force( panel, cog, press, forces );
    }
}


/*****************************************************************************
 * Gravitational forces
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_gravitational_forces( int body_id, cusfloat* forces )
{
    // Gravity acts at the body CoG, so there are no moment contributions.
    // F_heave = -mass * g  (downward, opposing positive-z convention)

    const cusfloat mass = this->_input->bodies[body_id]->mass;
    const cusfloat g    = this->_input->grav_acc;

    forces[2] -= mass * g;   // heave force only; all other DOFs unchanged
}


/*****************************************************************************
 * Mesh update
 ****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_compute_body_vel_bc( 
                                                    cusfloat t, 
                                                    std::vector<cusfloat>& bc,
                                                    std::vector<cusfloat>& bc2,
                                                    std::vector<cusfloat>& bc_kin,
                                                    std::vector<cusfloat>& bc_wave
                                                )
{
    const int np = this->_kernel->get_n_panels( );
    std::cout << "  [DBG _vel_bc] entry np=" << np << "\n" << std::flush;
    bc      .assign( np, static_cast<cusfloat>( 0.0 ) );
    bc2     .assign( np, static_cast<cusfloat>( 0.0 ) );
    bc_kin  .assign( np, static_cast<cusfloat>( 0.0 ) );
    bc_wave .assign( np, static_cast<cusfloat>( 0.0 ) );

    // The incident wave (parameters, dispersion and soft-start ramp) is owned by
    // _inc_wave — the single source of truth shared with the kernel.
    const bool has_wave = this->_inc_wave.active( );

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
                bc[j]       += vel[id] * n[id];
                bc2[j]      += acc[id] * n[id];
                bc_kin[j]   += vel[id] * n[id];
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
            const cusfloat x = panel->center[0];
            const cusfloat y = panel->center[1];
            const cusfloat z = panel->center[2];

            // Diffraction condition: minus the incident-wave normal velocity
            // (and its time derivative for the acceleration BC).  Both are
            // soft-start-ramped inside _inc_wave, consistently with the
            // Froude-Krylov pressure assembled in the kernel.
            const cusfloat wave_contrib = this->_inc_wave.normal_velocity( x, y, z, n, t );
            bc[j]       -= wave_contrib;
            bc_wave[j]  -= wave_contrib;
            bc2[j]      -= this->_inc_wave.normal_acceleration( x, y, z, n, t );
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
 * Debug BEM VTU/PVD output helpers
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_init_debug_bem_output( )
{
    namespace fs = std::filesystem;

    const fs::path results_dir   = fs::path( this->_input->folder_path )
                                   / fs::path( RESULTS_FOLDER_NAME );
    const fs::path debug_bem_dir = results_dir / fs::path( RESULTS_DEBUG_BEM_FOLDER_NAME );

    if ( !fs::exists( results_dir   ) ) { fs::create_directory( results_dir   ); }
    if ( !fs::exists( debug_bem_dir ) ) { fs::create_directory( debug_bem_dir ); }

    this->_debug_bem_dir = debug_bem_dir.string( );
    this->_debug_bem_pvd_entries.clear( );

    std::cout << " -> Debug BEM output: " << this->_debug_bem_dir << std::endl;
}


template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_finalize_debug_bem_output( )
{
    if ( this->_debug_bem_pvd_entries.empty( ) ) { return; }

    namespace fs = std::filesystem;

    const std::string pvd_path = ( fs::path( this->_debug_bem_dir ) / "debug_bem.pvd" ).string( );
    if ( write_pvd( pvd_path, this->_debug_bem_pvd_entries ) )
    {
        std::cout << " -> Debug BEM PVD written: " << pvd_path << std::endl;
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
 * GA RHS debug CSV helpers
 *****************************************************************************/

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_init_ga_debug_csv( )
{
    const int bodies_np = this->_input->bodies_np;
    const int dofs_np   = this->_input->dofs_np;

    // Enable debug tracing on every active GA integrator
    for ( int ib=0; ib<bodies_np; ib++ )
    {
        if ( this->_gen_alpha[ib] != nullptr )
            this->_gen_alpha[ib]->debug_rhs = true;
    }

    // Open CSV in <folder>/1_results/
    namespace fs = std::filesystem;
    const fs::path results_dir = fs::path( this->_input->folder_path ) / "1_results";
    fs::create_directories( results_dir );
    const fs::path csv_path = results_dir / "ga_rhs_debug.csv";

    this->_ga_debug_csv.open( csv_path.string( ) );
    if ( !this->_ga_debug_csv.is_open( ) )
    {
        std::cerr << "[WARN] Could not open GA RHS debug CSV: " << csv_path << std::endl;
        return;
    }

    // Write CSV header
    this->_ga_debug_csv << "step,time,body,dof";
    for ( int id=0; id<dofs_np; id++ )
        this->_ga_debug_csv << ",fext_" << id;
    for ( int id=0; id<dofs_np; id++ )
        this->_ga_debug_csv << ",Ku_" << id;
    for ( int id=0; id<dofs_np; id++ )
        this->_ga_debug_csv << ",SDVv_" << id;
    for ( int id=0; id<dofs_np; id++ )
        this->_ga_debug_csv << ",SDAa_" << id;
    for ( int id=0; id<dofs_np; id++ )
        this->_ga_debug_csv << ",rhs_" << id;
    this->_ga_debug_csv << "\n";

    std::cout << " -> GA RHS debug CSV: " << csv_path << std::endl;
}

template<std::size_t N, int NGPT>
void TimeSolver<N, NGPT>::_finalize_ga_debug_csv( )
{
    if ( this->_ga_debug_csv.is_open( ) )
    {
        this->_ga_debug_csv.flush( );
        this->_ga_debug_csv.close( );
        std::cout << " -> GA RHS debug CSV closed." << std::endl;
    }
}


/*****************************************************************************
 * Output step
 ****************************************************************************/

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

    // ---- Build geometry arrays for VTU ----
    // Pressure component arrays (_panel_phi_dt_comp, _panel_kinetic_comp,
    // _panel_hydrostatic_comp, _panel_pressure) are already populated each step
    // by _compute_hydro_forces, so no recomputation is needed here.
    std::vector<cusfloat>   nodes_x, nodes_y, nodes_z;
    std::vector<int32_t>    connectivity, offsets;
    std::vector<uint8_t>    types;

    const std::size_t snp = static_cast<std::size_t>( np );

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
                                this->_panel_pressure.data( ),
                                this->_panel_phi_dt_comp.data( ),
                                this->_panel_phi_dt_rad_comp.data( ),
                                this->_panel_phi_dt_wave_comp.data( ),
                                this->_panel_kinetic_comp.data( ),
                                this->_panel_hydrostatic_comp.data( ),
                                this->_kernel->get_sigma( )
                            );

    // PVD entries are resolved relative to the PVD file's directory. The PVD
    // lives next to the VTU files inside the paraview dir, so the entry is
    // just the basename.
    this->_pvd_entries.emplace_back( static_cast<double>( t ), vtu_name );

    // ---- HDF5 time-series append ----
    if ( this->_hdf5_exporter != nullptr && this->_hdf5_exporter->is_open( ) )
    {
        // Integrate panel pressure components into per-body 6-DOF forces.
        const int n_bodies = this->_input->bodies_np;

        std::vector<std::array<cusfloat, 6>> body_force_total            ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_phi_dt           ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_radiation        ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_radiation_static ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_radiation_memory ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_wave             ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_kinetic          ( static_cast<std::size_t>( n_bodies ) );
        std::vector<std::array<cusfloat, 6>> body_force_hydrostatic      ( static_cast<std::size_t>( n_bodies ) );

        for ( int ib = 0; ib < n_bodies; ib++ )
        {
            body_force_total[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_phi_dt[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_radiation[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_radiation_static[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_radiation_memory[ib].fill( static_cast<cusfloat>( 0.0 ) );
            body_force_wave[ib].fill( static_cast<cusfloat>( 0.0 ) );
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

                _add_panel_pressure_force( panel, cog, this->_panel_phi_dt_comp[sip],            body_force_phi_dt[ib].data()             );
                _add_panel_pressure_force( panel, cog, this->_panel_phi_dt_rad_comp[sip],        body_force_radiation[ib].data()          );
                _add_panel_pressure_force( panel, cog, this->_panel_phi_dt_rad_static_comp[sip], body_force_radiation_static[ib].data()   );
                _add_panel_pressure_force( panel, cog, this->_panel_phi_dt_rad_memory_comp[sip], body_force_radiation_memory[ib].data()   );
                _add_panel_pressure_force( panel, cog, this->_panel_phi_dt_wave_comp[sip],       body_force_wave[ib].data()               );
                _add_panel_pressure_force( panel, cog, this->_panel_kinetic_comp[sip],           body_force_kinetic[ib].data()            );
                _add_panel_pressure_force( panel, cog, this->_panel_hydrostatic_comp[sip],       body_force_hydrostatic[ib].data()        );
                _add_panel_pressure_force( panel, cog, this->_panel_pressure[sip],               body_force_total[ib].data()              );
            }
        }

        this->_hdf5_exporter->append_step(
            t,
            this->_panel_pressure.data( ),
            this->_panel_phi_dt_comp.data( ),
            this->_panel_kinetic_comp.data( ),
            this->_panel_hydrostatic_comp.data( ),
            np,
            this->_body_pos,
            this->_body_vel,
            this->_body_acc,
            body_force_total,
            body_force_phi_dt,
            body_force_radiation,
            body_force_radiation_static,
            body_force_radiation_memory,
            body_force_wave,
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

    // Determine the maximum number of sigma history steps to retain.
    // If duhamel_hist_time > 0 the user has specified an explicit time window;
    // otherwise the full simulation history is kept (worst-case memory).
    const cusfloat hist_time  = this->_input->duhamel_hist_time;
    const int      max_hist   = ( hist_time > static_cast<cusfloat>( 0.0 ) )
                                ? std::max( 1, static_cast<int>( hist_time / dt ) )
                                : n_steps;

    std::cout << " -> Duhamel history: " << max_hist << " steps"
              << "  (" << static_cast<cusfloat>( max_hist ) * dt << " s)"
              << std::endl;

    // Reserve the circular buffer — O(1) insertions per time step, bounded memory.
    this->_sigma_hist.reserve( max_hist );

    // Print a progress line roughly every 1 % of the total steps
    const int print_interval = std::max( 1, n_steps / 100 );

    // Prepare ParaView output directory (no-op when out_pressure is false)
    this->_init_paraview_output( );

    // Prepare HDF5 time-series file (no-op when out_pressure is false or HDF5 absent)
    this->_init_hdf5_output( );

    // Prepare debug-BEM output directory (always on; the per-step VTU dump below
    // is unconditional, so its target folder must always exist).
    this->_init_debug_bem_output( );

    if constexpr ( GA_DEBUG_ON )
        this->_init_ga_debug_csv( );

    for ( int step=0; step<n_steps; step++ )
    {
        const cusfloat t = static_cast<cusfloat>( step ) * dt;

        // -----------------------------------------------------------------
        // 1. Update mesh positions (rigid-body kinematics)
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (1) _update_mesh_positions\n" << std::flush;
        this->_update_mesh_positions( );

        // -----------------------------------------------------------------
        // 1b. Periodically refresh the infinite-frequency added mass on the
        //     current (moved) geometry and re-fold it into the structural mass.
        //     Disabled when added_mass_update_niter == 0 (A∞ frozen at init).
        //     Runs right after the mesh/kernel rebuild so A∞ reflects the
        //     current pose; the same step's build_rhs/solve restore the kernel
        //     working state that _compute_infinite_added_mass clobbers.
        // -----------------------------------------------------------------
        const int amu_niter = this->_input->added_mass_update_niter;
        if ( amu_niter > 0 && step > 0 && ( step % amu_niter ) == 0 )
        {
            std::cout << "[DBG] step=" << step << "  (1b) _update_added_mass\n" << std::flush;
            this->_update_added_mass( t );
        }

        // -----------------------------------------------------------------
        // 2. Build RHS using Duhamel convolution over sigma history
        //    + body kinematic BC + incident wave diffraction BC
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (2a) _compute_body_vel_bc\n" << std::flush;
        this->_compute_body_vel_bc( 
                                        t, 
                                        this->_body_vel_bc, 
                                        this->_body_acc_bc,
                                        this->_body_kin_bc,
                                        this->_wave_bc
                                    );

        std::cout << "[DBG] step=" << step << "  (2b) build_rhs\n" << std::flush;
        this->_kernel->build_rhs( 
                                        t, 
                                        this->_sigma_hist, 
                                        dt, 
                                        this->_body_vel_bc.data( ),
                                        this->_body_acc_bc.data( ),
                                        this->_body_kin_bc.data( ),
                                        this->_wave_bc.data( )
                                );

        // -----------------------------------------------------------------
        // 3. Solve for current source intensities sigma(t)
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (3) solve\n" << std::flush;
        this->_kernel->solve( );

        // Store current sigma in history as a per-physical-face snapshot.
        // Multiple BEM panels can share a parent face (FS-refined sub-panels
        // of an FS-intersecting original face); aggregate their σ by an
        // area-weighted average so the entry represents the mean σ on the
        // wet portion of that face at this time step.
        std::cout << "[DBG] step=" << step << "  (3b) sigma_hist.push_front\n" << std::flush;
        {
            const int       np      = this->_kernel->get_n_panels( );
            const cusfloat* sigma   = this->_kernel->get_sigma( );

            SigmaFaceMap                                snapshot;
            std::unordered_map<SigmaFaceKey, cusfloat>  weight_sum;
            snapshot   .reserve( static_cast<std::size_t>( np ) );
            weight_sum .reserve( static_cast<std::size_t>( np ) );

            for ( int j=0; j<np; j++ )
            {
                PanelGeom* panel = this->_mesh_gp->panels[j];
                const SigmaFaceKey key =
                    make_face_key( panel->parent_body_id, panel->parent_face_id );
                // Fallback: a panel that somehow missed identity tagging
                // (e.g. a stationary obstacle whose mesh path bypasses the
                // remesher) is keyed by its current BEM index in a separate
                // namespace to avoid collisions.  This keeps σ_hist usable
                // even when the new identity machinery is incomplete.
                const SigmaFaceKey safe_key =
                    ( panel->parent_face_id >= 0 )
                        ? key
                        : make_face_key( -1, j );

                const cusfloat w   = ( panel->area > static_cast<cusfloat>( 0.0 ) )
                                        ? panel->area
                                        : static_cast<cusfloat>( 1.0 );
                snapshot   [safe_key] += w * sigma[j];
                weight_sum [safe_key] += w;
            }
            // Normalise the area-weighted sums into proper averages.
            for ( auto& kv : snapshot )
            {
                const cusfloat ws = weight_sum[kv.first];
                if ( ws > static_cast<cusfloat>( 0.0 ) )
                    kv.second /= ws;
            }

            this->_sigma_hist.push_front( std::move( snapshot ) );
        }

        // -----------------------------------------------------------------
        // 4. Compute potential derivatives dt, dx, dy and dz. They will be  
        //    used to calculate the pressure over the panels and then the 
        //    hydrodynamic forces.
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (4) compute_potential_derivatives\n" << std::flush;
        this->_kernel->compute_potential_derivatives( );

        // Debug: export RHS and sigma (source) fields to VTU for ParaView,
        // and record an entry for the aggregating PVD file written at run end.
        {
            std::ostringstream oss;
            oss << "debug_bem_step_" << std::setw(6) << std::setfill('0') << step << ".vtu";
            const std::string vtu_name = oss.str( );

            namespace fs = std::filesystem;
            const std::string vtu_path = ( fs::path( this->_debug_bem_dir ) / vtu_name ).string( );
            this->_kernel->export_vtu( vtu_path );

            // PVD entries are basenames (resolved relative to the PVD's own dir).
            this->_debug_bem_pvd_entries.emplace_back( static_cast<double>( t ), vtu_name );
        }

        // -----------------------------------------------------------------
        // 5. Compute hydrodynamic forces from sigma into the member vectors
        //    (_hydro_forces[ib] is referenced directly by _fext_structs[ib])
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (5) compute_hydro_forces\n" << std::flush;
        for ( int ib=0; ib<bodies_np; ib++ )
        {
            this->_compute_hydro_forces( ib, this->_hydro_forces[ib].data( ) );
            this->_compute_gravitational_forces( ib, this->_hydro_forces[ib].data( ) );
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

            if constexpr ( GA_DEBUG_ON )
                if ( this->_ga_debug_csv.is_open( ) && ga.debug_rhs )
                {
                    const cusfloat t_step = static_cast<cusfloat>( step ) * this->_input->dt;
                    this->_ga_debug_csv << std::scientific << std::setprecision( 10 )
                        << step << "," << t_step << "," << ib << "," << dofs_np;
                    for ( int id=0; id<dofs_np; id++ ) this->_ga_debug_csv << "," << ga.dbg_fext[id];
                    for ( int id=0; id<dofs_np; id++ ) this->_ga_debug_csv << "," << ga.dbg_Ku[id];
                    for ( int id=0; id<dofs_np; id++ ) this->_ga_debug_csv << "," << ga.dbg_SDVv[id];
                    for ( int id=0; id<dofs_np; id++ ) this->_ga_debug_csv << "," << ga.dbg_SDAa[id];
                    for ( int id=0; id<dofs_np; id++ ) this->_ga_debug_csv << "," << ga.rhs[id];
                    this->_ga_debug_csv << "\n";
                }
        }

        // -----------------------------------------------------------------
        // 7. Write output
        // -----------------------------------------------------------------
        std::cout << "[DBG] step=" << step << "  (7) _output_step\n" << std::flush;
        this->_output_step( t, step );

        // -----------------------------------------------------------------
        // 7b. Periodic HDF5 close+reopen to force OS-level commit.
        //     H5Fflush alone is not enough on some platforms (notably
        //     Windows, where the file size visible to other processes
        //     does not update for open handles).  Cadence is controlled
        //     by the OutputFlushNIter input field; 0 disables it.
        // -----------------------------------------------------------------
        if ( this->_input->output_flush_niter > 0
             && this->_hdf5_exporter != nullptr
             && this->_hdf5_exporter->is_open( )
             && ( ( step + 1 ) % this->_input->output_flush_niter ) == 0 )
        {
            this->_hdf5_exporter->reopen( );
        }

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

    // Write the debug-BEM PVD collection file (no-op when no steps were recorded)
    this->_finalize_debug_bem_output( );

    if constexpr ( GA_DEBUG_ON )
        this->_finalize_ga_debug_csv( );
}
