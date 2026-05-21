
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

// Include general usage libraries
#include <iostream>

// Include local modules
#include "formulation_kernel_backend_t.hpp"
#include "../../green/source.hpp"
#include "../../math/math_tools.hpp"
#include "../../tools.hpp"
#include "../../waves/wave_dispersion_base_fo.hpp"
#include "../../waves/waves_common.hpp"

// MKL include
#include "mkl.h"

// Auxiliary macros for ScaLAPACK local index mapping (mirror freq-domain style)
#define _T_COL_MAJOR_INDEX( index_cm, row_count, col_count, num_rows_local ) \
    index_cm = (col_count) * (num_rows_local) + (row_count);


template<std::size_t N, int NGPT>
FormulationKernelBackendT<N, NGPT>::~FormulationKernelBackendT( )
{
    if ( this->_solver != nullptr )
    {
        delete this->_solver;
        this->_solver = nullptr;
    }
    if ( this->_sysmat_steady != nullptr )
    {
        mkl_free( this->_sysmat_steady );
        this->_sysmat_steady = nullptr;
    }
    if ( this->_sysmat != nullptr )
    {
        mkl_free( this->_sysmat );
        this->_sysmat = nullptr;
    }
    if ( this->_rhs != nullptr )
    {
        mkl_free( this->_rhs );
        this->_rhs = nullptr;
    }
    if ( this->_rhs_dt != nullptr )
    {
        mkl_free( this->_rhs_dt );
        this->_rhs_dt = nullptr;
    }
    if ( this->_sigma != nullptr )
    {
        mkl_free( this->_sigma );
        this->_sigma = nullptr;
    }
    if ( this->_sigma_dt != nullptr )
    {
        mkl_free( this->_sigma_dt );
        this->_sigma_dt = nullptr;
    }
    if ( this->_phi_dt != nullptr )
    {
        mkl_free( this->_phi_dt );
        this->_phi_dt = nullptr;
    }
    if ( this->_phi_dx != nullptr )
    {
        mkl_free( this->_phi_dx );
        this->_phi_dx = nullptr;
    }
    if ( this->_phi_dy != nullptr )
    {
        mkl_free( this->_phi_dy );
        this->_phi_dy = nullptr;
    }
    if ( this->_phi_dz != nullptr )
    {
        mkl_free( this->_phi_dz );
        this->_phi_dz = nullptr;
    }
    if ( this->_acc_phi_dt != nullptr )
    {
        mkl_free( this->_acc_phi_dt );
        this->_acc_phi_dt = nullptr;
    }
    if ( this->_acc_phi_dx != nullptr )
    {
        mkl_free( this->_acc_phi_dx );
        this->_acc_phi_dx = nullptr;
    }
    if ( this->_acc_phi_dy != nullptr )
    {
        mkl_free( this->_acc_phi_dy );
        this->_acc_phi_dy = nullptr;
    }
    if ( this->_acc_phi_dz != nullptr )
    {
        mkl_free( this->_acc_phi_dz );
        this->_acc_phi_dz = nullptr;
    }
    if ( this->_steady_field_points != nullptr )
    {
        mkl_free( this->_steady_field_points );
        this->_steady_field_points = nullptr;
    }
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::initialize(
                                                        MeshGroup*  mesh_gp,
                                                        InputT*     input,
                                                        MpiConfig*  mpi_config
                                                    )
{
    this->_input      = input;
    this->_mesh_gp    = mesh_gp;
    this->_mpi_config = mpi_config;
    this->_n_panels   = mesh_gp->panels_tnp;

    const int np = this->_n_panels;

    // Instantiate ScaLAPACK distributed real solver (1 RHS column)
    this->_solver = new SclReal(
                                    static_cast<MKL_INT>( np ),
                                    static_cast<MKL_INT>( 1 ),
                                    static_cast<MKL_INT>( this->_mpi_config->procs_total ),
                                    static_cast<MKL_INT>( this->_mpi_config->proc_rank ),
                                    static_cast<MKL_INT>( this->_mpi_config->proc_root ),
                                    MPI_COMM_WORLD
                                );

    // Local block dimensions from ScaLAPACK distribution
    const int rows_local = this->_solver->num_rows_local;   // = np  (1 row of procs)
    const int cols_local = this->_solver->num_cols_local;

    // Allocate clean steady sysmat (never passed to Solve) + working copy
    const std::size_t sysmat_np = static_cast<std::size_t>( rows_local ) * static_cast<std::size_t>( cols_local );
    this->_sysmat_steady = static_cast<cusfloat*>( mkl_calloc( sysmat_np, sizeof(cusfloat), 64 ) );
    this->_sysmat        = static_cast<cusfloat*>( mkl_calloc( sysmat_np, sizeof(cusfloat), 64 ) );

    // RHS and solution: full vector on every process (num_rows_local = np for 1×P grid)
    this->_rhs          = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );
    this->_rhs_dt       = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );
    this->_sigma        = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );
    this->_sigma_dt     = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );

    // Potential gradient Duhamel integrals (direct eval, no additional solve)
    this->_phi_dt   = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );
    this->_phi_dx   = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );
    this->_phi_dy   = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );
    this->_phi_dz   = static_cast<cusfloat*>( mkl_calloc( np, sizeof(cusfloat), 64 ) );

    // Scratch buffers for compute_potential_derivatives() — steady contribution
    this->_acc_phi_dt          = generate_empty_vector<cusfloat>( np );
    this->_acc_phi_dx          = generate_empty_vector<cusfloat>( np );
    this->_acc_phi_dy          = generate_empty_vector<cusfloat>( np );
    this->_acc_phi_dz          = generate_empty_vector<cusfloat>( np );
    this->_steady_field_points = generate_empty_vector<cusfloat>( 3 * np );

    // Build the steady BEM matrix (fills _sysmat_steady; factorization on first solve)
    this->_build_steady_matrix( );
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::_build_steady_matrix( )
{
    const int   np          = this->_n_panels;
    const int   ndim        = 3;
    int         col_count   = 0;
    int         row_count   = 0;
    int         index_cm    = 0;

    cusfloat    field_point_i[ndim];    clear_vector( ndim, field_point_i );
    PanelGeom*  panel_i                 = nullptr;
    PanelGeom*  panel_mirror_i          = nullptr;
    cusfloat    vel_0[ndim];            clear_vector( ndim, vel_0 );
    cusfloat    vel_1[ndim];            clear_vector( ndim, vel_1 );
    cusfloat    vel_total_sf[ndim];     clear_vector( ndim, vel_total_sf );
    cusfloat    pot_0 = 0.0, pot_1 = 0.0;

    const cusfloat water_depth = this->_input->water_depth;

    // Collect all panel centres (replicated on every process)
    cusfloat* field_points = generate_empty_vector<cusfloat>( 3 * np );
    for ( int i=0; i<np; i++ )
    {
        copy_vector( 3, this->_mesh_gp->panels[i]->center, &(field_points[3*i]) );
    }

    // Outer loop: local source columns owned by this process
    col_count = 0;
    for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
    {
        panel_i        = this->_mesh_gp->panels[i];
        panel_mirror_i = this->_mesh_gp->panels_mirror[i];

        // Inner loop: all collocation rows (num_rows_local = np for 1×P grid)
        row_count = 0;
        for ( int j=0; j<np; j++ )
        {
            clear_vector( ndim, vel_0 );
            clear_vector( ndim, vel_1 );
            clear_vector( ndim, vel_total_sf );
            pot_0 = 0.0; pot_1 = 0.0;

            // r0 source
            calculate_source_newman( panel_i, &(field_points[3*j]), 0, 0, vel_0, pot_0 );

            // r1 source mirror
            field_point_i[0] = field_points[3*j  ];
            field_point_i[1] = field_points[3*j+1];
            field_point_i[2] = field_points[3*j+2];
            calculate_source_newman( panel_mirror_i, field_point_i, 0, 0, vel_1, pot_1 );

            // Compose total normal velocity (regular finite-depth)
            vel_total_sf[0] = - ( vel_0[0] - vel_1[0] );
            vel_total_sf[1] = - ( vel_0[1] - vel_1[1] );
            vel_total_sf[2] = - ( vel_0[2] - vel_1[2] );

            // Normal projection (influence of source panel i on collocation j)
            cusfloat int_dn;
            if ( i == j )
            {
                int_dn = static_cast<cusfloat>( 0.5 );
            }
            else
            {
                int_dn = (
                              this->_mesh_gp->source_nodes[j]->normal_vec[0] * vel_total_sf[0]
                            + this->_mesh_gp->source_nodes[j]->normal_vec[1] * vel_total_sf[1]
                            + this->_mesh_gp->source_nodes[j]->normal_vec[2] * vel_total_sf[2]
                          ) / static_cast<cusfloat>( 4.0 * PI );
            }

            // Store in column-major local block
            _T_COL_MAJOR_INDEX( index_cm, row_count, col_count, this->_solver->num_rows_local )
            this->_sysmat_steady[index_cm] = int_dn;

            row_count++;
        }

        col_count++;
    }

    mkl_free( field_points );

    MPI_Barrier( MPI_COMM_WORLD );
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::build_rhs(
                                                        cusfloat                                            t_current,
                                                        const std::vector<std::vector<cusfloat>>&           sigma_hist,
                                                        cusfloat                                            dt,
                                                        const cusfloat*                                     body_vel_bc,
                                                        const cusfloat*                                     body_acc_bc
                                                   )
{
    const int np     = this->_n_panels;
    const int n_hist = static_cast<int>( sigma_hist.size( ) );

    // Cache current time for use by compute_potential_derivatives()
    this->_t_current = t_current;

    // Initialize RHS to zero (partial contributions from this process)
    clear_vector( np, this->_rhs        );
    clear_vector( np, this->_rhs_dt     );
    clear_vector( np, this->_phi_dt );
    clear_vector( np, this->_phi_dx );
    clear_vector( np, this->_phi_dy );
    clear_vector( np, this->_phi_dz );

    // -----------------------------------------------------------------------
    // Body kinematic BC: rhs[j] += n_j · vel_body
    // body_vel_bc[j] = normal velocity of the body surface at collocation j
    // -----------------------------------------------------------------------
    if ( body_vel_bc != nullptr )
    {
        for ( int j=0; j<np; j++ )
        {
            this->_rhs[j]       += body_vel_bc[j];
            this->_rhs_dt[j]    += body_acc_bc[j];
        }
    }

    // -----------------------------------------------------------------------
    // Duhamel convolution: subtract influence of past source intensities
    // Each process handles its local source columns [start_col_0, end_col_0)
    // Contributions from all processes are summed via MPI_Allreduce.
    // -----------------------------------------------------------------------
    for ( int k=0; k<n_hist; k++ )
    {
        const cusfloat t_lag_end   = t_current - static_cast<cusfloat>( k     ) * dt;
        const cusfloat t_lag_start = t_current - static_cast<cusfloat>( k + 1 ) * dt;

        if ( t_lag_end <= t_lag_start ) { continue; }

        const std::vector<cusfloat>& sigma_k = sigma_hist[k];

        // Outer: local source panels (columns)
        for ( int src=this->_solver->start_col_0; src<this->_solver->end_col_0; src++ )
        {
            SourceNode* source_src = this->_mesh_gp->source_nodes[src];

            // Inner: all observation rows
            for ( int obs=0; obs<np; obs++ )
            {
                SourceNode* source_obs = this->_mesh_gp->source_nodes[obs];

                // Set up time-domain Green's function integrator:
                //   source_i = integration panel (src)
                //   source_j = position source   (obs)
                this->_gwtfcns_interf.set_source_i( source_src, static_cast<cusfloat>( 1.0 ) );
                this->_gwtfcns_interf.set_source_j( source_obs );

                cusfloat dtn_val   = 0.0;
                cusfloat dtx_val   = 0.0;
                cusfloat dty_val   = 0.0;
                cusfloat dtz_val   = 0.0;
                cusfloat dtt_val   = 0.0;
                cusfloat dttn_val  = 0.0;
                cusfloat dttx_val  = 0.0;
                cusfloat dtty_val  = 0.0;
                cusfloat dttz_val  = 0.0;

                quadrature_panel_time_t<PanelGeom, GWTFcnsInterfaceT<N*N>, N, NGPT>(
                                                                                        source_src->panel,
                                                                                        this->_gwtfcns_interf,
                                                                                        t_lag_start,
                                                                                        t_lag_end,
                                                                                        dtn_val,
                                                                                        dtx_val,
                                                                                        dty_val,
                                                                                        dtz_val,
                                                                                        dtt_val,
                                                                                        dttn_val,
                                                                                        dttx_val,
                                                                                        dtty_val,
                                                                                        dttz_val,
                                                                                        false
                                                                                     );

                // Accumulate Duhamel contribution (partial sum for local columns)
                this->_rhs[obs]         -= sigma_k[src] * dtn_val  / static_cast<cusfloat>( 4.0 * PI );
                this->_rhs_dt[obs]      -= sigma_k[src] * dttn_val / static_cast<cusfloat>( 4.0 * PI );
                this->_phi_dt[obs]  -= sigma_k[src] * dtt_val / static_cast<cusfloat>( 4.0 * PI );
                this->_phi_dx[obs]  -= sigma_k[src] * dtx_val / static_cast<cusfloat>( 4.0 * PI );
                this->_phi_dy[obs]  -= sigma_k[src] * dty_val / static_cast<cusfloat>( 4.0 * PI );
                this->_phi_dz[obs]  -= sigma_k[src] * dtz_val / static_cast<cusfloat>( 4.0 * PI );
            }
        }
    }

    // Sum partial RHS contributions from all processes
    MPI_Allreduce(
                    MPI_IN_PLACE,
                    this->_rhs,
                    np,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                 );

    MPI_Allreduce(
                    MPI_IN_PLACE,
                    this->_rhs_dt,
                    np,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                 );

    MPI_Allreduce(
                    MPI_IN_PLACE,
                    this->_phi_dt,
                    np,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                 );

    MPI_Allreduce(
                    MPI_IN_PLACE,
                    this->_phi_dx,
                    np,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                 );

    MPI_Allreduce(
                    MPI_IN_PLACE,
                    this->_phi_dy,
                    np,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                 );

    MPI_Allreduce(
                    MPI_IN_PLACE,
                    this->_phi_dz,
                    np,
                    mpi_cusfloat,
                    MPI_SUM,
                    MPI_COMM_WORLD
                 );
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::solve( )
{
    const int np = this->_n_panels;

    // Copy clean steady sysmat into working buffer (pgesv overwrites it with LU factors)
    const std::size_t sysmat_np = static_cast<std::size_t>( this->_solver->num_rows_local )
                                * static_cast<std::size_t>( this->_solver->num_cols_local );
    copy_vector( static_cast<int>( sysmat_np ), this->_sysmat_steady, this->_sysmat );

    // Copy rhs into sigma (pgesv overwrites it with the solution).
    // After this call _sysmat holds the LU factorization and _solver->ipiv
    // holds the pivot vector — both are reused for the sigma_dt solve below.
    copy_vector( np, this->_rhs, this->_sigma );

    // Distributed solve: A * sigma = rhs  (ScaLAPACK pgesv: factorize + back-sub)
    this->_solver->Solve( this->_sysmat, this->_sigma );

    // Broadcast sigma to all processes for consistency
    MPI_Bcast(
                this->_sigma,
                np,
                mpi_cusfloat,
                this->_mpi_config->proc_root,
                MPI_COMM_WORLD
             );

    // -----------------------------------------------------------------
    // Solve for sigma_dt: A * sigma_dt = rhs_dt
    // The LU factorization of A is already stored in _sysmat from the
    // pgesv call above, so we use pgetrs (back-substitution only) to
    // avoid a second O(n³) factorization.
    // -----------------------------------------------------------------
    copy_vector( np, this->_rhs_dt, this->_sigma_dt );

    this->_solver->SolveWithStoredLU( this->_sysmat, this->_sigma_dt );

    MPI_Bcast(
                this->_sigma_dt,
                np,
                mpi_cusfloat,
                this->_mpi_config->proc_root,
                MPI_COMM_WORLD
             );
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::compute_potential_derivatives( )
{
    // _phi_dt/dx/dy/dz already contain the radiation (Duhamel) contribution
    // assembled in build_rhs().  Here we add:
    //   (a) the incident wave contribution  (phi_wave)
    //   (b) the steady Green's function contribution  (phi_steady)
    // so that the total velocity potential derivatives
    //   phi_d* = phi_rad_d* + phi_wave_d* + phi_steady_d*
    // are available at each collocation point for pressure and force evaluation.

    const int np = this->_n_panels;

    // ----------------------------------------------------------------
    // (a) Incident wave contributions (skipped when no wave is active)
    // ----------------------------------------------------------------
    {
        const cusfloat w  = this->_input->ang_freq;
        const cusfloat aw = this->_input->wave_amp;
        const cusfloat h  = this->_input->water_depth;
        const cusfloat g  = this->_input->grav_acc;
        const cusfloat mu = this->_input->head * PI / static_cast<cusfloat>( 180.0 );
        const cusfloat t  = this->_t_current;

        if ( aw > static_cast<cusfloat>( 0.0 ) && w > static_cast<cusfloat>( 0.0 ) )
        {
            const cusfloat k = w2k( w, h, g );
            for ( int j=0; j<np; j++ )
            {
                const cusfloat x = this->_mesh_gp->panels[j]->center[0];
                const cusfloat y = this->_mesh_gp->panels[j]->center[1];
                const cusfloat z = this->_mesh_gp->panels[j]->center[2];

                this->_phi_dt[j] += wave_potential_fo_time_dt( aw, w, k, h, g, x, y, z, mu, t );
                this->_phi_dx[j] += wave_potential_fo_time_dx( aw, w, k, h, g, x, y, z, mu, t );
                this->_phi_dy[j] += wave_potential_fo_time_dy( aw, w, k, h, g, x, y, z, mu, t );
                this->_phi_dz[j] += wave_potential_fo_time_dz( aw, w, k, h, g, x, y, z, mu, t );
            }
        }
    }

    // ----------------------------------------------------------------
    // (b) Steady Green's function contributions
    //
    // For each collocation point j, accumulate the matrix-vector product:
    //   phi_dx[j] += Σ_i  sigma[i] * vel_total_sf[0](j,i) / 4π
    //   phi_dy[j] += Σ_i  sigma[i] * vel_total_sf[1](j,i) / 4π
    //   phi_dz[j] += Σ_i  sigma[i] * vel_total_sf[2](j,i) / 4π
    //   phi_dt[j] += Σ_i  sigma[i] * ( pot_0 - pot_1 )(j,i) / 4π
    //
    // The outer loop runs over the local ScaLAPACK column slice; partial sums
    // are reduced across all MPI processes via MPI_Allreduce before being added
    // to the output vectors.
    //
    // Velocity and potential kernels mirror _build_steady_matrix:
    //   vel_total_sf[k] = -(vel_0[k] - vel_1[k])   (direct - mirror-image source)
    //   pot_total       =   pot_0 - pot_1
    // Diagonal term (i==j): velocity uses the solid-angle free-surface condition
    //   vel_total_sf[k] = 0.5 * normal_vec[k];
    //   pot_total = 0 (singular; no contribution).
    // ----------------------------------------------------------------
    {
        const int   ndim = 3;
        cusfloat    fp_i[ndim];
        cusfloat    vel_0[ndim], vel_1[ndim], vel_sf[ndim];
        cusfloat    pot_0, pot_1;

        // Clear preallocated partial-sum accumulators before use
        clear_vector( np, this->_acc_phi_dt );
        clear_vector( np, this->_acc_phi_dx );
        clear_vector( np, this->_acc_phi_dy );
        clear_vector( np, this->_acc_phi_dz );

        // Refresh panel centres (mesh may have moved since last step)
        for ( int i=0; i<np; i++ )
        {
            copy_vector( 3, this->_mesh_gp->panels[i]->center, &(this->_steady_field_points[3*i]) );
        }

        // Outer loop: local source columns owned by this MPI process
        for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
        {
            PanelGeom* panel_i_s   = this->_mesh_gp->panels[i];
            PanelGeom* panel_mir_i = this->_mesh_gp->panels_mirror[i];
            const cusfloat sigma_i = this->_sigma[i];

            // Inner loop: all collocation rows
            for ( int j=0; j<np; j++ )
            {
                clear_vector( ndim, vel_0 );
                clear_vector( ndim, vel_1 );
                clear_vector( ndim, vel_sf );
                pot_0 = 0.0;  pot_1 = 0.0;

                // r0: direct Rankine source evaluated at collocation point j
                calculate_source_newman( panel_i_s, &(this->_steady_field_points[3*j]), 0, 0, vel_0, pot_0 );

                // r1: mirror-image source at the same field point
                fp_i[0] = this->_steady_field_points[3*j  ];
                fp_i[1] = this->_steady_field_points[3*j+1];
                fp_i[2] = this->_steady_field_points[3*j+2];
                calculate_source_newman( panel_mir_i, fp_i, 0, 0, vel_1, pot_1 );

                // Compose total velocity components (mirrors _build_steady_matrix)
                vel_sf[0] = - ( vel_0[0] - vel_1[0] );
                vel_sf[1] = - ( vel_0[1] - vel_1[1] );
                vel_sf[2] = - ( vel_0[2] - vel_1[2] );

                if ( i == j )
                {
                    // Diagonal: free-surface solid-angle condition
                    const cusfloat half = static_cast<cusfloat>( 0.5 );
                    this->_acc_phi_dx[j] += sigma_i * half * panel_i_s->normal_vec[0];
                    this->_acc_phi_dy[j] += sigma_i * half * panel_i_s->normal_vec[1];
                    this->_acc_phi_dz[j] += sigma_i * half * panel_i_s->normal_vec[2];
                    // pot_total is singular on the diagonal: no phi_dt contribution
                }
                else
                {
                    const cusfloat inv4pi = static_cast<cusfloat>( 1.0 )
                                         / static_cast<cusfloat>( 4.0 * PI );
                    this->_acc_phi_dx[j] += sigma_i * vel_sf[0] * inv4pi;
                    this->_acc_phi_dy[j] += sigma_i * vel_sf[1] * inv4pi;
                    this->_acc_phi_dz[j] += sigma_i * vel_sf[2] * inv4pi;
                    this->_acc_phi_dt[j] += sigma_i * ( pot_0 - pot_1 ) * inv4pi;
                }
            }
        }

        // Reduce partial sums across all MPI processes
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dt, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dx, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dy, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dz, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );

        // Accumulate into the output vectors
        for ( int j=0; j<np; j++ )
        {
            this->_phi_dt[j] += this->_acc_phi_dt[j];
            this->_phi_dx[j] += this->_acc_phi_dx[j];
            this->_phi_dy[j] += this->_acc_phi_dy[j];
            this->_phi_dz[j] += this->_acc_phi_dz[j];
        }
    }
}
