
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
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

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
    if ( this->_rhs_body_kin != nullptr )
    {
        mkl_free( this->_rhs_body_kin );
        this->_rhs_body_kin = nullptr;
    }
    if ( this->_rhs_duhamel != nullptr )
    {
        mkl_free( this->_rhs_duhamel );
        this->_rhs_duhamel = nullptr;
    }
    if ( this->_rhs_wave != nullptr )
    {
        mkl_free( this->_rhs_wave );
        this->_rhs_wave = nullptr;
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
    // Radiation split arrays
    if ( this->_phi_dt_rad != nullptr )
    {
        mkl_free( this->_phi_dt_rad );
        this->_phi_dt_rad = nullptr;
    }
    if ( this->_phi_dx_rad != nullptr )
    {
        mkl_free( this->_phi_dx_rad );
        this->_phi_dx_rad = nullptr;
    }
    if ( this->_phi_dy_rad != nullptr )
    {
        mkl_free( this->_phi_dy_rad );
        this->_phi_dy_rad = nullptr;
    }
    if ( this->_phi_dz_rad != nullptr )
    {
        mkl_free( this->_phi_dz_rad );
        this->_phi_dz_rad = nullptr;
    }
    // Incident-wave split arrays
    if ( this->_phi_dt_wave != nullptr )
    {
        mkl_free( this->_phi_dt_wave );
        this->_phi_dt_wave = nullptr;
    }
    if ( this->_phi_dx_wave != nullptr )
    {
        mkl_free( this->_phi_dx_wave );
        this->_phi_dx_wave = nullptr;
    }
    if ( this->_phi_dy_wave != nullptr )
    {
        mkl_free( this->_phi_dy_wave );
        this->_phi_dy_wave = nullptr;
    }
    if ( this->_phi_dz_wave != nullptr )
    {
        mkl_free( this->_phi_dz_wave );
        this->_phi_dz_wave = nullptr;
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
    const int rows_local        = this->_solver->num_rows_local;   // = np  (1 row of procs)
    const int cols_local        = this->_solver->num_cols_local;

    // Allocate clean steady sysmat (never passed to Solve) + working copy
    const std::size_t sysmat_np = static_cast<std::size_t>( rows_local ) * static_cast<std::size_t>( cols_local );
    this->_sysmat_steady        = generate_empty_vector<cusfloat>( sysmat_np );
    this->_sysmat               = generate_empty_vector<cusfloat>( sysmat_np );

    // RHS and solution: full vector on every process (num_rows_local = np for 1×P grid)
    this->_rhs                  = generate_empty_vector<cusfloat>( np );
    this->_rhs_dt               = generate_empty_vector<cusfloat>( np );
    this->_sigma                = generate_empty_vector<cusfloat>( np );
    this->_sigma_dt             = generate_empty_vector<cusfloat>( np );

    // Decomposed RHS contributions (debug export)
    this->_rhs_body_kin         = generate_empty_vector<cusfloat>( np );
    this->_rhs_duhamel          = generate_empty_vector<cusfloat>( np );
    this->_rhs_wave             = generate_empty_vector<cusfloat>( np );

    // Potential gradient: total (= rad + wave), updated at end of compute_potential_derivatives
    this->_phi_dt               = generate_empty_vector<cusfloat>( np );
    this->_phi_dx               = generate_empty_vector<cusfloat>( np );
    this->_phi_dy               = generate_empty_vector<cusfloat>( np );
    this->_phi_dz               = generate_empty_vector<cusfloat>( np );

    // Split arrays: radiation (Duhamel + steady Rankine)
    this->_phi_dt_rad           = generate_empty_vector<cusfloat>( np );
    this->_phi_dx_rad           = generate_empty_vector<cusfloat>( np );
    this->_phi_dy_rad           = generate_empty_vector<cusfloat>( np );
    this->_phi_dz_rad           = generate_empty_vector<cusfloat>( np );

    // Split arrays: incident wave
    this->_phi_dt_wave          = generate_empty_vector<cusfloat>( np );
    this->_phi_dx_wave          = generate_empty_vector<cusfloat>( np );
    this->_phi_dy_wave          = generate_empty_vector<cusfloat>( np );
    this->_phi_dz_wave          = generate_empty_vector<cusfloat>( np );

    // Scratch buffers for compute_potential_derivatives() — steady contribution
    this->_acc_phi_dt           = generate_empty_vector<cusfloat>( np );
    this->_acc_phi_dx           = generate_empty_vector<cusfloat>( np );
    this->_acc_phi_dy           = generate_empty_vector<cusfloat>( np );
    this->_acc_phi_dz           = generate_empty_vector<cusfloat>( np );
    this->_steady_field_points  = generate_empty_vector<cusfloat>( 3 * np );

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
                                                        const CircularBuffer<std::vector<cusfloat>>&        sigma_hist,
                                                        cusfloat                                            dt,
                                                        const cusfloat*                                     body_vel_bc,
                                                        const cusfloat*                                     body_acc_bc,
                                                        const cusfloat*                                     body_kin_bc,
                                                        const cusfloat*                                     wave_bc
                                                   )
{
    const int np     = this->_n_panels;
    const int n_hist = static_cast<int>( sigma_hist.size( ) );

    std::cout << "  [DBG build_rhs] entry np=" << np << " n_hist=" << n_hist << "\n" << std::flush;

    // Cache current time for use by compute_potential_derivatives()
    this->_t_current = t_current;

    // Initialize RHS to zero (partial contributions from this process)
    clear_vector( np, this->_rhs        );
    clear_vector( np, this->_rhs_dt     );
    // Clear decomposed contribution arrays
    clear_vector( np, this->_rhs_body_kin );
    clear_vector( np, this->_rhs_duhamel  );
    clear_vector( np, this->_rhs_wave     );
    // Clear radiation-split arrays; wave and total arrays are set in compute_potential_derivatives
    clear_vector( np, this->_phi_dt_rad );
    clear_vector( np, this->_phi_dx_rad );
    clear_vector( np, this->_phi_dy_rad );
    clear_vector( np, this->_phi_dz_rad );

    // -----------------------------------------------------------------------
    // Body kinematic BC: rhs[j] += n_j · vel_body
    // body_vel_bc[j] = combined normal velocity (body motion + wave diffraction)
    // body_kin_bc[j] = body motion contribution only  (optional, for debug export)
    // wave_bc[j]     = wave diffraction contribution only (optional, for debug export)
    // -----------------------------------------------------------------------
    if ( body_vel_bc != nullptr )
    {
        for ( int j=0; j<np; j++ )
        {
            this->_rhs[j]       += body_vel_bc[j];
            this->_rhs_dt[j]    += body_acc_bc[j];
        }
    }
    // Store individual contributions if provided
    if ( body_kin_bc != nullptr )
    {
        for ( int j=0; j<np; j++ )
        {
            this->_rhs_body_kin[j] = body_kin_bc[j];
        }
    }
    if ( wave_bc != nullptr )
    {
        for ( int j=0; j<np; j++ )
        {
            this->_rhs_wave[j] = wave_bc[j];
        }
    }

    // -----------------------------------------------------------------------
    // Duhamel convolution: subtract influence of past source intensities
    // Each process handles its local source columns [start_col_0, end_col_0)
    // Contributions from all processes are summed via MPI_Allreduce.
    // -----------------------------------------------------------------------
    if ( this->_input->use_duhamel )
    {
        std::cout << "  [DBG build_rhs] Duhamel loop start ("
                  << ( this->_input->use_gk_duhamel ? "GK+sigma-interp" : "trapezoidal" )
                  << ")\n" << std::flush;
        const auto _t_duhamel_start = std::chrono::high_resolution_clock::now();

        if ( this->_input->use_gk_duhamel )
        {
            this->_accumulate_duhamel_gk_sigma( sigma_hist, t_current, dt );
        }
        else
        {
            this->_accumulate_duhamel_trapezoidal( sigma_hist, t_current, dt );
        }

        const double _t_duhamel_s = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - _t_duhamel_start ).count();

        this->_gwtfcns_interf.print_timers();

        const auto _t_mpi_start = std::chrono::high_resolution_clock::now();
        std::cout << "  [DBG build_rhs] Duhamel loop done  (" << _t_duhamel_s << " s), MPI_Allreduce start\n" << std::flush;
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

        // Also reduce the Duhamel-only contribution (distributed across columns)
        MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_rhs_duhamel,
                        np,
                        mpi_cusfloat,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

        MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_phi_dt_rad,
                        np,
                        mpi_cusfloat,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

        MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_phi_dx_rad,
                        np,
                        mpi_cusfloat,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

        MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_phi_dy_rad,
                        np,
                        mpi_cusfloat,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

        MPI_Allreduce(
                        MPI_IN_PLACE,
                        this->_phi_dz_rad,
                        np,
                        mpi_cusfloat,
                        MPI_SUM,
                        MPI_COMM_WORLD
                    );

        const double _t_mpi_s = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - _t_mpi_start ).count();
        std::cout << "  [DBG build_rhs] MPI_Allreduce done (" << _t_mpi_s << " s)\n" << std::flush;

    } // end if ( use_duhamel )

    std::cout << "  [DBG build_rhs] done\n" << std::flush;
}


// =============================================================================
// Duhamel accumulator — trapezoidal (history-step constant-sigma) rule
// =============================================================================

template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::_accumulate_duhamel_trapezoidal(
                                                                            const CircularBuffer<std::vector<cusfloat>>&  sigma_hist,
                                                                            cusfloat                                      t_current,
                                                                            cusfloat                                      dt
                                                                         )
{
    const int np     = this->_n_panels;
    const int n_hist = static_cast<int>( sigma_hist.size( ) );

    for ( int k=0; k<n_hist; k++ )
    {
        // Integrate the kernel in the LAG variable s = t_current - tau.
        // Segment k corresponds to lag in [k*dt, (k+1)*dt] (k=0 is most recent).
        // set_time_diff() receives the lag, which is what the free-surface
        // kernel expects (beta = sqrt(g/R) * lag).
        const cusfloat t_lag_start = static_cast<cusfloat>( k     ) * dt;
        const cusfloat t_lag_end   = static_cast<cusfloat>( k + 1 ) * dt;

        if ( t_lag_end <= t_lag_start ) { continue; }

        const std::vector<cusfloat>& sigma_k = sigma_hist[k];

        // Outer loop: local source panels (columns owned by this MPI process)
        for ( int src=this->_solver->start_col_0; src<this->_solver->end_col_0; src++ )
        {
            SourceNode* source_src = this->_mesh_gp->source_nodes[src];
            const cusfloat sigma_src = sigma_k[src];   // constant over [t_lag_start, t_lag_end]

            // Inner loop: all observation / collocation rows
            for ( int obs=0; obs<np; obs++ )
            {
                SourceNode* source_obs = this->_mesh_gp->source_nodes[obs];

                this->_gwtfcns_interf.set_source_i( source_src, static_cast<cusfloat>( 1.0 ) );
                this->_gwtfcns_interf.set_source_j( source_obs );

                cusfloat dtn_val  = 0.0;
                cusfloat dtx_val  = 0.0;
                cusfloat dty_val  = 0.0;
                cusfloat dtz_val  = 0.0;
                cusfloat dtt_val  = 0.0;
                cusfloat dttn_val = 0.0;
                cusfloat dttx_val = 0.0;
                cusfloat dtty_val = 0.0;
                cusfloat dttz_val = 0.0;

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

                const cusfloat inv4pi = static_cast<cusfloat>( 1.0 / ( 4.0 * PI ) );
                const cusfloat contrib_rhs = sigma_src * dtn_val * inv4pi;
                this->_rhs[obs]         += contrib_rhs;
                this->_rhs_duhamel[obs] += contrib_rhs;
                this->_rhs_dt[obs]      += sigma_src * dttn_val * inv4pi;
                this->_phi_dt_rad[obs]  += sigma_src * dtt_val  * inv4pi;
                this->_phi_dx_rad[obs]  += sigma_src * dtx_val  * inv4pi;
                this->_phi_dy_rad[obs]  += sigma_src * dty_val  * inv4pi;
                this->_phi_dz_rad[obs]  += sigma_src * dtz_val  * inv4pi;
            }
        }
    }
}


// =============================================================================
// Duhamel accumulator — GK + piecewise-linear sigma interpolation
// =============================================================================

template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::_accumulate_duhamel_gk_sigma(
                                                                          const CircularBuffer<std::vector<cusfloat>>&  sigma_hist,
                                                                          cusfloat                                      t_current,
                                                                          cusfloat                                      dt
                                                                       )
{
    const int np     = this->_n_panels;
    const int n_hist = static_cast<int>( sigma_hist.size( ) );
    const int n_seg  = std::max( 1, this->_input->gk_n_seg );  // composite GK macro-intervals

    if ( n_hist == 0 ) { return; }

    // Divide the full Duhamel history window into n_seg equal macro-intervals.
    // Each macro-interval [t_seg_start, t_seg_end] is integrated with a single
    // G7K15 Gauss-Kronrod call (15 points).  sigma(t) at any quadrature point
    // is obtained by piecewise-linear interpolation from the discrete history.
    const cusfloat T_hist = static_cast<cusfloat>( n_hist ) * dt;   // total history span [s]
    const cusfloat T_seg  = T_hist / static_cast<cusfloat>( n_seg ); // macro-interval width [s]

    for ( int seg = 0; seg < n_seg; seg++ )
    {
        // Integrate the kernel in the LAG variable s = t_current - tau.
        // Segment seg covers lag in [seg*T_seg, (seg+1)*T_seg] (seg=0 is most recent).
        // set_time_diff() receives the lag inside the quadrature, matching the
        // kernel convention (beta = sqrt(g/R) * lag).
        const cusfloat t_seg_start = static_cast<cusfloat>( seg     ) * T_seg;
        const cusfloat t_seg_end   = static_cast<cusfloat>( seg + 1 ) * T_seg;

        if ( t_seg_end <= t_seg_start ) { continue; }

        for ( int src = this->_solver->start_col_0; src < this->_solver->end_col_0; src++ )
        {
            SourceNode* source_src = this->_mesh_gp->source_nodes[src];

            // sigma at lag s: piecewise-linear interpolation from the discrete history.
            //   k_exact = s / dt              (fractional history index)
            //   sigma_hist[k_lo] is sigma at t_current - k_lo*dt  (newer side)
            //   sigma_hist[k_hi] is sigma at t_current - k_hi*dt  (older side)
            const auto sigma_at = [&sigma_hist, n_hist, dt, src]( cusfloat t_lag ) -> cusfloat
            {
                const cusfloat k_real = t_lag / dt;
                const int      k_lo   = static_cast<int>( std::floor( k_real ) );
                const cusfloat alpha  = k_real - static_cast<cusfloat>( k_lo );
                const int      k_lo_c = std::max( 0, std::min( k_lo,     n_hist - 1 ) );
                const int      k_hi_c = std::max( 0, std::min( k_lo + 1, n_hist - 1 ) );
                return ( static_cast<cusfloat>( 1.0 ) - alpha ) * sigma_hist[k_lo_c][src]
                     +   alpha                                  * sigma_hist[k_hi_c][src];
            };

            for ( int obs = 0; obs < np; obs++ )
            {
                SourceNode* source_obs = this->_mesh_gp->source_nodes[obs];

                this->_gwtfcns_interf.set_source_i( source_src, static_cast<cusfloat>( 1.0 ) );
                this->_gwtfcns_interf.set_source_j( source_obs );

                cusfloat dtn_val  = 0.0, dtx_val  = 0.0, dty_val  = 0.0, dtz_val  = 0.0;
                cusfloat dtt_val  = 0.0, dttn_val = 0.0, dttx_val = 0.0;
                cusfloat dtty_val = 0.0, dttz_val = 0.0;

                // sigma already folded into the quadrature weights;
                // outputs = ∫ sigma(t) · G_kernel(t) dt  over [t_seg_start, t_seg_end].
                quadrature_panel_time_t_sigma_gl<PanelGeom, GWTFcnsInterfaceT<N*N>, N, 15>(
                                                                                              source_src->panel,
                                                                                              this->_gwtfcns_interf,
                                                                                              t_seg_start,
                                                                                              t_seg_end,
                                                                                              sigma_at,
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

                const cusfloat inv4pi      = static_cast<cusfloat>( 1.0 / ( 4.0 * PI ) );
                const cusfloat contrib_rhs = dtn_val * inv4pi;
                this->_rhs[obs]         += contrib_rhs;
                this->_rhs_duhamel[obs] += contrib_rhs;
                this->_rhs_dt[obs]      += dttn_val * inv4pi;
                this->_phi_dt_rad[obs]  += dtt_val  * inv4pi;
                this->_phi_dx_rad[obs]  += dtx_val  * inv4pi;
                this->_phi_dy_rad[obs]  += dty_val  * inv4pi;
                this->_phi_dz_rad[obs]  += dtz_val  * inv4pi;
            }
        }
    }
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::solve( )
{
    const int np = this->_n_panels;
    std::cout << "  [DBG solve] entry np=" << np << "\n" << std::flush;

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
    std::cout << "  [DBG solve] Solve(sigma) done\n" << std::flush;

    // Broadcast sigma to all processes for consistency
    MPI_Bcast(
                this->_sigma,
                np,
                mpi_cusfloat,
                this->_mpi_config->proc_root,
                MPI_COMM_WORLD
             );
    std::cout << "  [DBG solve] Bcast(sigma) done\n" << std::flush;

    // -----------------------------------------------------------------
    // Solve for sigma_dt: A * sigma_dt = rhs_dt
    // The LU factorization of A is already stored in _sysmat from the
    // pgesv call above, so we use pgetrs (back-substitution only) to
    // avoid a second O(n³) factorization.
    // -----------------------------------------------------------------
    copy_vector( np, this->_rhs_dt, this->_sigma_dt );

    this->_solver->SolveWithStoredLU( this->_sysmat, this->_sigma_dt );
    std::cout << "  [DBG solve] SolveWithStoredLU(sigma_dt) done\n" << std::flush;

    MPI_Bcast(
                this->_sigma_dt,
                np,
                mpi_cusfloat,
                this->_mpi_config->proc_root,
                MPI_COMM_WORLD
             );
    std::cout << "  [DBG solve] done\n" << std::flush;
}


template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::compute_potential_derivatives( )
{
    // _phi_dt_rad/dx/dy/dz already hold the Duhamel (memory-kernel) contribution
    // assembled in build_rhs().  Here we:
    //   (a) fill _phi_dt_wave/dx/dy/dz from the analytical incident wave potential,
    //   (b) add the steady Rankine (σ × Green) contribution to the radiation arrays,
    //   (c) set the total arrays: _phi_dt = _phi_dt_rad + _phi_dt_wave  (and likewise
    //       for dx, dy, dz).
    // This keeps radiation and incident-wave contributions separate for force
    // decomposition in the time solver.

    const int np = this->_n_panels;
    std::cout << "  [DBG cpd] entry np=" << np << "\n" << std::flush;

    // ----------------------------------------------------------------
    // (a) Incident wave contributions (skipped when no wave is active)
    //     Fills the wave-split arrays; wave arrays are zeroed first.
    // ----------------------------------------------------------------
    {
        clear_vector( np, this->_phi_dt_wave );
        clear_vector( np, this->_phi_dx_wave );
        clear_vector( np, this->_phi_dy_wave );
        clear_vector( np, this->_phi_dz_wave );

        const cusfloat w  = this->_input->ang_freq;
        const cusfloat aw = this->_input->wave_amp;
        const cusfloat h  = this->_input->water_depth;
        const cusfloat g  = this->_input->grav_acc;
        const cusfloat mu = this->_input->wave_heading * PI / static_cast<cusfloat>( 180.0 );
        const cusfloat t  = this->_t_current;

        if ( aw > static_cast<cusfloat>( 0.0 ) && w > static_cast<cusfloat>( 0.0 ) )
        {
            const cusfloat k = w2k( w, h, g );
            for ( int j=0; j<np; j++ )
            {
                const cusfloat x = this->_mesh_gp->panels[j]->center[0];
                const cusfloat y = this->_mesh_gp->panels[j]->center[1];
                const cusfloat z = this->_mesh_gp->panels[j]->center[2];

                this->_phi_dt_wave[j] += wave_potential_fo_time_dt( aw, w, k, h, g, x, y, z, mu, t );
                this->_phi_dx_wave[j] += wave_potential_fo_time_dx( aw, w, k, h, g, x, y, z, mu, t );
                this->_phi_dy_wave[j] += wave_potential_fo_time_dy( aw, w, k, h, g, x, y, z, mu, t );
                this->_phi_dz_wave[j] += wave_potential_fo_time_dz( aw, w, k, h, g, x, y, z, mu, t );
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
        std::cout << "  [DBG cpd] steady Green loop start\n" << std::flush;
        for ( int i=this->_solver->start_col_0; i<this->_solver->end_col_0; i++ )
        {
            PanelGeom* panel_i_s        = this->_mesh_gp->panels[i];
            PanelGeom* panel_mir_i      = this->_mesh_gp->panels_mirror[i];
            const cusfloat sigma_i      = this->_sigma[i];
            const cusfloat sigma_i_dt   = this->_sigma_dt[i];

            // Inner loop: all collocation rows
            for ( int j=0; j<np; j++ )
            {
                clear_vector( ndim, vel_0 );
                clear_vector( ndim, vel_1 );
                clear_vector( ndim, vel_sf );
                pot_0 = 0.0;  pot_1 = 0.0;
                PanelGeom* panel_j_s   = this->_mesh_gp->panels[j];

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
                    this->_acc_phi_dx[j] += sigma_i * half * panel_j_s->normal_vec[0];
                    this->_acc_phi_dy[j] += sigma_i * half * panel_j_s->normal_vec[1];
                    this->_acc_phi_dz[j] += sigma_i * half * panel_j_s->normal_vec[2];
                    // pot_total is singular on the diagonal: no phi_dt contribution
                }
                else
                {
                    const cusfloat inv4pi = static_cast<cusfloat>( 1.0 )
                                         / static_cast<cusfloat>( 4.0 * PI );
                    this->_acc_phi_dx[j] += sigma_i    * vel_sf[0] * inv4pi;
                    this->_acc_phi_dy[j] += sigma_i    * vel_sf[1] * inv4pi;
                    this->_acc_phi_dz[j] += sigma_i    * vel_sf[2] * inv4pi;
                    this->_acc_phi_dt[j] += sigma_i_dt * ( pot_0 - pot_1 ) * inv4pi;
                }
            }
        }

        // Reduce partial sums across all MPI processes
        std::cout << "  [DBG cpd] steady loop done, MPI_Allreduce start\n" << std::flush;
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dt, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dx, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dy, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );
        MPI_Allreduce( MPI_IN_PLACE, this->_acc_phi_dz, np, mpi_cusfloat, MPI_SUM, MPI_COMM_WORLD );

        // Accumulate steady Rankine into the radiation arrays and compute totals
        for ( int j=0; j<np; j++ )
        {
            // Steady Rankine (σ-driven) is part of the radiation potential
            this->_phi_dt_rad[j] += this->_acc_phi_dt[j];
            this->_phi_dx_rad[j] += this->_acc_phi_dx[j];
            this->_phi_dy_rad[j] += this->_acc_phi_dy[j];
            this->_phi_dz_rad[j] += this->_acc_phi_dz[j];
            // Total = radiation + incident wave
            this->_phi_dt[j] = this->_phi_dt_rad[j] + this->_phi_dt_wave[j];
            this->_phi_dx[j] = this->_phi_dx_rad[j] + this->_phi_dx_wave[j];
            this->_phi_dy[j] = this->_phi_dy_rad[j] + this->_phi_dy_wave[j];
            this->_phi_dz[j] = this->_phi_dz_rad[j] + this->_phi_dz_wave[j];
        }
    }

    std::cout << "  [DBG cpd] done\n" << std::flush;
}


// =============================================================================
// VTU export (debug visualisation in ParaView)
// =============================================================================

template<std::size_t N, int NGPT>
void FormulationKernelBackendT<N, NGPT>::export_vtu( const std::string& filename ) const
{
    // Only root rank writes; all other ranks hold identical data after solve()
    if ( this->_mpi_config->proc_rank != this->_mpi_config->proc_root ) { return; }

    const int np = this->_n_panels;

    // Count total points across all panels (panels may duplicate shared corners)
    int total_points = 0;
    for ( int ip = 0; ip < np; ip++ )
    {
        total_points += this->_mesh_gp->panels[ip]->num_nodes;
    }

    std::ofstream f( filename );
    if ( !f.is_open( ) )
    {
        std::cerr << "[export_vtu] Cannot open file: " << filename << "\n";
        return;
    }

    f << std::scientific << std::setprecision( 10 );

    // ---- VTU header --------------------------------------------------------
    f << "<?xml version=\"1.0\"?>\n";
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    f << "  <UnstructuredGrid>\n";
    f << "    <Piece NumberOfPoints=\"" << total_points
      << "\" NumberOfCells=\"" << np << "\">\n";

    // ---- Points (panel corner coordinates) ---------------------------------
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        PanelGeom* panel = this->_mesh_gp->panels[ip];
        for ( int k = 0; k < panel->num_nodes; k++ )
        {
            f << "          "
              << static_cast<double>( panel->x[k] ) << " "
              << static_cast<double>( panel->y[k] ) << " "
              << static_cast<double>( panel->z[k] ) << "\n";
        }
    }
    f << "        </DataArray>\n";
    f << "      </Points>\n";

    // ---- Cells -------------------------------------------------------------
    f << "      <Cells>\n";

    // connectivity
    f << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    {
        int pt_idx = 0;
        for ( int ip = 0; ip < np; ip++ )
        {
            const int nn = this->_mesh_gp->panels[ip]->num_nodes;
            f << "         ";
            for ( int k = 0; k < nn; k++ )
            {
                f << " " << ( pt_idx + k );
            }
            f << "\n";
            pt_idx += nn;
        }
    }
    f << "        </DataArray>\n";

    // offsets (cumulative node count per cell)
    f << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    {
        int offset = 0;
        for ( int ip = 0; ip < np; ip++ )
        {
            offset += this->_mesh_gp->panels[ip]->num_nodes;
            f << "          " << offset << "\n";
        }
    }
    f << "        </DataArray>\n";

    // cell types: 5 = VTK_TRIANGLE, 9 = VTK_QUAD
    f << "        <DataArray type=\"Int8\" Name=\"types\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        const int nn       = this->_mesh_gp->panels[ip]->num_nodes;
        const int vtk_type = ( nn == 3 ) ? 5 : 9;
        f << "          " << vtk_type << "\n";
    }
    f << "        </DataArray>\n";
    f << "      </Cells>\n";

    // ---- CellData ----------------------------------------------------------
    f << "      <CellData>\n";

    // rhs (total)
    f << "        <DataArray type=\"Float64\" Name=\"rhs\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        f << "          " << static_cast<double>( this->_rhs[ip] ) << "\n";
    }
    f << "        </DataArray>\n";

    // rhs_body_kin (body motion BC contribution)
    f << "        <DataArray type=\"Float64\" Name=\"rhs_body_kin\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        f << "          " << static_cast<double>( this->_rhs_body_kin[ip] ) << "\n";
    }
    f << "        </DataArray>\n";

    // rhs_duhamel (Duhamel convolution contribution)
    f << "        <DataArray type=\"Float64\" Name=\"rhs_duhamel\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        f << "          " << static_cast<double>( this->_rhs_duhamel[ip] ) << "\n";
    }
    f << "        </DataArray>\n";

    // rhs_wave (incident wave diffraction BC contribution)
    f << "        <DataArray type=\"Float64\" Name=\"rhs_wave\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        f << "          " << static_cast<double>( this->_rhs_wave[ip] ) << "\n";
    }
    f << "        </DataArray>\n";

    // sigma (source intensity)
    f << "        <DataArray type=\"Float64\" Name=\"sigma\" format=\"ascii\">\n";
    for ( int ip = 0; ip < np; ip++ )
    {
        f << "          " << static_cast<double>( this->_sigma[ip] ) << "\n";
    }
    f << "        </DataArray>\n";

    f << "      </CellData>\n";

    // ---- Footer ------------------------------------------------------------
    f << "    </Piece>\n";
    f << "  </UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    f.close( );
    std::cout << "[export_vtu] Written: " << filename << "\n" << std::flush;
}
