
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

// Include local modules
#include "gwtfcns_interface_t.hpp"
#include "../math/shape_functions.hpp"
#include "../math/math_tools.hpp"
#include "../math/math_interface.hpp"
#include "../green/time_domain_evaluator.hpp"
#include "../green/td_data/td_database.hpp"
#include "../green/time_domain_asymptotic.hpp"

// Timer shorthand used throughout this file
using _Clock = std::chrono::high_resolution_clock;


template<std::size_t N>
GWTFcnsInterfaceT<N>::GWTFcnsInterfaceT( 
                                            SourceNode* source_i,
                                            SourceNode* source_j,
                                            cusfloat    time_ref,
                                            cusfloat    grav_acc
                                        )
{
    this->_grav_acc     = grav_acc;
    this->_source_i     = source_i;
    this->_source_j     = source_j;
    this->_time_diff    = time_ref;
}


template<std::size_t N>
void        GWTFcnsInterfaceT<N>::operator()( 
                                                        cusfloat*   ,
                                                        cusfloat*   ,
                                                        cusfloat*   xi,
                                                        cusfloat*   eta,
                                                        cusfloat*   zeta,
                                                        bool        verbose
                                                    )
{
    ++_timer_call_count;

    // ----------------------------------------------------------------
    // Section 1 – Geometry: dX, dY, dZ, R, R2, R3
    // ----------------------------------------------------------------
    auto _t0 = _Clock::now();

    cusfloat source_x = this->_source_j->position[0];
    cusfloat source_y = this->_source_j->position[1];
    cusfloat source_z = this->_source_j->position[2];

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_dX[i]    = source_x - xi[i];
        this->_dY[i]    = source_y - eta[i];
        this->_z[i]     = source_z;
        this->_dZ[i]    = source_z - zeta[i];
        this->_dZp[i]   = source_z + zeta[i];
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_R[i] = ( pow2s( this->_dX[i] ) + pow2s( this->_dY[i] ) + pow2s( this->_dZ[i] ) );
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_R[i] = std::max( this->_R[i], 1e-12 );
    }

    lv_sqrt<cusfloat>( N, this->_R, this->_R );

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_R2[i] = this->_R[i] * this->_R[i];
        this->_R3[i] = this->_R2[i] * this->_R[i];
    }

    _timer_geometry += _Clock::now() - _t0;

    // ----------------------------------------------------------------
    // Section 2 – Leading term, beta and mu
    // ----------------------------------------------------------------
    _t0 = _Clock::now();

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_lt[i]  = 2.0 * std::sqrt( this->_grav_acc / this->_R3[i] );   // 2*sqrt(g/R^3)
        this->_lt2[i] = 2.0 * this->_grav_acc / this->_R2[i];                // 2*g/R^2 = d/dt( _lt * sqrt(g*R) )
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_beta[i] = std::sqrt( this->_grav_acc / this->_R[i] ) * this->_time_diff;   // sqrt(g/R)*time_diff
        this->_mu[i]   = - this->_dZp[i] / this->_R[i];
    }

    _timer_leading += _Clock::now() - _t0;

    // ----------------------------------------------------------------
    // Section 3 – Tabulated integrals (bilinear interpolation from embedded database)
    //
    // Five static bilinear interpolators (thread-safe C++11 magic statics) are
    // built once on the first call and reused thereafter.  They store the
    // COMPLETE table value (G0 + residual), so no separate G0 correction is
    // needed.  _timer_eval_dGdt accumulates the total bilinear evaluation time;
    // the individual sub-timers (dGdtx, dGdtt, etc.) are no longer incremented.
    // ----------------------------------------------------------------

    static const Bilinear2D<double> s_bilin_Gt   = td_db::make_bilinear_Gt();
    static const Bilinear2D<double> s_bilin_Gtx  = td_db::make_bilinear_Gtx();
    static const Bilinear2D<double> s_bilin_Gtt  = td_db::make_bilinear_Gtt();
    static const Bilinear2D<double> s_bilin_Gttx = td_db::make_bilinear_Gttx();
    static const Bilinear2D<double> s_bilin_Gttt = td_db::make_bilinear_Gttt();

    // β threshold above which the asymptotic expansion supersedes the database.
    // Results from the asymptotic functions are divided by 2 to match the
    // database scaling convention.
    constexpr double BETA_ASYMP_THRESHOLD = 50.0;
    constexpr cusfloat ASYMP_HALF = cusfloat(0.5);

    _t0 = _Clock::now();
    for ( std::size_t i = 0; i < N; ++i )
    {
        const double log_mu = std::log10( static_cast<double>( this->_mu[i]   ) );
        const double beta_i =             static_cast<double>( this->_beta[i] );

        if ( beta_i > BETA_ASYMP_THRESHOLD )
        {
            // Asymptotic expansion valid for β > 50 (divided by 2 to match database)
            // this->_ftab    [i] = 1.0;
            this->_ftab    [i] = ASYMP_HALF * dGdt_asymptotic  ( this->_beta[i], this->_mu[i] );
            this->_ftab_dmu[i] = ASYMP_HALF * dGdtx_asymptotic ( this->_beta[i], this->_mu[i] );
            this->_ftab_dt [i] = ASYMP_HALF * dGdtt_asymptotic ( this->_beta[i], this->_mu[i] );
            this->_ftab_dtmu[i]= ASYMP_HALF * dGdttx_asymptotic( this->_beta[i], this->_mu[i] );
            this->_ftab_dtt[i] = ASYMP_HALF * dGdttt_asymptotic( this->_beta[i], this->_mu[i] );
        }
        else
        {
            // Tabulated bilinear interpolation for β ≤ 50
            // this->_ftab    [i] = 1.0;
            this->_ftab    [i] = static_cast<cusfloat>( s_bilin_Gt  .eval( beta_i, log_mu ) );
            this->_ftab_dmu[i] = static_cast<cusfloat>( s_bilin_Gtx .eval( beta_i, log_mu ) );
            this->_ftab_dt [i] = static_cast<cusfloat>( s_bilin_Gtt .eval( beta_i, log_mu ) );
            this->_ftab_dtmu[i]= static_cast<cusfloat>( s_bilin_Gttx.eval( beta_i, log_mu ) );
            this->_ftab_dtt[i] = static_cast<cusfloat>( s_bilin_Gttt.eval( beta_i, log_mu ) );
        }
    }
    _timer_eval_dGdt += _Clock::now() - _t0;

    // ----------------------------------------------------------------
    // Section 4 – Cartesian derivatives
    // ----------------------------------------------------------------
    _t0 = _Clock::now();

    cusfloat a = 0.0;
    cusfloat c = 0.0;
    for ( std::size_t i=0; i<N; i++ )
    {
        // Calculate first order time derivative
        // dG_dt = lt * Gt
        // d(dG_dt)/d(src_k) = a * dK * c   (k = x,y; for k=z an extra term appears)
        // where:
        //   a = -lt / R^2
        //   c = 0.5*Gt + 0.5*beta*Gtt - Gtx*dZp/R   (= 0.5*Gt + 0.5*beta*Gtt + mu*Gtx)
        //   extra z-term: -lt * Gtx / R  (from d(mu)/d(sz) = -1/R - mu*dZ/R^2)
        a                   = - this->_lt[i] / this->_R2[i];
        c                   = 1.5 * this->_ftab[i] + 0.5 * this->_beta[i] * this->_ftab_dt[i] - this->_ftab_dmu[i] * this->_dZp[i] / this->_R[i];
        this->dG_dt[i]      = this->_lt[i] * this->_ftab[i];
        this->dG_dtx[i]     = a * this->_dX[i] * c;
        this->dG_dty[i]     = a * this->_dY[i] * c;
        this->dG_dtz[i]     = a * this->_dZ[i] * c - this->_lt[i] * this->_ftab_dmu[i] / this->_R[i];

        // Calculate second order time derivative
        a                   = - this->_lt2[i] * ( 
                                                    - 2.0 * this->_ftab_dt[i]
                                                    - 0.5 * this->_beta[i] * this->_ftab_dtt[i]
                                                    + this->_ftab_dtmu[i] * this->_dZp[i] / this->_R[i]
                                                ) / this->_R2[i];
        this->dG_dtt[i]     = - this->_lt2[i] * this->_ftab_dt[i];
        this->dG_dttx[i]    = a * this->_dX[i];
        this->dG_dtty[i]    = a * this->_dY[i];
        this->dG_dttz[i]    = a * this->_dZ[i] + this->_lt2[i] * this->_ftab_dtmu[i] / this->_R[i];

    }

    _timer_derivatives += _Clock::now() - _t0;

    // ----------------------------------------------------------------
    // Verbose check – numerical (central-difference) verification of
    // dG_dtx / dG_dty / dG_dtz against finite differences of dG_dt
    // with respect to the source position.
    // ----------------------------------------------------------------
    if ( false /*verbose*/ )
    {
        // Recompute dG_dt at a perturbed source position.
        // s_bilin_Gt, BETA_ASYMP_THRESHOLD and ASYMP_HALF are still in scope.
        auto compute_dGdt_at = [&]( cusfloat sx, cusfloat sy, cusfloat sz,
                                    std::size_t k ) -> cusfloat
        {
            const cusfloat ddX  = sx - xi[k];
            const cusfloat ddY  = sy - eta[k];
            const cusfloat ddZ  = sz - zeta[k];
            const cusfloat ddZp = sz + zeta[k];
            cusfloat r2 = ddX*ddX + ddY*ddY + ddZ*ddZ;
            r2 = std::max( r2, cusfloat(1e-12) );
            const cusfloat r    = std::sqrt( r2 );
            const cusfloat lt   = cusfloat(2) * std::sqrt( this->_grav_acc / r );
            const cusfloat beta = lt * this->_time_diff / cusfloat(2);
            const cusfloat mu_  = -ddZp / r;
            const double   lm   = std::log10( static_cast<double>( mu_ ) );
            const double   bd   = static_cast<double>( beta );
            cusfloat ftab;
            if ( bd > BETA_ASYMP_THRESHOLD )
                ftab = ASYMP_HALF * dGdt_asymptotic( beta, mu_ );
            else
                ftab = static_cast<cusfloat>( s_bilin_Gt.eval( bd, lm ) );
            return lt * ftab;
        };

        std::printf( "\n--- dG_dt Cartesian derivative check (central differences on source position) ---\n" );
        std::printf( "  %4s  %13s  %13s  %10s  %13s  %13s  %10s  %13s  %13s  %10s\n",
                     "i",
                     "dG_dtx_an", "dG_dtx_fd", "err_x",
                     "dG_dty_an", "dG_dty_fd", "err_y",
                     "dG_dtz_an", "dG_dtz_fd", "err_z" );

        for ( std::size_t i = 0; i < N; i++ )
        {
            const cusfloat h = std::max( cusfloat(1e-6) * this->_R[i], cusfloat(1e-8) );

            const cusfloat fd_x = ( compute_dGdt_at( source_x + h, source_y,     source_z,     i )
                                  - compute_dGdt_at( source_x - h, source_y,     source_z,     i ) )
                                / ( cusfloat(2) * h );
            const cusfloat fd_y = ( compute_dGdt_at( source_x,     source_y + h, source_z,     i )
                                  - compute_dGdt_at( source_x,     source_y - h, source_z,     i ) )
                                / ( cusfloat(2) * h );
            const cusfloat fd_z = ( compute_dGdt_at( source_x,     source_y,     source_z + h, i )
                                  - compute_dGdt_at( source_x,     source_y,     source_z - h, i ) )
                                / ( cusfloat(2) * h );

            std::printf( "  %4zu  %13.5e  %13.5e  %10.3e  %13.5e  %13.5e  %10.3e  %13.5e  %13.5e  %10.3e\n",
                         i,
                         static_cast<double>( this->dG_dtx[i] ),
                         static_cast<double>( fd_x ),
                         std::abs( static_cast<double>( this->dG_dtx[i] - fd_x ) ),
                         static_cast<double>( this->dG_dty[i] ),
                         static_cast<double>( fd_y ),
                         std::abs( static_cast<double>( this->dG_dty[i] - fd_y ) ),
                         static_cast<double>( this->dG_dtz[i] ),
                         static_cast<double>( fd_z ),
                         std::abs( static_cast<double>( this->dG_dtz[i] - fd_z ) ) );
        }
        std::printf( "----------------------------------------------------------------------------------\n\n" );
    }

    // ----------------------------------------------------------------
    // Section 5 – Normal derivatives
    // ----------------------------------------------------------------
    _t0 = _Clock::now();

    cusfloat nx_pf = this->_source_j->normal_vec[0];
    cusfloat ny_pf = this->_source_j->normal_vec[1];
    cusfloat nz_pf = this->_source_j->normal_vec[2];

    for ( std::size_t i=0; i<N; i++ )
    {
        this->dG_dtn[i] = (
                            + 
                            this->dG_dtx[i] * nx_pf
                            +
                            this->dG_dty[i] * ny_pf
                            +
                            this->dG_dtz[i] * nz_pf
                        );

        this->dG_dttn[i] = (
                            + 
                            this->dG_dttx[i] * nx_pf
                            +
                            this->dG_dtty[i] * ny_pf
                            +
                            this->dG_dttz[i] * nz_pf
                        );
    }

    _timer_normals += _Clock::now() - _t0;

}


template<std::size_t N>
void    GWTFcnsInterfaceT<N>::set_field_point(
                                            cusfloat*   fp
                                        )
{
    this->_field_point_j[0] = fp[0];
    this->_field_point_j[1] = fp[1];
    this->_field_point_j[2] = fp[2];
}


template<std::size_t N>
void    GWTFcnsInterfaceT<N>::set_source_i(
                                            SourceNode* source_node,
                                            cusfloat    source_value
                                    )
{
    this->_source_i     = source_node;
    this->_source_value = source_value;
}


template<std::size_t N>
void    GWTFcnsInterfaceT<N>::set_source_j(
                                                SourceNode* source_node
                                            )
{
    this->_source_j = source_node;
}


template<std::size_t N>
void    GWTFcnsInterfaceT<N>::set_time_diff(
                                                cusfloat    time_diff
                                            )
{
    this->_time_diff = time_diff;
}


template<std::size_t N>
void    GWTFcnsInterfaceT<N>::print_timers() const
{
    // Total accumulated time across all sections
    double t_geo   = _timer_geometry.count();
    double t_lead  = _timer_leading.count();
    double t_dGdt  = _timer_eval_dGdt.count();
    double t_dGdtx = _timer_eval_dGdtx.count();
    double t_dGdtt = _timer_eval_dGdtt.count();
    double t_dGdttx= _timer_eval_dGdttx.count();
    double t_dGdttt= _timer_eval_dGdttt.count();
    double t_tab   = t_dGdt + t_dGdtx + t_dGdtt + t_dGdttx + t_dGdttt;
    double t_der   = _timer_derivatives.count();
    double t_nor   = _timer_normals.count();
    double t_tot   = t_geo + t_lead + t_tab + t_der + t_nor;

    auto pct = [&]( double t ) -> double {
        return ( t_tot > 0.0 ) ? 100.0 * t / t_tot : 0.0;
    };

    std::printf( "\n--- GWTFcnsInterfaceT<%zu> operator() timing report ---\n",  N );
    std::printf( "  Calls              : %zu\n",      _timer_call_count );
    std::printf( "  Total              : %.6f s\n",   t_tot );
    std::printf( "  Geometry           : %.6f s  (%5.1f %%)\n", t_geo,    pct(t_geo)    );
    std::printf( "  Leading/beta/mu    : %.6f s  (%5.1f %%)\n", t_lead,   pct(t_lead)   );
    std::printf( "  Tabulated (total)  : %.6f s  (%5.1f %%)\n", t_tab,    pct(t_tab)    );
    std::printf( "    eval_dGdt        : %.6f s  (%5.1f %%)\n", t_dGdt,   pct(t_dGdt)   );
    std::printf( "    eval_dGdtx       : %.6f s  (%5.1f %%)\n", t_dGdtx,  pct(t_dGdtx)  );
    double t_dGdtt_logmu  = _timer_eval_dGdtt_logmu.count();
    double t_dGdtt_cheby  = _timer_eval_dGdtt_cheby.count();
    double t_dGdtt_G0     = _timer_eval_dGdtt_G0.count();
    double t_dGdttx_logmu = _timer_eval_dGdttx_logmu.count();
    double t_dGdttx_cheby = _timer_eval_dGdttx_cheby.count();
    double t_dGdttx_G0    = _timer_eval_dGdttx_G0.count();
    std::printf( "    eval_dGdtt       : %.6f s  (%5.1f %%)\n", t_dGdtt,       pct(t_dGdtt)       );
    std::printf( "      log_mu         : %.6f s  (%5.1f %%)\n", t_dGdtt_logmu, pct(t_dGdtt_logmu) );
    std::printf( "      cheby eval     : %.6f s  (%5.1f %%)\n", t_dGdtt_cheby, pct(t_dGdtt_cheby) );
    std::printf( "      G0 correction  : %.6f s  (%5.1f %%)\n", t_dGdtt_G0,    pct(t_dGdtt_G0)    );
    std::printf( "    eval_dGdttx      : %.6f s  (%5.1f %%)\n", t_dGdttx,      pct(t_dGdttx)      );
    std::printf( "      log_mu         : %.6f s  (%5.1f %%)\n", t_dGdttx_logmu,pct(t_dGdttx_logmu));
    std::printf( "      cheby eval     : %.6f s  (%5.1f %%)\n", t_dGdttx_cheby,pct(t_dGdttx_cheby));
    std::printf( "      G0 correction  : %.6f s  (%5.1f %%)\n", t_dGdttx_G0,   pct(t_dGdttx_G0)   );
    std::printf( "    eval_dGdttt      : %.6f s  (%5.1f %%)\n", t_dGdttt, pct(t_dGdttt) );
    std::printf( "  Derivatives        : %.6f s  (%5.1f %%)\n", t_der,    pct(t_der)    );
    std::printf( "  Normal derivs      : %.6f s  (%5.1f %%)\n", t_nor,    pct(t_nor)    );
    std::printf( "-------------------------------------------------------\n\n" );
}


template<std::size_t N>
void    GWTFcnsInterfaceT<N>::reset_timers()
{
    _timer_geometry     = Duration::zero();
    _timer_leading      = Duration::zero();
    _timer_eval_dGdt    = Duration::zero();
    _timer_eval_dGdtx   = Duration::zero();
    _timer_eval_dGdtt        = Duration::zero();
    _timer_eval_dGdtt_logmu  = Duration::zero();
    _timer_eval_dGdtt_cheby  = Duration::zero();
    _timer_eval_dGdtt_G0     = Duration::zero();
    _timer_eval_dGdttx       = Duration::zero();
    _timer_eval_dGdttx_logmu = Duration::zero();
    _timer_eval_dGdttx_cheby = Duration::zero();
    _timer_eval_dGdttx_G0    = Duration::zero();
    _timer_eval_dGdttt  = Duration::zero();
    _timer_derivatives  = Duration::zero();
    _timer_normals      = Duration::zero();
    _timer_call_count   = 0;
}