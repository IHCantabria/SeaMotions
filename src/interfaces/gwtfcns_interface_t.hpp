
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

// Include standard library
#include <atomic>
#include <chrono>
#include <cstdio>

// Include local modules
#include "../containers/source_node.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "../green/time_domain_evaluator.hpp"
#include "../math/math_tools.hpp"
#include "../waves/wave_dispersion_fo.hpp"


template<std::size_t N>
struct GWTFcnsInterfaceT
{
protected:
    // Define protected attributes
    cusfloat                    _beta[N];
    cusfloat                    _dX[N];
    cusfloat                    _dY[N];
    cusfloat                    _dZ[N];
    cusfloat                    _dZp[N];
    cusfloat                    _field_point_j[3]       = { 0.0, 0.0, 0.0 };
    cusfloat                    _ftab[N];
    cusfloat                    _ftab_dmu[N];
    cusfloat                    _ftab_dt[N];
    cusfloat                    _ftab_dtt[N];
    cusfloat                    _ftab_dtmu[N];
    cusfloat                    _grav_acc               = 9.81;
    cusfloat                    _mu[N];
    cusfloat                    _lt[N];
    cusfloat                    _lt2[N];
    cusfloat                    _R[N];
    cusfloat                    _R2[N];
    cusfloat                    _R3[N];
    SourceNode*                 _source_i               = nullptr;
    SourceNode*                 _source_j               = nullptr;
    cusfloat                    _source_value           = 0.0;
    cusfloat                    _time_diff              = 0.0;
    cusfloat                    _z[N];

public:
    // Define public attributes
    cusfloat    dG_dt[N];       // Green function first order time derivative
    cusfloat    dG_dtn[N];      // Green function first order time normal derivative
    cusfloat    dG_dtx[N];      // Green function first order time x derivative
    cusfloat    dG_dty[N];      // Green function first order time y derivative
    cusfloat    dG_dtz[N];      // Green function first order time z derivative
    cusfloat    dG_dtt[N];      // Green function second order time derivative
    cusfloat    dG_dttn[N];     // Green function second order time normal derivative
    cusfloat    dG_dttx[N];     // Green function second order time x derivative
    cusfloat    dG_dtty[N];     // Green function second order time y derivative
    cusfloat    dG_dttz[N];     // Green function second order time z derivative

    // Cumulative profiling timers (accumulated over all operator() calls)
    using Duration = std::chrono::duration<double>;
    Duration    _timer_geometry             = Duration::zero();  // R, dX, dY, dZ, R2, R3
    Duration    _timer_leading              = Duration::zero();  // leading term, beta, mu
    Duration    _timer_eval_dGdt            = Duration::zero();  // eval_dGdt_vec
    Duration    _timer_eval_dGdtx           = Duration::zero();  // eval_dGdtx_vec
    Duration    _timer_eval_dGdtt           = Duration::zero();  // eval_dGdtt_vec (total)
    Duration    _timer_eval_dGdtt_logmu     = Duration::zero();  //   |- compute_log_mu
    Duration    _timer_eval_dGdtt_cheby     = Duration::zero();  //   |- eval_time_residual_2d_vec
    Duration    _timer_eval_dGdtt_G0        = Duration::zero();  //   |- G0 correction loop
    Duration    _timer_eval_dGdttx          = Duration::zero();  // eval_dGdttx_vec (total)
    Duration    _timer_eval_dGdttx_logmu    = Duration::zero();  //   |- compute_log_mu
    Duration    _timer_eval_dGdttx_cheby    = Duration::zero();  //   |- eval_time_residual_2d_vec
    Duration    _timer_eval_dGdttx_G0       = Duration::zero();  //   |- G0 correction loop
    Duration    _timer_eval_dGdttt          = Duration::zero();  // eval_dGdttt_vec
    Duration    _timer_derivatives          = Duration::zero();  // Cartesian derivative loop
    Duration    _timer_normals              = Duration::zero();  // normal derivative loop
    std::size_t _timer_call_count           = 0;

    // Define constructors and destructors
    GWTFcnsInterfaceT( )   = default;

    GWTFcnsInterfaceT( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    time_ref,
                                    cusfloat    grav_acc
                );

    // Define class methods
    void        operator()( 
                                    cusfloat*   xi,
                                    cusfloat*   eta,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    bool        verbose=false
                           );

    void        set_field_point(
                                    cusfloat*   fp
                                );

    void        set_source_i(
                                    SourceNode* source_node,
                                    cusfloat    source_value

                            );

    void        set_source_j(
                                    SourceNode* source_node
                            );

    void        set_time_diff(
                                    cusfloat    time_diff
                            );

    void        print_timers() const;

    void        reset_timers();

    // Read-only accessors over the spatial Gauss-point cloud.
    // Used by the Duhamel-tracking diagnostic to summarise where each
    // quadrature node lives in (β, μ, R) space.
    cusfloat    min_mu()   const { cusfloat v = _mu[0];   for ( std::size_t i = 1; i < N; ++i ) if ( _mu[i]   < v ) v = _mu[i];   return v; }
    cusfloat    max_mu()   const { cusfloat v = _mu[0];   for ( std::size_t i = 1; i < N; ++i ) if ( _mu[i]   > v ) v = _mu[i];   return v; }
    cusfloat    min_beta() const { cusfloat v = _beta[0]; for ( std::size_t i = 1; i < N; ++i ) if ( _beta[i] < v ) v = _beta[i]; return v; }
    cusfloat    max_beta() const { cusfloat v = _beta[0]; for ( std::size_t i = 1; i < N; ++i ) if ( _beta[i] > v ) v = _beta[i]; return v; }
    cusfloat    min_R()    const { cusfloat v = _R[0];    for ( std::size_t i = 1; i < N; ++i ) if ( _R[i]    < v ) v = _R[i];    return v; }
    cusfloat    max_R()    const { cusfloat v = _R[0];    for ( std::size_t i = 1; i < N; ++i ) if ( _R[i]    > v ) v = _R[i];    return v; }

    // Raw per-Gauss-point read-only accessors.  Length = N for every array.
    // Used by the Duhamel tracker to dump the underlying geometry that goes
    // into beta/mu so we can spot anomalies (e.g. a GP crossing the waterline)
    // step by step.
    static constexpr std::size_t  gp_count() { return N; }
    const cusfloat*  gp_R()    const { return _R;   }
    const cusfloat*  gp_dX()   const { return _dX;  }
    const cusfloat*  gp_dY()   const { return _dY;  }
    const cusfloat*  gp_dZ()   const { return _dZ;  }
    const cusfloat*  gp_dZp()  const { return _dZp; }
    const cusfloat*  gp_beta() const { return _beta;}
    const cusfloat*  gp_mu()   const { return _mu;  }

};

// Include function definitions
#include "gwtfcns_interface_t.txx"
