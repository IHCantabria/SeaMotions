
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
    cusfloat    dG_dtn[N];      // Green function first order time normal derivative
    cusfloat    dG_dtx[N];      // Green function first order time x derivative
    cusfloat    dG_dty[N];      // Green function first order time y derivative
    cusfloat    dG_dtz[N];      // Green function first order time z derivative
    cusfloat    dG_dtt[N];      // Green function second order time derivative
    cusfloat    dG_dttn[N];     // Green function second order time normal derivative
    cusfloat    dG_dttx[N];     // Green function second order time x derivative
    cusfloat    dG_dtty[N];     // Green function second order time y derivative
    cusfloat    dG_dttz[N];     // Green function second order time z derivative

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

};

// Include function definitions
#include "gwtfcns_interface_t.txx"
