
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
#include "../waves/wave_dispersion_fo.hpp"


template<std::size_t N>
struct GWTFcnsInterfaceT
{
protected:
    // Define protected attributes
    cusfloat                    _beta[N];
    cusfloat                    _dG_dn[N];
    cusfloat                    _dX[N];
    cusfloat                    _dY[N];
    cusfloat                    _dZ[N];
    cusfloat                    _dZp[N];
    cusfloat                    _field_point_j[3]       = { 0.0, 0.0, 0.0 };
    cusfloat                    _ftab[N];
    cusfloat                    _ftab_dmu[N];
    cusfloat                    _ftab_dt[N];
    cusfloat                    _grav_acc               = 9.81;
    cusfloat                    _mu[N];
    cusfloat                    _lt[N];
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
    cusfloat    G[N];         // Green function value
    cusfloat    dG_dx[N];     // Green function x derivative
    cusfloat    dG_dy[N];     // Green function y derivative
    cusfloat    dG_dz[N];     // Green function z derivative

    // Define constructors and destructors
    GWTFcnsInterfaceT( )   = default;

    GWTFcnsInterfaceT( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    time_ref,
                                    cusfloat    grav_acc
                );

    // Define class methods
    template<auto Kernel>
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
                                    cuscomplex  source_value

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
