
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
#include "../green/integrals_db.hpp"
#include "../green/pulsating_fin_depth.hpp"
#include "../waves/wave_dispersion_fo.hpp"


template<std::size_t N>
struct GWFcnsInterfaceT
{
protected:
    // Define protected attributes
    BesselFactoryVecUpTo<N>     _bessel_factory;
    cuscomplex                  _dG_dR[N];
    cusfloat                    _dX[N];
    cusfloat                    _dY[N];
    cusfloat                    _grav_acc               = 0.0;
    cusfloat                    _field_point_j[3]       = { 0.0, 0.0, 0.0 };
    cusfloat                    _R[N];
    SourceNode*                 _source_i               = nullptr;
    SourceNode*                 _source_j               = nullptr;
    cuscomplex                  _source_value           = 0.0;
    cusfloat                    _z[N];
    cusfloat                    _water_depth            = 0.0;
    WaveDispersionFONK          _wave_data;

    // Define protected methods
    void    _initialize(
                                    cusfloat    ang_freq
                        );

public:
    // Define public attributes
    cuscomplex    G[N];         // Green function value
    cuscomplex    dG_dn_sf[N];  // Green function normal derivative
    cuscomplex    dG_dn_pf[N];  // Green function normal derivative
    cuscomplex    dG_dx[N];     // Green function x derivative
    cuscomplex    dG_dy[N];     // Green function y derivative
    cuscomplex    dG_dz[N];     // Green function z derivative
    cuscomplex    dG_dzeta[N];  // Green function zeta derivative

    // Define constructors and destructors
    GWFcnsInterfaceT( )   = default;

    GWFcnsInterfaceT( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    ang_freq,
                                    cusfloat    water_depth,
                                    cusfloat    grav_acc
                );

    // Define class methods
    void        initialize( 
                                    cusfloat    ang_freq,
                                    cusfloat    water_depth,
                                    cusfloat    grav_acc
                            );
    
    template<auto Kernel>
    void        operator()( 
                                    cusfloat*   xi,
                                    cusfloat*   eta,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    bool        verbose=false
                           );

    void        set_ang_freq(
                                    cusfloat    ang_freq
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

};

// Include function definitions
#include "gwfcns_interface_t.txx"
