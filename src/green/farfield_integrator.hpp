
#pragma once

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
#include "../config.hpp"
#include "../inout/input.hpp"
#include "triple_hankel_far_field.hpp"


template<QTFTypeE qtf_type>
class FarFieldIntegrator
{
private:
/*** Declare private member variables ***/
    cusfloat        _alpha              = 0.0;
    std::size_t     _freq_i             = 0;
    std::size_t     _freq_j             = 0;
    cusfloat        _grav_acc           = 9.81;
    Input*          _input              = nullptr;
    cusfloat        _omega              = 0.0;
    cusfloat        _sign               = 1.0;
    TripleHankelIO  _hankel_int_opts    ;
    cusfloat        _partition_circle   = 0.0;
    cusfloat        _wi                 = 0.0;
    cusfloat        _wj                 = 0.0;
    cusfloat        _wave_i             = 0.0;
    cusfloat        _wave_j             = 0.0;
    cusfloat        _wave_k             = 0.0;

public:
    /*** Declare class constructors ***/
    FarFieldIntegrator( ) = default;

    FarFieldIntegrator( 
                                        Input*      input,
                                        std::size_t freq_i,
                                        std::size_t freq_j,
                                        cusfloat    partition_circle
                        );

    /*** Declare class constructors ***/
    void    integrate( 
                                        std::size_t     n,
                                        std::size_t     N,
                                        cuscomplex      Ac,
                                        cuscomplex      As,
                                        cuscomplex      Bc,
                                        cuscomplex      Bs,
                                        cuscomplex*     Qc,
                                        cuscomplex*     Qs
                        );
    
    void    set_frequency_indices( 
                                        std::size_t freq_i,
                                        std::size_t freq_j
                                    );
    
    template<typename T>
    cuscomplex U(
                    T                       F, 
                    int                     l,
                    int                     m
                );

    template<typename T>
    cuscomplex V( 
                    T                       F, 
                    int                     l,
                    int                     m
                );

    template<typename T>
    cuscomplex W( 
                    T                       F, 
                    int                     l,
                    int                     m
                );

};

// Include template definitions
#include "farfield_integrator.txx"