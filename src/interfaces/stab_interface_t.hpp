
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
#include "../config.hpp"
#include "../static_tools.hpp"


template<std::size_t N>
struct StabInterfaceT
{
protected:
    /* Define protected attributes */
    cusfloat*   _ref_point      = nullptr;  // Reference point to take momentums from
    cusfloat    _grav_acc       = 0.0;      // Gravitational acceleration
    cusfloat    _water_density  = 0.0;      // Water density


public:
    // Define public attributes
    cuscomplex    force[N];     // Force modulus
    cuscomplex    mom_x[N];     // Moment X modulus
    cuscomplex    mom_y[N];     // Moment Y modulus
    cuscomplex    mom_z[N];     // Moment Z modulus

    // Define constructors and destructors
    StabInterfaceT( )   = default;

    StabInterfaceT( 
                        cusfloat    water_density,
                        cusfloat    grav_acc,
                        cusfloat*   cog
                    );


    void        operator()( 
                                    cusfloat*   xi,
                                    cusfloat*   eta,
                                    cusfloat*   x,
                                    cusfloat*   y,
                                    cusfloat*   z,
                                    bool        verbose=false
                           );

};

// Include function definitions
#include "stab_interface_t.txx"
