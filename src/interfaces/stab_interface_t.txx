
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
#include "stab_interface_t.hpp"


template<std::size_t N>
StabInterfaceT<N>::StabInterfaceT(
                                    cusfloat    water_density,
                                    cusfloat    grav_acc,
                                    cusfloat*   ref_point
                                )
{
    // Storage local modules
    this->_ref_point        = ref_point;
    this->_grav_acc         = grav_acc;
    this->_water_density    = water_density;
}


template<std::size_t N>
void    StabInterfaceT<N>::operator()( 
                                        cusfloat*   ,
                                        cusfloat*   ,
                                        cusfloat*   xi,
                                        cusfloat*   eta,
                                        cusfloat*   zeta,
                                        bool        
                                    )
{
    // Calculate local constants
    cusfloat rhog   = this->_water_density * this->_grav_acc;

    // Calculate force modulus
    LOOP_DEF( N, this->force[i] = rhog * zeta[i]; )
    
    // Calculate moment x
    LOOP_DEF( N, this->mom_x[i] = this->force[i] * ( xi[i] - this->_ref_point[0] ); )

    // Calculate moment y
    LOOP_DEF( N, this->mom_y[i] = this->force[i] * ( eta[i] - this->_ref_point[1] ); )

    // Calculate moment z
    LOOP_DEF( N, this->mom_z[i] = this->force[i] * ( zeta[i] - this->_ref_point[2] ); )
    
}