
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
#include "euler_transforms.hpp"


void    euler_rpy(
                                                                cusfloat    r,
                                                                cusfloat    p,
                                                                cusfloat    y,
                                                                cusfloat*   mat
                                            )
{
    // First row
    mat[0]  =   std::cos( p ) * std::cos( y );
    mat[1]  =  -std::cos( p ) * std::sin( y );
    mat[2]  =   std::sin( p );

    // Second row
    mat[3]  =   std::cos( r ) * std::sin( y ) + std::cos( y ) * std::sin( p ) * std::sin( r );
    mat[4]  =   std::cos( r ) * std::cos( y ) - std::sin( p ) * std::sin( r ) * std::sin( y );
    mat[5]  =  -std::cos( p ) * std::sin( r );

    // Third row
    mat[6]  =   std::sin( r ) * std::sin( y ) - std::cos( r ) * std::cos( y ) * std::sin( p );
    mat[7]  =   std::cos( y ) * std::sin( r ) + std::cos( r ) * std::sin( p ) * std::sin( y );
    mat[8]  =   std::cos( p ) * std::cos( r );
}


void    euler_ypr(
                                                                cusfloat    r,
                                                                cusfloat    p,
                                                                cusfloat    y,
                                                                cusfloat*   mat
                                            )
{
    // First row
    mat[0]  =   std::cos( p ) * std::cos( y );
    mat[1]  =   std::cos( y ) * std::sin( p ) * std::sin( r ) - std::cos( r ) * std::sin( y );
    mat[2]  =   std::sin( r ) * std::sin( y ) + std::cos( r ) * std::cos( y ) * std::sin( p );

    // Second row
    mat[3]  =   std::cos( p ) * std::sin( y );
    mat[4]  =   std::cos( r ) * std::cos( y ) + std::sin( p ) * std::sin( r ) * std::sin( y );
    mat[5]  =   std::cos( r ) * std::sin( p ) * std::sin( y ) - std::cos( y ) * std::sin(r);
    
    // Third row
    mat[6]  =  -std::sin( p );
    mat[7]  =   std::cos( p ) * std::sin( r );
    mat[8]  =   std::cos( p ) * std::cos( r );
}