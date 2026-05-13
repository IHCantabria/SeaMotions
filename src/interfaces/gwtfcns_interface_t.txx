
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
template<auto Kernel>
void        GWTFcnsInterfaceT<N>::operator()( 
                                                        cusfloat*   ,
                                                        cusfloat*   ,
                                                        cusfloat*   xi,
                                                        cusfloat*   eta,
                                                        cusfloat*   zeta,
                                                        bool        verbose
                                                    )
{
    // Calculate horizontal radius
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

    // Calculate leading term
    for ( std::size_t i=0; i<N; i++ )
    {
        this->_lt[i] = 2.0 * std::sqrt( this->_grav_acc / this->_R[i] );
    }

    // Calculate beta and mu
    for ( std::size_t i=0; i<N; i++ )
    {
        this->_beta[i] = this->_lt[i] * this->_time_diff / 2.0;
        this->_mu[i]   = - this->_dZp[i] / this->_R[i];
    }

    // Calculate tabulated integrals
    clear_vector<N>( this->_ftab );
    clear_vector<N>( this->_ftab_dmu );
    clear_vector<N>( this->_ftab_dt );

    eval_dGdt_vec<N, STATIC_LOOP_ON( N, this->_beta, this->_mu, this->_ftab );
    eval_dGdtx_vec<N, STATIC_LOOP_ON( N, this->_beta, this->_mu, this->_ftab_dmu );
    eval_dGdtt_vec<N, STATIC_LOOP_ON( N, this->_beta, this->_mu, this->_ftab_dt );

    // Calculate G
    for ( std::size_t i=0; i<N; i++ )
    {
        this->G[i] = this->_lt[i] * this->_ftab[i];
    }

    // Calculate X, Y and Z cartesian coordinates derivatives
    cusfloat dc = 0.0;
    for ( std::size_t i=0; i<N; i++ )
    {
        dc              = this->_lt[i] * (
                                                -( this->_ftab[i] + this->_beta[i] * this->_ftab_dt[i] ) / this->_R2[i]
                                                +
                                                this->_ftab_dmu[i] * ( this->_dZp[i] / this->_R3[i] );
                                            );
        this->dG_dx[i] = this->_dX[i] * dc;
        this->dG_dy[i] = this->_dY[i] * dc;
        this->dG_dz[i] = this->_dZ[i] * dc + this->_lt[i] * this->_ftab_dmu[i] / this->_R[i];
    }

    // Calculate normal derivate
    cusfloat nx_pf = this->_source_i->normal_vec[0];
    cusfloat ny_pf = this->_source_i->normal_vec[1];
    cusfloat nz_pf = this->_source_i->normal_vec[2];

    for ( std::size_t i=0; i<N; i++ )
    {
        this->dG_dn[i] = (
                            - 
                            this->dG_dx[i] * nx_pf
                            -
                            this->dG_dy[i] * ny_pf
                            +
                            this->dG_dz[i] * nz_pf
                        );
    }
    
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
                                            cuscomplex  source_value
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