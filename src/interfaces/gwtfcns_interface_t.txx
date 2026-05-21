
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
        this->_lt[i]  = 2.0 * std::sqrt( this->_grav_acc / this->_R[i] );
        this->_lt2[i] = 2.0 * this->_grav_acc / this->_R2[i];
    }

    // Calculate beta and mu
    for ( std::size_t i=0; i<N; i++ )
    {
        this->_beta[i] = this->_lt[i] * this->_time_diff / 2.0;
        this->_mu[i]   = - this->_dZp[i] / this->_R[i];
    }

    // Calculate tabulated integrals
    clear_vector<cusfloat, N>( this->_ftab );
    clear_vector<cusfloat, N>( this->_ftab_dmu );
    clear_vector<cusfloat, N>( this->_ftab_dt );

    eval_dGdt_vec<N, STATIC_LOOP_ON>( N, this->_beta, this->_mu, this->_ftab );
    eval_dGdtx_vec<N, STATIC_LOOP_ON>( N, this->_beta, this->_mu, this->_ftab_dmu );
    eval_dGdtt_vec<N, STATIC_LOOP_ON>( N, this->_beta, this->_mu, this->_ftab_dt );
    eval_dGdttx_vec<N, STATIC_LOOP_ON>( N, this->_beta, this->_mu, this->_ftab_dtmu );
    eval_dGdttt_vec<N, STATIC_LOOP_ON>( N, this->_beta, this->_mu, this->_ftab_dtt );

    // Calculate X, Y and Z cartesian coordinates derivatives
    cusfloat a = 0.0;
    cusfloat b = 0.0;
    cusfloat c = 0.0;
    for ( std::size_t i=0; i<N; i++ )
    {
        // Calculate first order time derivative
        a                   = - this->_lt[i] / this->_R2[i];
        b                   = -1.5 * this->_ftab[i];
        c                   = - 0.5 * this->_beta[i] * this->_ftab_dt[i] + this->_ftab_dmu[i] * this->_dZp[i] / this->_R[i];
        this->dG_dtx[i]     = a * ( b + this->_dX[i] * c );
        this->dG_dty[i]     = a * ( b + this->_dY[i] * c );
        this->dG_dtz[i]     = a * ( b + this->_dZ[i] * c ) + this->_lt[i] * this->_ftab_dmu[i] / this->_R[i];

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

    // Calculate normal derivate
    cusfloat nx_pf = this->_source_i->normal_vec[0];
    cusfloat ny_pf = this->_source_i->normal_vec[1];
    cusfloat nz_pf = this->_source_i->normal_vec[2];

    for ( std::size_t i=0; i<N; i++ )
    {
        this->dG_dtn[i] = (
                            - 
                            this->dG_dtx[i] * nx_pf
                            -
                            this->dG_dty[i] * ny_pf
                            +
                            this->dG_dtz[i] * nz_pf
                        );

        this->dG_dttn[i] = (
                            - 
                            this->dG_dttx[i] * nx_pf
                            -
                            this->dG_dtty[i] * ny_pf
                            +
                            this->dG_dttz[i] * nz_pf
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