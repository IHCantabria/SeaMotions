
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


/*******************************************/
/******* CONSTANTS DECLARATION  BLOCK*******/
/*******************************************/
const int MAX_LIN_NODES = 4;

/*******************************************/
/************* Define MACROS ***************/
/*******************************************/
#define JACMAT2D_DOT_PRODUCT( jac_mat, col_idx, coord, shape_der, NP, NV )      \
{                                                                               \
    for ( int i=0; i<NP; i++ )                                                  \
    {                                                                           \
        jac_mat[4*i+col_idx] = 0;                                               \
        for ( int j=0; j<NV; j++ )                                              \
        {                                                                       \
            jac_mat[4*i+col_idx] += coord[j] * shape_der[i*NV+j];               \
        }                                                                       \
    }                                                                           \
}                                                                               \

/*******************************************/
/******* FUNCTION DECLARATION  BLOCK********/
/*******************************************/
                                    cusfloat    jacobi_det_2d(
                                                                        int         np,
                                                                        cusfloat*   xn,
                                                                        cusfloat*   yn,
                                                                        cusfloat    xi,
                                                                        cusfloat    eta
                                                                );

template<int NV, int NP>    inline  void        jacobi_det_2d_vec(
                                                                        cusfloat*   xn,
                                                                        cusfloat*   yn,
                                                                        cusfloat*   xi,
                                                                        cusfloat*   eta,
                                                                        cusfloat*   jac_det
                                                                    );

                                    void        jacobi_mat_2d(
                                                                        int         np,
                                                                        cusfloat*   xn,
                                                                        cusfloat*   yn,
                                                                        cusfloat    xi,
                                                                        cusfloat    eta,
                                                                        cusfloat*   jac_mat
                                                            );

template<int NV, int NP>    inline  void    jacobi_mat_2d_vec(
                                                                        cusfloat*   xn,
                                                                        cusfloat*   yn,
                                                                        cusfloat*   xi,
                                                                        cusfloat*   eta,
                                                                        cusfloat*   jac_mat
                                                                );

                                    void        shape_fcn_2d( 
                                                                        int         np, 
                                                                        cusfloat    xi, 
                                                                        cusfloat    eta,
                                                                        cusfloat*   N
                                                            );

template<int NV, int NP>    inline  void        shape_fcn_2d_vec( 
                                                                        cusfloat*   xi,
                                                                        cusfloat*   eta,
                                                                        cusfloat*   N
                                                            );

                                    void        shape_fcn_der_2d(
                                                                        int         np,
                                                                        cusfloat    xi,
                                                                        cusfloat    eta,
                                                                        cusfloat*   dNdxi,
                                                                        cusfloat*   dNdeta
                                                            );

template<int NV, int NP>    inline  void        shape_fcn_der_2d_vec(
                                                                        cusfloat*   xi,
                                                                        cusfloat*   eta,
                                                                        cusfloat*   dNdxi,
                                                                        cusfloat*   dNdeta
                                                                    );

                                    void        shape_fcn_deta_2d( 
                                                                        int         np, 
                                                                        cusfloat    xi, 
                                                                        cusfloat*   N
                                                                );

template<int NV, int NP>    inline  void        shape_fcn_deta_2d_vec( 
                                                                        cusfloat*   xi, 
                                                                        cusfloat*   N
                                                                );

                                    void        shape_fcn_dxi_2d( 
                                                                        int         np, 
                                                                        cusfloat    eta,
                                                                        cusfloat*   N
                                                            );

template<int NV, int NP>    inline  void        shape_fcn_dxi_2d_vec( 
                                                                        cusfloat*   eta,
                                                                        cusfloat*   N
                                                                    );

                                    void        simplex_to_cartesian( 
                                                                        cusfloat    x, 
                                                                        cusfloat    y, 
                                                                        cusfloat&   xm,
                                                                        cusfloat&   ym
                                                                    );


#include "topology.txx"
