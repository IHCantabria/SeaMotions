
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
#include "topology.hpp"


template<int NV, int NP>
inline  void    jacobi_det_2d_vec(
                                    cusfloat*   xn,
                                    cusfloat*   yn,
                                    cusfloat*   xi,
                                    cusfloat*   eta,
                                    cusfloat*   jac_det
                                )
{
    /**
     * @brief Calculate jacobi determinant for 2D coordinate systems and using the shape 
     *        functions to represent the topological space. The jacobi determinant
     *        is calculated at a given point xi, eta.
     * 
     * @param xn        IN: First coordinate spatial points
     * @param yn        IN: Second coordinate spatial points
     * @param xi        IN: First coordinate evaluation points
     * @param eta       IN: Second coordinate evaluation points
     * 
     */

    // Calculate jacobi matrix
    static cusfloat jac_mat[4*NP];
    jacobi_mat_2d_vec<NV,NP>( xn, yn, xi, eta, jac_mat );

    // Calculate the determinant
    for ( int i=0; i<NP; i++ )
    {
        jac_det[i] = jac_mat[4*i+0] * jac_mat[4*i+3] - jac_mat[4*i+1] * jac_mat[4*i+2];
    }

}


template<int NV, int NP>
inline  void    jacobi_mat_2d_vec(
                                    cusfloat*   xn,
                                    cusfloat*   yn,
                                    cusfloat*   xi,
                                    cusfloat*   eta,
                                    cusfloat*   jac_mat
                                )
{
    /**
     * @brief Calculate jacobi matrix for 2D coordinate systems and using the shape 
     *        functions to represent the topological space. The jacobi matrix
     *        is calculated at a given point xi, eta.
     * 
     * @param xn        IN: First coordinate spatial points
     * @param yn        IN: Second coordinate spatial points
     * @param xi        IN: First coordinate evaluation point
     * @param eta       IN: Second coordinate evaluation point
     * @param jac_mat   OUT: Jacobi matrix outputed
     * 
     */
    // Get shape functions derivative values for the required points
    static cusfloat  dNdxi[NV*NP]  = { };
    static cusfloat  dNdeta[NV*NP] = { };
    
    shape_fcn_der_2d_vec<NV, NP>( xi, eta, dNdxi, dNdeta );

    // Calculate jacobi matrix entries    
    JACMAT2D_DOT_PRODUCT( jac_mat, 0, xn, dNdxi, NP, NV )
    JACMAT2D_DOT_PRODUCT( jac_mat, 1, xn, dNdeta, NP, NV )
    JACMAT2D_DOT_PRODUCT( jac_mat, 2, yn, dNdxi, NP, NV )
    JACMAT2D_DOT_PRODUCT( jac_mat, 3, yn, dNdeta, NP, NV )

}


template<int NV, int NP>
inline  void    shape_fcn_2d_vec( 
                                    cusfloat*   xi,
                                    cusfloat*   eta,
                                    cusfloat*   N
                                )
{
    /**
     * @brief Calculate value of the 2d shape functions for a given point in the local
     *        frame of reference
     * 
     * @param xi    IN: First dimension value
     * @param eta   IN: Second dimension value
     * @param N     OUT: Shape functions value.
     * 
     */

    static_assert(NV == 3 || NV == 4, "Only triangle and quadrilateral panels supported.");

    // Switch in between elements type
    if constexpr( NV == 3 )
    {
        for ( int i=0; i<NP; i++ )
        {
            N[3*i+0] = ( 1 - xi[i] ) * ( 1 - eta[i] ) / 4.0;
            N[3*i+1] = ( 1 + xi[i] ) * ( 1 - eta[i] ) / 4.0;
            N[3*i+2] = ( 1 + eta[i] ) / 2.0;
        }
    }

    if constexpr( NV == 4 )
    {
        for ( int i=0; i<NP; i++ )
        {
            N[4*i+0] = ( 1 - xi[i] ) * ( 1 - eta[i] ) / 4.0;
            N[4*i+1] = ( 1 + xi[i] ) * ( 1 - eta[i] ) / 4.0;
            N[4*i+2] = ( 1 + xi[i] ) * ( 1 + eta[i] ) / 4.0;
            N[4*i+3] = ( 1 - xi[i] ) * ( 1 + eta[i] ) / 4.0;
        }
    }
}



template<int NV, int NP>
inline  void    shape_fcn_der_2d_vec(
                                        cusfloat*   xi,
                                        cusfloat*   eta,
                                        cusfloat*   dNdxi,
                                        cusfloat*   dNdeta
                                    )
{
    /**
     * @brief Calculate value of the 2d shape functions second dimension derivative
     *        for a given point in the local frame of reference
     * 
     * @param xi        IN: First dimension value
     * @param eta       IN: Second dimension value
     * @param dNdxi     OUT: Shape functions derivative value w.r.t first dimension.
     * @param dNdeta    OUT: Shape functions derivative value w.r.t second dimension.
     * 
     */

    // Calculate derivative with respect to the first dimension
    shape_fcn_dxi_2d_vec<NV, NP>( eta, dNdxi );

    // Calculate derivative with respect to the first dimension
    shape_fcn_deta_2d_vec<NV, NP>( xi, dNdeta );
    
}


template<int NV, int NP>
inline  void    shape_fcn_deta_2d_vec( 
                                        cusfloat*   xi, 
                                        cusfloat*   N
                                    )
{
    /**
     * @brief Calculate value of the 2d shape functions second dimension derivative
     *        for a given point in the local frame of reference
     * 
     * @param xi    IN: First dimension value
     * @param N     OUT: Shape functions derivative value.
     * 
     */

    // Switch in between elements type
    if constexpr( NV == 3 )
    {
        for ( int i=0; i<NP; i++ )
        {
            N[3*i+0] = -( 1 - xi[i] ) / 4.0;
            N[3*i+1] = -( 1 + xi[i] ) / 4.0;
            N[3*i+2] =  0.5;
        }
    }

    // Calculate shape function values for rectangular regions
    if constexpr( NV == 4 )
    {
        for ( int i=0; i<NP; i++ )
        {
            N[4*i+0] = -( 1 - xi[i] ) / 4.0;
            N[4*i+1] = -( 1 + xi[i] ) / 4.0;
            N[4*i+2] =  ( 1 + xi[i] ) / 4.0;
            N[4*i+3] =  ( 1 - xi[i] ) / 4.0;
        }
    }
}


template<int NV, int NP>
inline  void    shape_fcn_dxi_2d_vec( 
                                        cusfloat*   eta,
                                        cusfloat*   N
                                    )
{
    /**
     * @brief Calculate value of the 2d shape functions first dimension derivative
     *        for a given point in the local frame of reference
     * 
     * @param eta   IN: Second dimension value
     * @param N     OUT: Shape functions derivative value.
     * 
     */

    // Calculate shape function values for simplex regions
    if constexpr( NV == 3 )
    {
        for ( int i=0; i<NP; i++ )
        {
            N[3*i+0] = -( 1 - eta[i] ) / 4.0;
            N[3*i+1] =  ( 1 - eta[i] ) / 4.0;
            N[3*i+2] = 0.0;
        }
    }

    // Calculate shape function values for rectangular regions
    if constexpr( NV == 4 )
    {
        for ( int i=0; i<NP; i++ )
        {
            N[4*i+0] = -( 1 - eta[i] ) / 4.0;
            N[4*i+1] =  ( 1 - eta[i] ) / 4.0;
            N[4*i+2] =  ( 1 + eta[i] ) / 4.0;
            N[4*i+3] = -( 1 + eta[i] ) / 4.0;
        }
    }
}