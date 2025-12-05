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
#include "gauss_t.hpp"
#include "topology.hpp"
#include "topology_gauss.hpp"


template<int NV, int NGP>
inline  void    jacobi_det_2d_gp(
                                    cusfloat*   xn,
                                    cusfloat*   yn,
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
    static cusfloat jac_mat[4*NGP*NGP];
    jacobi_mat_2d_gp<NV,NGP>( xn, yn, jac_mat );

    // Calculate the determinant
    for ( int i=0; i<NGP*NGP; i++ )
    {
        jac_det[i] = jac_mat[4*i+0] * jac_mat[4*i+3] - jac_mat[4*i+1] * jac_mat[4*i+2];
    }

}


template<int NV, int NGP>
inline  void    jacobi_mat_2d_gp(
                                    cusfloat*   xn,
                                    cusfloat*   yn,
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
    static cusfloat  dNdxi[NV*NGP*NGP]  = { };
    static cusfloat  dNdeta[NV*NGP*NGP] = { };
    
    shape_fcn_der_2d_gp<NV, NGP>( dNdxi, dNdeta );

    // Calculate jacobi matrix entries    
    JACMAT2D_DOT_PRODUCT( jac_mat[0], xn, dNdxi, NGP*NGP, NV )
    JACMAT2D_DOT_PRODUCT( jac_mat[1], xn, dNdeta, NGP*NGP, NV )
    JACMAT2D_DOT_PRODUCT( jac_mat[2], yn, dNdxi, NGP*NGP, NV )
    JACMAT2D_DOT_PRODUCT( jac_mat[3], yn, dNdeta, NGP*NGP, NV )

}


template<int NV, int NGP>
inline  void    shape_fcn_2d_gp( 
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
    if ( constexpr( NV == 3 ) )
    {
        for ( int i=0; i<NGP*NGP; i++ )
        {
            N[3*i+0] = ( 1 - GaussPointsT<1,NGP>::roots_x[i] ) * ( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[3*i+1] = ( 1 + GaussPointsT<1,NGP>::roots_x[i] ) * ( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[3*i+2] = ( 1 + GaussPointsT<1,NGP>::roots_y[i] ) / 2.0;
        }
    }

    if ( constexpr( NV == 4 ) )
    {
        for ( int i=0; i<NGP*NGP; i++ )
        {
            N[4*i+0] = ( 1 - GaussPointsT<1,NGP>::roots_x[i] ) * ( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[4*i+1] = ( 1 + GaussPointsT<1,NGP>::roots_x[i] ) * ( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[4*i+2] = ( 1 + GaussPointsT<1,NGP>::roots_x[i] ) * ( 1 + GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[4*i+3] = ( 1 - GaussPointsT<1,NGP>::roots_x[i] ) * ( 1 + GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
        }
    }
}



template<int NV, int NGP>
inline  void    shape_fcn_der_2d_gp(
                                        cusfloat*   dNdxi,
                                        cusfloat*   dNdeta
                                    )
{
    /**
     * @brief Calculate value of the 2D shape functions' derivatives 
     *        for a given point in the local frame of reference.
     *
     * This function evaluates the partial derivatives of shape functions 
     * with respect to the local coordinates (`xi`, `eta`) at a Gauss point.
     * The template parameters allow compile-time specification of the number 
     * of vertices and Gauss points.
     * 
     * @tparam NV        Number of vertices of the element.
     * @tparam NGP       Number of Gauss points.
     * 
     * @param dNdxi      OUT: Derivative of shape functions w.r.t. ξ.
     * @param dNdeta     OUT: Derivative of shape functions w.r.t. η.
     */

    // Calculate derivative with respect to the first dimension
    shape_fcn_dxi_2d_gp<NV, NGP>( dNdxi );

    // Calculate derivative with respect to the first dimension
    shape_fcn_deta_2d_gp<NV, NGP>( dNdeta );
    
}


template<int NV, int NGP>
inline  void    shape_fcn_deta_2d_gp( 
                                        cusfloat*   N
                                    )
{
    /**
     * @brief Calculate value of the 2d shape functions second dimension derivative
     *        for a given number of gauss points in the local frame of reference
     * 
     * @param N     OUT: Shape functions derivative value.
     * 
     */

    // Switch in between elements type
    if ( constexpr( NV == 3 ) )
    {
        for ( int i=0; i<NGP*NGP; i++ )
        {
            N[3*i+0] = -( 1 - GaussPointsT<1,NGP>::roots_x[i] ) / 4.0;
            N[3*i+1] = -( 1 + GaussPointsT<1,NGP>::roots_x[i] ) / 4.0;
            N[3*i+2] =  0.5;
        }
    }

    // Calculate shape function values for rectangular regions
    if ( constexpr( NV == 4 ) )
    {
        for ( int i=0; i<NGP*NGP; i++ )
        {
            N[4*i+0] = -( 1 - GaussPointsT<1,NGP>::roots_x[i] ) / 4.0;
            N[4*i+1] = -( 1 + GaussPointsT<1,NGP>::roots_x[i] ) / 4.0;
            N[4*i+2] =  ( 1 + GaussPointsT<1,NGP>::roots_x[i] ) / 4.0;
            N[4*i+3] =  ( 1 - GaussPointsT<1,NGP>::roots_x[i] ) / 4.0;
        }
    }
}


template<int NV, int NGP>
inline  void    shape_fcn_dxi_2d_gp( 
                                        cusfloat*   N
                                    )
{
    /**
     * @brief Calculate value of the 2d shape functions first dimension derivative
     *        for a given number of gauss points in the local frame of reference
     * 
     * @param N     OUT: Shape functions derivative value.
     * 
     */

    // Calculate shape function values for simplex regions
    if ( constexpr( NV == 3 ) )
    {
        for ( int i=0; i<NGP*NGP; i++ )
        {
            N[3*i+0] = -( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[3*i+1] =  ( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[3*i+2] = 0.0;
        }
    }

    // Calculate shape function values for rectangular regions
    if ( constexpr( NV == 4 ) )
    {
        for ( int i=0; i<NGP*NGP; i++ )
        {
            N[4*i+0] = -( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[4*i+1] =  ( 1 - GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[4*i+2] =  ( 1 + GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
            N[4*i+3] = -( 1 + GaussPointsT<1,NGP>::roots_y[i] ) / 4.0;
        }
    }
}