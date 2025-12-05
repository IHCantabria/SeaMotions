
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
#include "math_interface.hpp"
#include "math_tools.hpp"
#include "topology.hpp"


cusfloat    jacobi_det_2d(
                            int         np,
                            cusfloat*   xn,
                            cusfloat*   yn,
                            cusfloat    xi,
                            cusfloat    eta
                        )
{
    /**
     * @brief Calculate jacobi determinant for 2D coordinate systems and using the shape 
     *        functions to represent the topological space. The jacobi determinant
     *        is calculated at a given point xi, eta.
     * 
     * @param np        IN: Number of shape functions
     * @param xn        IN: First coordinate spatial points
     * @param yn        IN: Second coordinate spatial points
     * @param xi        IN: First coordinate evaluation point
     * @param eta       IN: Second coordinate evaluation point
     * 
     * @return              Jacobi matrix outputed
     * 
     */

    // Calculate jacobi matrix
    cusfloat jac_mat[4] = { 0.0, 0.0, 0.0, 0.0 };
    cusfloat jac_det    = 0.0;

    jacobi_mat_2d(  
                    np,
                    xn,
                    yn,
                    xi,
                    eta,
                    jac_mat
                );

    // Calculate the determinant
    jac_det = jac_mat[0] * jac_mat[3] - jac_mat[1] * jac_mat[2];

    return jac_det;
}


void jacobi_mat_2d(
                        int         np,
                        cusfloat*   xn,
                        cusfloat*   yn,
                        cusfloat    xi,
                        cusfloat    eta,
                        cusfloat*   jac_mat
                    )
{
    /**
     * @brief Calculate jacobi matrix for 2D coordinate systems and using the shape 
     *        functions to represent the topological space. The jacobi matrix
     *        is calculated at a given point xi, eta.
     * 
     * @param np        IN: Number of shape functions
     * @param xn        IN: First coordinate spatial points
     * @param yn        IN: Second coordinate spatial points
     * @param xi        IN: First coordinate evaluation point
     * @param eta       IN: Second coordinate evaluation point
     * @param jac_mat   OUT: Jacobi matrix outputed
     * 
     */
    // Set integer constants to use blas functions
    const int npc   = np;
    const int icnx  = 1;

    // Get shape functions derivative values for the required points
    cusfloat* dNdxi  = generate_empty_vector<cusfloat>( np );
    cusfloat* dNdeta = generate_empty_vector<cusfloat>( np );
    
    shape_fcn_der_2d( np, xi, eta, dNdxi, dNdeta );

    // Calculate jacobi matrix entries
    jac_mat[0] = cblas_dot<cusfloat>( npc, xn, icnx, dNdxi, icnx );
    jac_mat[1] = cblas_dot<cusfloat>( npc, xn, icnx, dNdeta, icnx );
    jac_mat[2] = cblas_dot<cusfloat>( npc, yn, icnx, dNdxi, icnx );
    jac_mat[3] = cblas_dot<cusfloat>( npc, yn, icnx, dNdeta, icnx );

    // Deallocate heap memory
    mkl_free( dNdxi );
    mkl_free( dNdeta );

}


void shape_fcn_2d( 
                    int         np, 
                    cusfloat    xi, 
                    cusfloat    eta,
                    cusfloat*   N
                )
{
    /**
     * @brief Calculate value of the 2d shape functions for a given point in the local
     *        frame of reference
     * 
     * @param np    IN: Number of vertexes of the polygon
     * @param xi    IN: First dimension value
     * @param eta   IN: Second dimension value
     * @param N     OUT: Shape functions value.
     * 
     */

    // Switch in between elements type
    switch ( np )
    {
        case 3:
            // Calculate shape function values for simplex regions
            N[0] = ( 1 - xi ) * ( 1 - eta ) / 4.0;
            N[1] = ( 1 + xi ) * ( 1 - eta ) / 4.0;
            N[2] = ( 1 + eta ) / 2.0;
            break;

        case 4:
            // Calculate shape function values for rectangular regions
            N[0] = ( 1 - xi ) * ( 1 - eta ) / 4.0;
            N[1] = ( 1 + xi ) * ( 1 - eta ) / 4.0;
            N[2] = ( 1 + xi ) * ( 1 + eta ) / 4.0;
            N[3] = ( 1 - xi ) * ( 1 + eta ) / 4.0;
            break;
        
        default:
            break;
    }
}


void shape_fcn_der_2d(
                        int         np,
                        cusfloat    xi,
                        cusfloat    eta,
                        cusfloat*   dNdxi,
                        cusfloat*   dNdeta
                    )
{
    /**
     * @brief Calculate value of the 2d shape functions second dimension derivative
     *        for a given point in the local frame of reference
     * 
     * @param np        IN: Number of vertexes of the polygon
     * @param xi        IN: First dimension value
     * @param eta       IN: Second dimension value
     * @param dNdxi     OUT: Shape functions derivative value w.r.t first dimension.
     * @param dNdeta    OUT: Shape functions derivative value w.r.t second dimension.
     * 
     */

    // Calculate derivative with respect to the first dimension
    shape_fcn_dxi_2d( np, eta, dNdxi );

    // Calculate derivative with respect to the first dimension
    shape_fcn_deta_2d( np, xi, dNdeta );
    
}


void shape_fcn_deta_2d( 
                            int         np, 
                            cusfloat    xi, 
                            cusfloat*   N
                        )
{
    /**
     * @brief Calculate value of the 2d shape functions second dimension derivative
     *        for a given point in the local frame of reference
     * 
     * @param np    IN: Number of vertexes of the polygon
     * @param xi    IN: First dimension value
     * @param N     OUT: Shape functions derivative value.
     * 
     */

    // Switch in between elements type
    switch ( np )
    {
    case 3:
        // Calculate shape function values for simplex regions
        N[0] = -( 1 - xi ) / 4.0;
        N[1] = -( 1 + xi ) / 4.0;
        N[2] =  0.5;
        break;

    case 4:
        // Calculate shape function values for rectangular regions
        N[0] = -( 1 - xi ) / 4.0;
        N[1] = -( 1 + xi ) / 4.0;
        N[2] =  ( 1 + xi ) / 4.0;
        N[3] =  ( 1 - xi ) / 4.0;
    
    default:
        break;
    }
}


void shape_fcn_dxi_2d( 
                        int         np, 
                        cusfloat    eta,
                        cusfloat*   N
                    )
{
    /**
     * @brief Calculate value of the 2d shape functions first dimension derivative
     *        for a given point in the local frame of reference
     * 
     * @param np    IN: Number of vertexes of the polygon
     * @param eta   IN: Second dimension value
     * @param N     OUT: Shape functions derivative value.
     * 
     */

    // Switch in between elements type
    switch ( np )
    {
    case 3:
        // Calculate shape function values for simplex regions
        N[0] = -( 1 - eta ) / 4.0;
        N[1] =  ( 1 - eta ) / 4.0;
        N[2] = 0.0;
        break;

    case 4:
        // Calculate shape function values for rectangular regions
        N[0] = -( 1 - eta ) / 4.0;
        N[1] =  ( 1 - eta ) / 4.0;
        N[2] =  ( 1 + eta ) / 4.0;
        N[3] = -( 1 + eta ) / 4.0;
    
    default:
        break;
    }
}


void simplex_to_cartesian( 
                            cusfloat    x, 
                            cusfloat    y, 
                            cusfloat&   xm,
                            cusfloat&   ym
                        )
{
    /**
     * @brief Convet from collapsed coordinate system to cartesian coordinates
     * 
     * @param x     IN: First dimension in collapsed coordinate system
     * @param y     IN: Second dimension in collapsed coordinate system
     * @param xm    OUT: First dimension in cartesian coordinate system
     * @param ym    OUT: Second dimensino in cartesian coordinate system
     * 
     */

    xm = ( 1 + x ) * ( 1 - y ) / 2.0 - 1.0;
    ym = y;

}