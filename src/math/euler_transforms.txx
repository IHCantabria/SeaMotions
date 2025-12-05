
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

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "euler_transforms.hpp"
#include "math_interface.hpp"
#include "math_tools.hpp"


template<typename T>
void    euler_local_to_global(
                                        T           alpha,
                                        T           beta,
                                        T           gamma,
                                        T*          rot_mat
                            )
{

    // Get rotation matrixes for X, Y, Z axes
    cusfloat    alpha_f = 1.0;
    cusfloat    beta_f  = 1.0;
    const int   rows_np = 3;
    const int   rotm_np = pow2s( rows_np );
    T           rotm_a[rotm_np]; clear_vector( rotm_np, rotm_a );
    T           rotm_x[rotm_np]; clear_vector( rotm_np, rotm_x );
    T           rotm_y[rotm_np]; clear_vector( rotm_np, rotm_y );
    T           rotm_z[rotm_np]; clear_vector( rotm_np, rotm_z );

    _rot_x( rotm_x, alpha );
    _rot_y( rotm_y, beta );
    _rot_z( rotm_z, gamma );

    // Calculate Ry · Rx
    cblas_gemm<T>(
                            CblasRowMajor, 
                            CblasNoTrans, 
                            CblasNoTrans, 
                            rows_np, 
                            rows_np, 
                            rows_np, 
                            &alpha_f, 
                            rotm_y, 
                            rows_np, 
                            rotm_x, 
                            rows_np, 
                            &beta_f, 
                            rotm_a, 
                            rows_np
                        );

    // Calcualte Rz · ( Ry · Rx )
    cblas_gemm<T>(
                            CblasRowMajor, 
                            CblasNoTrans, 
                            CblasNoTrans, 
                            rows_np, 
                            rows_np, 
                            rows_np, 
                            &alpha_f, 
                            rotm_z, 
                            rows_np, 
                            rotm_a, 
                            rows_np, 
                            &beta_f, 
                            rot_mat, 
                            rows_np
                        );

}


template<typename T>
void    euler_local_to_global_disp(
                                        T*          dofs_trans,
                                        T*          dofs_rot,
                                        cusfloat*   radius,
                                        T*          displacement
                                    )
{
    // Clear displacement vector just in clase
    clear_vector( 3, displacement );
    
    // Get local to global transformation matrix
    cusfloat    alpha   = 1.0;
    cusfloat    beta    = 1.0;
    int         icnx    = 1;
    int         icny    = 1;
    const int   rows_np = 3;
    const int   rotm_np = pow2s( rows_np );
    T*          rot_mat = generate_empty_vector<T>( rotm_np );

    euler_local_to_global(
                            dofs_rot[0],
                            dofs_rot[1],
                            dofs_rot[2],
                            rot_mat
                        );

    // Get local radius value in global coorindates
    cblas_gemv<T>(
                            CblasRowMajor,
                            CblasNoTrans,
                            rows_np,
                            rows_np,
                            &alpha,
                            rot_mat,
                            rows_np,
                            radius,
                            icnx,
                            &beta,
                            displacement,
                            icny
                        );

    // Add translation displacements in order to get the 
    // point displacement in global coorindates
    for ( int i=0; i<3; i++ )
    {
        displacement[i] += dofs_trans[i] - radius[i];
    }

    // Delete local heap memory
    mkl_free( rot_mat );
}


template<typename T>
void    _rot_x(
                                        T*          mat,
                                        T           alpha
                )
{
    mat[0] = 1.0;
    mat[1] = 0.0;
    mat[2] = 0.0;

    mat[3] = 0.0;
    mat[4] = std::cos( alpha );
    mat[5] = -std::sin( alpha );
    
    mat[6] = 0.0;
    mat[7] = std::sin( alpha );
    mat[8] = std::cos( alpha );
}


template<typename T>
void    _rot_y(
                                        T*          mat,
                                        T           beta
                )
{
    mat[0] = std::cos( beta );
    mat[1] = 0.0;
    mat[2] = std::sin( beta );

    mat[3] = 0.0;
    mat[4] = 1.0;
    mat[5] = 0.0;

    mat[6] = -std::sin( beta );
    mat[7] = 0.0;
    mat[8] = std::cos( beta );
}


template<typename T>
void    _rot_z(
                                        T*          mat,
                                        T           gamma
                )
{
    mat[0] = std::cos( gamma );
    mat[1] = -std::sin( gamma );
    mat[2] = 0.0;

    mat[3] = std::sin( gamma );
    mat[4] = std::cos( gamma );
    mat[5] = 0.0;

    mat[6] = 0.0;
    mat[7] = 0.0;
    mat[8] = 1.0;
}