
// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "euler_transforms.hpp"
#include "math_interface.hpp"
#include "math_tools.hpp"


void    euler_local_to_global(
                                        cusfloat    alpha,
                                        cusfloat    beta,
                                        cusfloat    gamma,
                                        cusfloat*   rot_mat
                            )
{

    // Get rotation matrixes for X, Y, Z axes
    const int   rows_np = 3;
    const int   rotm_np = pow2s( rows_np );
    cusfloat    rotm_a[rotm_np]; clear_vector( rotm_np, rotm_a );
    cusfloat    rotm_x[rotm_np]; clear_vector( rotm_np, rotm_x );
    cusfloat    rotm_y[rotm_np]; clear_vector( rotm_np, rotm_y );
    cusfloat    rotm_z[rotm_np]; clear_vector( rotm_np, rotm_z );

    _rot_x( rotm_x, alpha );
    _rot_y( rotm_y, beta );
    _rot_z( rotm_z, gamma );

    // Calculate Ry · Rx
    cblas_gemm<cusfloat>(
                            CblasRowMajor, 
                            CblasNoTrans, 
                            CblasNoTrans, 
                            rows_np, 
                            rows_np, 
                            rows_np, 
                            1.0, 
                            rotm_y, 
                            rows_np, 
                            rotm_x, 
                            rows_np, 
                            1.0, 
                            rotm_a, 
                            rows_np
                        );

    // Calcualte Rz · ( Ry · Rx )
    cblas_gemm<cusfloat>(
                            CblasRowMajor, 
                            CblasNoTrans, 
                            CblasNoTrans, 
                            rows_np, 
                            rows_np, 
                            rows_np, 
                            1.0, 
                            rotm_z, 
                            rows_np, 
                            rotm_a, 
                            rows_np, 
                            1.0, 
                            rot_mat, 
                            rows_np
                        );

}


void    euler_local_to_global_disp(
                                        cusfloat*   dofs_trans,
                                        cusfloat*   dofs_rot,
                                        cusfloat*   radius,
                                        cusfloat*   displacement
                                    )
{
    // Clear displacement vector just in clase
    clear_vector( 3, displacement );
    
    // Get local to global transformation matrix
    const int   rows_np = 3;
    const int   rotm_np = pow2s( rows_np );
    cusfloat*   rot_mat = generate_empty_vector<cusfloat>( rotm_np );

    euler_local_to_global(
                            dofs_rot[0],
                            dofs_rot[1],
                            dofs_rot[2],
                            rot_mat
                        );

    // Get local radius value in global coorindates
    cblas_gemv<cusfloat>(
                            CblasRowMajor,
                            CblasNoTrans,
                            rows_np,
                            rows_np,
                            1.0,
                            rot_mat,
                            rows_np,
                            radius,
                            1,
                            1.0,
                            displacement,
                            1
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


void    _rot_x(
                                        cusfloat*   mat,
                                        cusfloat    alpha
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


void    _rot_y(
                                        cusfloat*   mat,
                                        cusfloat    beta
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


void    _rot_z(
                                        cusfloat*   mat,
                                        cusfloat    gamma
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