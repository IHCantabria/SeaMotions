
// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "euler_rotations.hpp"
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
    const int   rotm_np = 9;
    cusfloat    rotm_a[rotm_np]; clear_vector( rotm_np, rotm_a );
    cusfloat    rotm_x[rotm_np]; clear_vector( rotm_np, rotm_x );
    cusfloat    rotm_y[rotm_np]; clear_vector( rotm_np, rotm_y );
    cusfloat    rotm_z[rotm_np]; clear_vector( rotm_np, rotm_z );

    _rot_x( rotm_x, alpha );
    _rot_y( rotm_y, alpha );
    _rot_z( rotm_z, alpha );

    // Calculate Ry · Rx
    cblas_gemm<cusfloat>(
                            CblasRowMajor, 
                            CblasNoTrans, 
                            CblasNoTrans, 
                            rotm_np, 
                            rotm_np, 
                            rotm_np, 
                            1.0, 
                            rotm_y, 
                            rotm_np, 
                            rotm_x, 
                            rotm_np, 
                            1.0, 
                            rotm_a, 
                            rotm_np
                        );

    // Calcualte Rz · ( Ry · Rx )
    cblas_gemm<cusfloat>(
                            CblasRowMajor, 
                            CblasNoTrans, 
                            CblasNoTrans, 
                            rotm_np, 
                            rotm_np, 
                            rotm_np, 
                            1.0, 
                            rotm_z, 
                            rotm_np, 
                            rotm_a, 
                            rotm_np, 
                            1.0, 
                            rot_mat, 
                            rotm_np
                        );

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