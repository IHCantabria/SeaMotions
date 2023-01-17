
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_interface.hpp"
#include "../../src/tools.hpp"


template<typename T>
void gemv(int N, T* mat, T* vec, T* sol)
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            sol[i] += mat[i*N+j]*vec[j];
        }
    }
}


int main(void)
{
    int num_reps = 1e5;
    cusfloat mat[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    cusfloat vec[3] = {1.0, 2.0, 3.0};
    cusfloat sol[3] = {0.0, 0.0, 0.0};

    double t0 = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        cblas_gemv<cusfloat>(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, mat, 3, vec, 1, 0, sol, 1);
    }
    double t1 = get_cpu_time();
    double blas_time = t1-t0;

    double t2 = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        gemv(3, mat, vec, sol);
    }
    double t3 = get_cpu_time();
    double sea_time = t3-t2;

    std::cout << "Blas time [s]: " << blas_time << std::endl;
    std::cout << "Sea time [s]: " << sea_time << std::endl;

    return 0;
}