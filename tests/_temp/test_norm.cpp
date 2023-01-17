
#include <iostream>
#include "mkl.h"

#include "../../src/config.hpp"
#include "../../src/math/math_interface.hpp"


void print_vector(cusfloat (&v)[3])
{
    std::cout << "v[0]: " << v[0];
    for (int i=1; i<3; i++)
    {
        std::cout << " - v[" << i <<"]: " << v[i];
    }
    std::cout << std::endl;
}


void print_vector_2(cusfloat* v, int num_points)
{
    std::cout << "v[0]: " << v[0];
    for (int i=1; i<num_points; i++)
    {
        std::cout << " - v[" << i <<"]: " << v[i];
    }
    std::cout << std::endl;
}


int main(void)
{
    cusfloat u[3] = {1.0, 2.0, 3.0};
    cusfloat v[3] = {0.0, 5.0, 0.0};
    cusfloat w[3];
    std::cout << "Pass by reference..." << std::endl;
    print_vector(v);
    std::cout << "Pass by pointer" << std::endl;
    print_vector_2(v, 3);

    MKL_INT N = 3;
    std::cout << cblas_ddot(N, v, 1, v, 1) << std::endl;
    std::cout << cblas_dnrm2(N, v, 1) << std::endl;
    std::cout << cblas_nrm2<cusfloat>(N, v, 1) << std::endl;
    
    vdAdd(3, u, v, w);
    print_vector(w);


    return 0;
}