
#include <cmath>
#include <iostream>
#include "mkl.h"
#include "../../src/tools.hpp"


template<typename T>
void cvAdd(int N, T* u, T* v, T* w)
{
    for(int i=0; i<N; i++)
    {
        w[i] = u[i] + v[i];
    }
}


template<typename T>
T nrm2(int N, T* v)
{
    T a = 0;
    for (int i=0; i<N; i++)
    {
        a += pow(v[i], 2.0);
    }

    return sqrt(a);
}


int main(void)
{
    // Create vectors to sum up
    int num_reps = 1e6;
    const int N = 40;
    double u[N], v[N], w[N];
    for (int i=0; i<N; i++)
    {
        u[i] = i;
        v[i] = 2*i;
    }

    double t2 = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        // vdAdd(N, u, v, w);
        cblas_dnrm2(N, u, 1);
    }
    double t3 = get_cpu_time();
    double vm_time = (t3-t2);

    double t0 = get_cpu_time();
    for (int i=0; i<num_reps; i++)
    {
        // cvAdd(N, u, v, w);
        nrm2(N, u);
    }
    double t1 = get_cpu_time();
    double compiler_time = (t1-t0);


    std::cout << "Compiler Time [s]: " << compiler_time << std::endl;
    std::cout << "VM Time [s]: " << vm_time << std::endl;


    return 0;
}