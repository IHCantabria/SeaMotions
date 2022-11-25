
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include <cmath>


template<typename T>
inline int assert_scalar_equality(T &u, T &v, T epsilon)
{
    int pass = 1;
    if (std::abs(u-v)>epsilon)
    {
        pass = 0;
    }

    return pass;
}


template<typename T>
inline int assert_vector_equality(int N, T* u, T* v, T epsilon)
{
    int pass = 1;
    for (int i=0; i<N; i++)
    {
        if (std::abs(u[i]-v[i])>epsilon)
        {
            pass = 0;
        }
    }

    return pass;
}


template<>
inline int assert_vector_equality<int>(int N, int* u, int* v, int epsilon)
{
    int pass = 1;
    for (int i=0; i<N; i++)
    {
        if ((u[i]-v[i]) != epsilon)
        {
            pass = 0;
        }
    }

    return pass;
}


template<typename T>
inline void cross(T (&u)[3], T (&v)[3], T (&w)[3])
{
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}


template<typename T>
inline T pow2s(T x)
{
    return x*x;
}


template<typename T>
inline T pow3s(T x)
{
    return x*x*x;
}


template<typename T>
inline void sv_add(int n, T* u, T* v, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i] + v[i];
    }
}


template<typename T>
inline void svs_add(int n, T* u, T s, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = s + u[i];
    }
}


template<typename T>
inline void sv_div(int n, T* u, T* v, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i]/v[i];
    }
}


template<typename T>
inline void sv_inv(int n, T s, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = s/u[i];
    }
}


template<typename T>
inline void sv_mod(int n, T* u, T &mod)
{
    mod = 0.0;
    for (int i=0; i<n; i++)
    {
        mod += u[i]*u[i];
    }
    mod = std::sqrt(mod);
}


template<typename T>
void sv_mult(int n, T* u, T* v, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i] * v[i];
    }
}


template<typename T>
inline void svs_mult(int n, T* u, T s, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = s*u[i];
    }
}


template<typename T>
inline void sv_pow2(int n, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i]*u[i];
    }
}


template<typename T>
inline void sv_sqrt(int n, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = std::sqrt(u[i]);
    }
}


template<typename T>
inline void sv_sub(int n, T* u, T* v, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i] - v[i];
    }
}


template<typename T>
inline void svs_sub(int n, T* u, T s, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i] - s;
    }
}


template<typename T>
inline void svs_sub(int n, T s, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = s - u[i];
    }
}


template<typename T>
inline void print_vector(int n, T* v, int mode)
{
    std::cout << "v[0]: " << v[0];
    if (mode == 0)
    {
        for (int i=1; i<n; i++)
        {
            std::cout << " - v[" << i <<"]: " << v[i];
        }
    }
    else if (mode == 1)
    {
        for (int i=1; i<n; i++)
        {
            std::cout << "\nv[" << i <<"]: " << v[i];
        }
    }
    else
    {
        throw std::runtime_error("Mode specified is not correct. It must be 0: horizontal | 1: vertical.");
    }
    std::cout << std::endl;
    
}