
// Include general usage libraries
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "tools.hpp"


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
inline void print_matrix(int num_rows, int num_cols, T* mat, int precision, int align, int scient_not)
{
    // Look for the biggest number in the matrix in order to calculate the width of the columns
    T max = mat[0];
    for (int i=0; i<num_rows*num_cols; i++)
    {
        if (mat[i]>max)
        {
            max = mat[i];
        }
    }

    size_t col_width = 0;
    if (max > 1.0)
    {
        max = std::floor(std::log10(max));
        col_width = static_cast<int>(max)+precision+6;
    }
    else
    {
        col_width = static_cast<int>(precision+6);
    }

    if (scient_not == 1)
    {
        col_width += 4;
    }
    else if (scient_not != 0)
    {
        std::string err_message("Scient notation flag must be: 0 - No Scientist notation | 1: Scientist notation.");
        std::cerr << err_message << std::endl;
        throw std::runtime_error(err_message);
    }

    size_t row_header_width = static_cast<int>(std::floor(std::log10(num_rows))+10);

    // Print table data
    std::stringstream ss_row_header;
    for (int i=0; i<num_rows; i++)
    {
        ss_row_header.str("");
        ss_row_header << "R [" << i << "]: ";
        std::cout << align_str(ss_row_header.str(), row_header_width, 1) << align_num(mat[num_cols*i], col_width, precision, align, scient_not);

        //  (mat[num_cols*i], col_width, precision, scient_not);
        for (int j=1; j<num_cols; j++)
        {
            std::cout << " | " << align_num(mat[num_cols*i+j], col_width, precision, align, scient_not);
        }
        std::cout << std::endl;
    }
}


template<typename T>
inline void print_vector(int n, T* v, int mode, int precision)
{
    std::cout << "v[0]: " << v[0];
    if (mode == 0)
    {
        for (int i=1; i<n; i++)
        {
            std::cout << " - v[" << i <<"]: " << std::setprecision(precision) << v[i];
        }
    }
    else if (mode == 1)
    {
        for (int i=1; i<n; i++)
        {
            std::cout << "\nv[" << i <<"]: " << std::setprecision(precision) << v[i];
        }
    }
    else
    {
        throw std::runtime_error("Mode specified is not correct. It must be 0: horizontal | 1: vertical.");
    }
    std::cout << std::endl;
    
}


template<typename T>
inline void slice_vector(T* parent_vector, int first_pos, int second_pos, T* slice_vector)
{
    int count = 0;
    for (int i=first_pos; i<second_pos; i++)
    {
        slice_vector[count] = parent_vector[i];
        count++;
    }
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