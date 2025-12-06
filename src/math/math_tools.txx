
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

// Include general usage libraries
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <sstream>

// Include general usage scientific libraries
#include <cmath>
#include "mkl.h"

// Include local modules
#include "../tools.hpp"


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
inline int assert_scalar_equality(T &u, T &v, T abs_eps, T rel_eps)
{
    int pass = 1;
    if ( std::abs( u - v ) > ( std::abs( u * rel_eps ) + abs_eps ) )
    {
        pass = 0;
    }

    return pass;
}


template<typename T>
inline int assert_vector_equality(int N, T* u, T* v, cusfloat epsilon)
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

template<typename T>
inline int assert_vector_equality(int N, T* u, T* v, cusfloat abs_eps, cusfloat rel_eps)
{
    int pass = 1;
    for (int i=0; i<N; i++)
    {
        if ( std::abs( u[i] - v[i] ) > ( std::abs( rel_eps*u[i] ) + abs_eps ) )
        {
            pass = 0;
        }
    }

    return pass;
}


template<typename T>
inline int assert_vector_equality(int N, int* u, int* v, int epsilon)
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
void clear_vector( int n, T* vec )
{
    for ( int i=0; i<n; i++ )
    {
        vec[i] = 0.0;
    }
}


template<typename T, std::size_t N>
void clear_vector( T* vec )
{
    for ( std::size_t i=0; i<N; i++ )
    {
        vec[i] = 0.0;
    }
}

template<typename T>
void copy_vector(int n, T* reference_vector, T* target_vector)
{
    for(int i=0; i<n; i++)
    {
        target_vector[i] = reference_vector[i];
    }
}


template<typename T, std::size_t N> 
void  copy_vector( T* reference_vector, T* target_vector )
{
    for ( std::size_t i=0; i<N; i++ )
    {
        target_vector[i] = reference_vector[i];
    }
}


template<typename T>
inline void cross( T* u, T* v, T* w )
{
    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}


template<typename T>
inline T eucledian_dist( int np, T* v0, T* v1 )
{
    T dist = 0.0;
    for ( int i=0; i<np; i++ )
    {
        dist += pow2s( v1[i] - v0[i] );
    }
    dist = std::sqrt( dist );

    return dist;
}


template<typename T>
T* generate_empty_vector(int size)
{
    return (T*)mkl_calloc( size, sizeof(T), FLOATING_PRECISION );
}


template<typename T>
T l2_norm( int n, T* x, T* y )
{
    T diff  = 0.0;
    T sum   = 0.0;
    for ( int i=0; i<n; i++ )
    {
        diff    = x[i] - y[i];
        sum     += diff * diff;
    }
    return std::sqrt( sum );
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


template <typename T> int sign(T val) {
    return (T(0) < val) - (val < T(0));
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
inline void sv_add(int n, T alpha, T* u, T beta, T* v, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = alpha*u[i] + beta*v[i];
    }
}


template<typename T>
inline void sv_cbrt(int n, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = std::cbrt(u[i]);
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
inline void svs_div(int n, T* u, T s, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i]/s;
    }
}


template<typename T>    
inline  T   sv_dot(
                        int n,
                        T* u,
                        T* v
                    )
{
    T value = 0;
    for ( int i=0; i<n; i++ )
    {
        value += u[i] * v[i];
    }

    return value;
}


template<typename T>
inline void sv_dot_vm(
                        int n,
                        T*  vec,
                        T*  mat,
                        T*  vec_out
                    )
{
    for ( int i=0; i<n; i++ )
    {
        vec_out[i]  = 0.0;
        for ( int j=0; j<n; j++ )
        {
            vec_out[i] += vec[j] * mat[j*n+i];
        }
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
inline void sv_pow(int n, T* u, T* v, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = std::pow( u[i], v[i] );
    }
}


template<typename T>
inline void svs_pow(int n, T* u, T s, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = std::pow( u[i], s );
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
inline void sv_pow3(int n, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = u[i]*u[i]*u[i];
    }
}


template<typename T>
inline void sv_sqrt(int n, T* u, T* w)
{
    for (int i=0; i<n; i++)
    {
        w[i] = std::sqrt( u[i] );
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