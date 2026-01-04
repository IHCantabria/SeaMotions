
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

#pragma once

// Include general usage libraries
#include <functional>
#include <string>
#include <optional>

// Include local modules
#include "../config.hpp"

// Define local type to compactify expressions
typedef signed long long sll_type;

//////////////////////////////////////////////
////// MATHEMATICAL CONSTANTS BLOCK //////////
//////////////////////////////////////////////
const cusfloat PI = 3.141592653589793;
const cusfloat EULERGAMMA = 0.577215664901533;
const cusfloat LOG2_GAMMA = std::log( 2.0 ) - EULERGAMMA;


//////////////////////////////////////////////
////////// MACRO DEFINITION BLOCK ////////////
//////////////////////////////////////////////
#define POW2S(x) ((x)*(x))


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
                                            cusfloat    angfreq_to_freq( 
                                                                                    cusfloat angfreq 
                                                                        );

                                            cusfloat    angfreq_to_period( 
                                                                                    cusfloat angfreq 
                                                                        );
                                                                    
                                            int         assert_complex_equality(
                                                                                    cuscomplex  u, 
                                                                                    cuscomplex  v, 
                                                                                    cusfloat    epsilon
                                                                                );
        
template<typename T>                inline  int         assert_scalar_equality(         
                                                                                    T           &u, 
                                                                                    T           &v, 
                                                                                    T           epsilon
                                                                                );

template<typename T>                inline  int         assert_scalar_equality(         
                                                                                    T           &u, 
                                                                                    T           &v, 
                                                                                    T           abs_eps,
                                                                                    T           rel_eps
                                                                                );

template<typename T>                inline  int         assert_vector_equality(         
                                                                                    int         N, 
                                                                                    T*          u, 
                                                                                    T*          v, 
                                                                                    cusfloat    epsilon
                                                                                );

template<typename T>                inline  int         assert_vector_equality(         
                                                                                    int N, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    cusfloat abs_eps,
                                                                                    cusfloat rel_eps
                                                                                );

template<typename T>                inline  int         assert_vector_equality(
                                                                                    int N, 
                                                                                    int* u, 
                                                                                    int* v, 
                                                                                    int epsilon
                                                                                );
        
                                            void        bisection(
                                                                                    std::function<cusfloat(cusfloat)> f_def, 
                                                                                    cusfloat a, 
                                                                                    cusfloat b, 
                                                                                    cusfloat abs_prec, 
                                                                                    cusfloat rel_prec, 
                                                                                    int max_iter, 
                                                                                    bool verbose,
                                                                                    cusfloat &sol, 
                                                                                    int &info
                                                                    );

template <typename T>               inline  void        bisection_false(
                                                                                    T&          f_def, 
                                                                                    cusfloat    a, 
                                                                                    cusfloat    b, 
                                                                                    cusfloat    fabs_tol, 
                                                                                    cusfloat    xrel_tol, 
                                                                                    int         max_iter, 
                                                                                    bool        verbose,
                                                                                    cusfloat    &sol,
                                                                                    int         &info
                                                                        );

                                            cusfloat    check_zero_eps( 
                                                                                    cusfloat value,
                                                                                    cusfloat eps
                                                                    );
        
template<typename T>                inline  void        clear_vector(                    
                                                                                    int n, 
                                                                                    T* vec 
                                                                    );

template<typename T, std::size_t N> inline  void        clear_vector(                    
                                                                                    T* vec 
                                                                    );
        
template<typename T>                inline  void        copy_vector(                    
                                                                                    int n, 
                                                                                    T* reference_vector, 
                                                                                    T* target_vector
                                                                    );

template<typename T, std::size_t N> inline  void        copy_vector(
                                                                                    T* reference_vector, 
                                                                                    T* target_vector
                                                                   );

                                            void        conj_vector(
                                                                                    int             n,
                                                                                    cuscomplex*     u,
                                                                                    cuscomplex*     v
                                                                    );

                                            cusfloat    cos3_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                    );

                                            cusfloat    cos2sin_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                        );

                                            cusfloat    cossin2_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                        );

template<typename T>                inline  void        cross(
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                            );

                                            cusfloat    deg_to_rad(    
                                                                                    cusfloat deg 
                                                                    );

template<typename T>                inline  T           eucledian_dist(
                                                                                    int     np,
                                                                                    T*      v0,
                                                                                    T*      v1
                                                                        );
        
                                            sll_type    factorial(
                                                                                    int n
                                                                );
                                
                                            cusfloat    freq_to_angfreq( 
                                                                                    cusfloat freq 
                                                                        );

template<typename T>                        inline  T*  generate_empty_vector(          
                                                                                    int size
                                                                            );

template<typename T>                        inline  T   l2_norm( 
                                                                                    int N,
                                                                                    T* x,
                                                                                    T* y
                                                                );

                                            void        newton_raphson(
                                                                                    std::function<cusfloat(cusfloat)> f_def, 
                                                                                    std::function<cusfloat(cusfloat)> f_der_def,
                                                                                    cusfloat x0, 
                                                                                    cusfloat fabs_tol, 
                                                                                    cusfloat xrel_tol, 
                                                                                    int max_iter, 
                                                                                    bool verbose, 
                                                                                    cusfloat &sol, 
                                                                                    int &info
                                                                        );

                                            cusfloat    period_to_angfreq( 
                                                                                    cusfloat period 
                                                                        );
        
template<typename T>                inline  T           pow2s(
                                                                                    T x
                                                            );
        
template<typename T>                inline  T           pow3s(
                                                                                    T x
                                                                );
        
template<typename T>                inline  void        print_matrix(
                                                                                    int num_rows, 
                                                                                    int num_cols, 
                                                                                    T* mat, 
                                                                                    int precision, 
                                                                                    int align, 
                                                                                    int scient_not
                                                                    );
        
template<typename T>                inline  void        print_vector(
                                                                                    int n, 
                                                                                    T* v, 
                                                                                    int mode, 
                                                                                    int precision
                                                                    );

/**
 * @brief Calculate geometric properties of a quadrilateral.
 * 
 * Computes the area, centroid, first-order moments, and second-order moments
 * of a quadrilateral defined by four vertices. Optionally, moments can be computed
 * in a custom reference system.
 * 
 * @param[in]  x          Array of 4 x-coordinates of quadrilateral vertices
 * @param[in]  y          Array of 4 y-coordinates of quadrilateral vertices
 * @param[in]  z          Array of 4 z-coordinates of quadrilateral vertices
 * @param[out] area       Quadrilateral area
 * @param[out] centroid   Array[3] containing the centroid coordinates (x, y, z)
 * @param[out] moments_fo Array[3] containing first-order moments about the centroid
 * @param[out] moments_so Array[3] containing second-order moments (Ixx, Iyy, Ixy) about the centroid
 * @param[in]  ref_sys    Optional reference system origin for moment calculation.
 *                        If not provided, moments are computed about the centroid.
 */
                                            void        quad_geom_properties( 
                                                                                    const   cusfloat*   x,
                                                                                    const   cusfloat*   y,
                                                                                    const   cusfloat*   z,
                                                                                            cusfloat&   area,
                                                                                            cusfloat*   centroid,
                                                                                            cusfloat*   moments_fo,
                                                                                            cusfloat*   moments_so,
                                                                                    const   std::optional<cusfloat*>& ref_sys = std::nullopt
                                                                                );
        
template <typename T>               inline  int         sign(
                                                                                    T val
                                                                );

                                            cusfloat    sin_alpha( 
                                                                                    int         alpha,
                                                                                    cusfloat    theta
                                                                );

                                            cusfloat    sin3_int_0_2PI( 
                                                                                    int m,
                                                                                    int n,
                                                                                    int p
                                                                        );
        
template<typename T>                inline  void        slice_vector(
                                                                                    T* parent_vector, 
                                                                                    int first_pos, 
                                                                                    int second_pos, 
                                                                                    T* slice_vector
                                                                    );
        
template<typename T>                inline  void        sv_add(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );

template<typename T>        inline  void                sv_add(
                                                                                    int n, 
                                                                                    T alpha, 
                                                                                    T* u, 
                                                                                    T beta, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_add(                                
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T s, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        sv_cbrt(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_div(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        svs_div(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T s,
                                                                                    T* w
                                                                );

template<typename T>                inline  T           sv_dot(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T* v
                                                                );

template<typename T>                inline  void        sv_dot_vm(
                                                                                    int n,
                                                                                    T*  vec,
                                                                                    T*  mat,
                                                                                    T*  vec_out
                                                                );
        
template<typename T>                inline  void        sv_inv(
                                                                                    int n, 
                                                                                    T s, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_mod(                     int n, 
                                                                                    T* u, 
                                                                                    T &mod
                                                                );
        
template<typename T>                inline  void        sv_mult(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_mult(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T s, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        sv_pow(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T* v,
                                                                                    T* w
                                                                );

template<typename T>                inline  void        svs_pow(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T s,
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_pow2(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );

template<typename T>                inline  void        sv_pow3(
                                                                                    int n,
                                                                                    T* u,
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_sqrt(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        sv_sub(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T* v, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_sub(
                                                                                    int n, 
                                                                                    T* u, 
                                                                                    T s, 
                                                                                    T* w
                                                                );
        
template<typename T>                inline  void        svs_sub(
                                                                                    int n, 
                                                                                    T s, 
                                                                                    T* u, 
                                                                                    T* w
                                                                );



/**
 * @brief Calculate geometric properties of a triangle.
 * 
 * Computes the area, centroid, first-order moments, and second-order moments
 * of a triangle defined by three vertices. Optionally, moments can be computed
 * in a custom reference system.
 * 
 * @param[in]  x          Array of 3 x-coordinates of triangle vertices
 * @param[in]  y          Array of 3 y-coordinates of triangle vertices
 * @param[in]  z          Array of 3 z-coordinates of triangle vertices
 * @param[out] area       Triangle area
 * @param[out] centroid   Array[3] containing the centroid coordinates (x, y, z)
 * @param[out] moments_fo Array[3] containing first-order moments about the centroid
 * @param[out] moments_so Array[3] containing second-order moments (Ixx, Iyy, Ixy) about the centroid
 * @param[in]  ref_sys    Optional reference system origin for moment calculation.
 *                        If not provided, moments are computed about the centroid.
 */
                                            void        triangle_geom_properties( 
                                                                                    const   cusfloat*   x,
                                                                                    const   cusfloat*   y,
                                                                                    const   cusfloat*   z,
                                                                                            cusfloat&   area,
                                                                                            cusfloat*   centroid,
                                                                                            cusfloat*   moments_fo,
                                                                                            cusfloat*   moments_so,
                                                                                    const   std::optional<cusfloat*>& ref_sys = std::nullopt
                                                                                );


#include "math_tools.txx"
