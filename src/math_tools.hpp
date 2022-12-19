
#ifndef __math_tools_hpp
#define __math_tools_hpp


// Include general usage libraries
#include <string>

// Include local modules
#include "config.hpp"

//////////////////////////////////////////////
////// MATHEMATICAL CONSTANTS BLOCK //////////
//////////////////////////////////////////////
const cusfloat PI = 3.141592653589793;
const cusfloat EULERGAMMA = 0.577215664901533;


//////////////////////////////////////////////
/////// FUNCTION DEFINITION BLOCK ////////////
//////////////////////////////////////////////
template<typename T> inline int assert_scalar_equality(T &u, T &v, T epsilon);
template<typename T> inline int assert_vector_equality(int N, T* u, T* v, T epsilon);
template<> inline int assert_vector_equality<int>(int N, int* u, int* v, int epsilon);
cusfloat* generate_empty_vector(int size);
template<typename T> void copy_vector(int n, T* reference_vector, T* target_vector);
template<typename T> inline void cross(T (&u)[3], T (&v)[3], T (&w)[3]);
signed long long factorial(int n);
template<typename T> inline T pow2s(T x);
template<typename T> inline T pow3s(T x);
template<typename T> inline void print_matrix(int num_rows, int num_cols, T* mat, int precision, int align, int scient_not);
template<typename T> inline void print_vector(int n, T* v, int mode, int precision);
template<typename Functor> cusfloat romberg_quadrature(Functor f, cusfloat a, cusfloat b, double precision);
template<typename T> inline void slice_vector(T* parent_vector, int first_pos, int second_pos, T* slice_vector);
template<typename T> inline void sv_add(int n, T* u, T* v, T* w);
template<typename T> inline void svs_add(int n, T* u, T s, T* w);
template<typename T> inline void sv_div(int n, T* u, T* v, T* w);
template<typename T> inline void sv_inv(int n, T s, T* u, T* w);
template<typename T> inline void sv_mod(int n, T* u, T &mod);
template<typename T> inline void sv_mult(int n, T* u, T* v, T* w);
template<typename T> inline void svs_mult(int n, T* u, T s, T* w);
template<typename T> inline void sv_pow2(int n, T* u, T* w);
template<typename T> inline void sv_sqrt(int n, T* u, T* w);
template<typename T> inline void sv_sub(int n, T* u, T* v, T* w);
template<typename T> inline void svs_sub(int n, T* u, T s, T* w);
template<typename T> inline void svs_sub(int n, T s, T* u, T* w);


#include "math_tools.txx"

#endif