
#ifndef __math_tools_hpp
#define __math_tools_hpp

#include "config.hpp"


template<typename T> inline int assert_scalar_equality(T &u, T &v, T epsilon);
template<typename T> inline int assert_vector_equality(int N, T* u, T* v, T epsilon);
template<> inline int assert_vector_equality<int>(int N, int* u, int* v, int epsilon);
cusfloat* generate_empty_vector(int size);
template<typename T> inline void cross(T (&u)[3], T (&v)[3], T (&w)[3]);
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
template<typename T> inline void print_vector(int n, T* v, int mode);


#include "math_tools.txx"

#endif