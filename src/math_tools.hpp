
#ifndef __math_tools_hpp
#define __math_tools_hpp

#include "config.hpp"


template<typename T> inline int assert_vector_equality(int N, T* u, T* v, T epsilon);
template<> inline int assert_vector_equality<int>(int N, int* u, int* v, int epsilon);
cusfloat* generate_empty_vector(int size);
template<typename T> void cross(T (&u)[3], T (&v)[3], T (&w)[3]);
template<typename T> void vAdd(int N, T* u, T* v, T* w);
template<typename T> void vsAdd(int N, T* u, T s, T* w);
template<typename T> void vSub(int N, T* u, T* v, T* w);
template<typename T> void vsSub(int N, T* u, T s, T* w);
template<typename T> void vsSub(int N, T s, T* u, T* w);

#include "math_tools.txx"

#endif