
#ifndef __tools_hpp
#define __tools_hpp

#include "config.hpp"


template<typename T> inline int assert_vector_equality(int N, T* u, T* v, T epsilon);
template<> inline int assert_vector_equality<int>(int N, int* u, int* v, int epsilon);
cusfloat* generate_empty_vector(int size);
template<typename T> void cross(T (&u)[3], T (&v)[3], T (&w)[3]);

#include "tools.txx"

#endif