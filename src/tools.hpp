
#ifndef __tools_hpp
#define __tools_hpp

#include "config.hpp"


template<typename T> int assert_vector_equality(T (&u)[3], T (&v)[3], T epsilon);
template<> int assert_vector_equality<int>(int (&u)[3], int (&v)[3], int epsilon);
cusfloat* generate_empty_vector(int size);
template<typename T> void cross(T (&u)[3], T (&v)[3], T (&w)[3]);

#include "tools.txx"

#endif