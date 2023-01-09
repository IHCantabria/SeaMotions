
#ifndef __chebyshev_hpp
#define __chebyshev_hpp

#include "config.hpp"

cusfloat chebyshev_poly_raw(int n, cusfloat x);
cusfloat chebyshev_poly_der_raw(int n, cusfloat x);
void chebyshev_poly_roots(int num_points, cusfloat* roots);
void get_gauss_chebyshev(int num_points, cusfloat* weights, cusfloat* roots);

#endif