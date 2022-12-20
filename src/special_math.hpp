
#ifndef __special_hpp
#define __special_hpp

// Include local modules
#include "config.hpp"


cusfloat besseli0(cusfloat x);
cusfloat besseli1(cusfloat x);
cusfloat besselj0(cusfloat x);
cusfloat besselj1(cusfloat x);
cusfloat besselk0(cusfloat x);
cusfloat besselk1(cusfloat x);
cusfloat bessely0(cusfloat x);
cusfloat bessely1(cusfloat x);
cusfloat chebyshev_poly_raw(int n, cusfloat x);
cusfloat expint_i(cusfloat x);
cusfloat legendre_poly_raw(int n, cusfloat x);
cusfloat polynomial_f0(cusfloat x);
cusfloat polynomial_f1(cusfloat x);
cusfloat polynomial_th0(cusfloat x);
cusfloat polynomial_th1(cusfloat x);
cusfloat psi_fun(int n);
cusfloat rational_fraction_f0(cusfloat x);
cusfloat rational_fraction_f1(cusfloat x);
cusfloat rational_fraction_th0(cusfloat x);
cusfloat rational_fraction_th1(cusfloat x);
cusfloat struve0(cusfloat x);
cusfloat struve1(cusfloat x);

#endif