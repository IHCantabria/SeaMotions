
#ifndef __config_hpp
#define __config_hpp

#ifdef SIMPLE_PREC
typedef float cusfloat;
constexpr int FLOATING_PRECISION = 32;
constexpr cusfloat EPS_PRECISION = 1e-6;
constexpr cusfloat EPS_PRECISION_ORDER = -6;
#else
typedef double cusfloat;
constexpr int FLOATING_PRECISION = 32;
constexpr cusfloat EPS_PRECISION = 1e-14;
constexpr cusfloat EPS_PRECISION_ORDER = -14;
#endif

#endif