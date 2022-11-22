
// Include local modules
#include "config.hpp"

#ifdef SIMPLE_PREC
int FLOATING_PRECISION = 32;
cusfloat EPS_PRECISION = 1e-6;
#else
int FLOATING_PRECISION = 64;
cusfloat EPS_PRECISION = 1e-14;
#endif