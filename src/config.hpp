
#ifndef __config_hpp
#define __config_hpp

//#define PRECISION "SIMPLE"
#ifdef SIMPLE_PREC
int FLOATING_PRECISION = 32;
typedef float cusfloat;
#else
int FLOATING_PRECISION = 64;
typedef double cusfloat;
#endif

#endif