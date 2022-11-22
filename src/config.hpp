
#ifndef __config_hpp
#define __config_hpp

#ifdef SIMPLE_PREC
extern int FLOATING_PRECISION;
typedef float cusfloat;
extern cusfloat EPS_PRECISION;
#else
extern int FLOATING_PRECISION;
typedef double cusfloat;
extern cusfloat EPS_PRECISION;
#endif

#endif