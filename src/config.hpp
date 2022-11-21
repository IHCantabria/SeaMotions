
#ifndef __config_hpp
#define __config_hpp

#ifdef SIMPLE_PREC
extern int FLOATING_PRECISION;
typedef float cusfloat;
#else
extern int FLOATING_PRECISION;
typedef double cusfloat;
#endif

#endif