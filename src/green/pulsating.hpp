
#ifndef __pulsating_hpp
#define __pulsating_hpp

#include "../../src/config.hpp"

cusfloat expint_inf_depth_num(cusfloat X, cusfloat Y, cusfloat t);
cusfloat expint_inf_depth_num_dx(cusfloat X, cusfloat Y, cusfloat t);
cusfloat wave_term_inf_depth(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num(cusfloat X, cusfloat Y);
cusfloat wave_term_inf_depth_num_dx(cusfloat X, cusfloat Y);

#endif