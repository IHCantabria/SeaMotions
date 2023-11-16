
#ifndef __euler_transfomrs_hpp
#define __euler_transfomrs_hpp

// Include local modules
#include "../config.hpp"


void    _rot_x(
                    cusfloat*   mat,
                    cusfloat    alpha
                );


void    _rot_y(
                    cusfloat*   mat,
                    cusfloat    beta
                );


void    _rot_z(
                    cusfloat*   mat,
                    cusfloat    gamma
                );

#endif