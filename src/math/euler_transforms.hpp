
#ifndef __euler_transfomrs_hpp
#define __euler_transfomrs_hpp

// Include local modules
#include "../config.hpp"


void    euler_local_to_global(
                                        cusfloat    alpha,
                                        cusfloat    beta,
                                        cusfloat    gamma,
                                        cusfloat*   rot_mat
                            );


void    euler_local_to_global_disp(
                                        cusfloat*   dofs_trans,
                                        cusfloat*   dofs_rot,
                                        cusfloat*   radius,
                                        cusfloat*   displacement
                                    );


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