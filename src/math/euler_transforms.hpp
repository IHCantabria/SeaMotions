
#ifndef __euler_transfomrs_hpp
#define __euler_transfomrs_hpp

// Include local modules
#include "../config.hpp"



template<typename T>    void    euler_local_to_global(
                                                                T           alpha,
                                                                T           beta,
                                                                T           gamma,
                                                                T*          rot_mat
                                                    );


template<typename T>    void    euler_local_to_global_disp(
                                                                T*          dofs_trans,
                                                                T*          dofs_rot,
                                                                cusfloat*   radius,
                                                                T*          displacement
                                                            );


template<typename T>    void    _rot_x(
                                                                T*          mat,
                                                                T           alpha
                                        );


template<typename T>    void    _rot_y(
                                                                T*          mat,
                                                                T           beta
                                        );


template<typename T>    void    _rot_z(
                                                                T*          mat,
                                                                T           gamma
                                        );


#include "euler_transforms.txx"

#endif