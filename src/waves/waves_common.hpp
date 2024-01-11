
#ifndef __waves_common_hpp
#define __waves_common_hpp

// Include local modules
#include "../config.hpp"
#include "../math/math_tools.hpp"


// Declare module functions
cusfloat    k2w( 
                                            cusfloat k, 
                                            cusfloat h, 
                                            cusfloat g 
                );

cusfloat    w2k(
                                            cusfloat w, 
                                            cusfloat h, 
                                            cusfloat g
                );

void        w2ki(
                                            cusfloat w, 
                                            cusfloat h, 
                                            cusfloat g, 
                                            int n, 
                                            cusfloat* kn
                );

cuscomplex  wave_potential_airy_space( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                    );

cuscomplex  wave_potential_airy_space_dx( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        );

cuscomplex  wave_potential_airy_space_dy( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        );

cuscomplex  wave_potential_airy_space_dz( 
                                            cusfloat aw,
                                            cusfloat w,
                                            cusfloat k,
                                            cusfloat h,
                                            cusfloat g,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z,
                                            cusfloat mu
                                        );

#endif