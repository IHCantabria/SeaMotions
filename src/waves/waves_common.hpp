
#ifndef __waves_common_hpp
#define __waves_common_hpp

// Include local modules
#include "../config.hpp"
#include "../math/math_tools.hpp"
#include "../waves/wave_dispersion_so.hpp"


// Declare module functions
cuscomplex  wave_potential_fo_space( 
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

cuscomplex  wave_potential_fo_space_dx( 
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

cuscomplex  wave_potential_fo_space_dy( 
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

cuscomplex  wave_potential_fo_space_dz( 
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

cuscomplex  wave_potential_so_space( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cuscomplex  wave_potential_so_space_dx( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cuscomplex  wave_potential_so_space_dy( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cuscomplex  wave_potential_so_space_dz( 
                                            cusfloat            x,
                                            cusfloat            y,
                                            cusfloat            z,
                                            WaveDispersionSO*   wd,
                                            int                 qtf_type
                                    );

cusfloat    wave_vertical_profile_fo(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                    );

cusfloat    wave_vertical_profile_fo_dz(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                        );

cusfloat    wave_vertical_profile_mod_fo(
                                            cusfloat    k,
                                            cusfloat    h,
                                            cusfloat    z
                                        );

#endif