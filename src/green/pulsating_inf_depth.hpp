
#ifndef __pulsating_inf_depth_hpp
#define __pulsating_inf_depth_hpp

// Include general usage libraries
#include <tuple>

// Include local modules
#include "integrals_db.hpp"

cusfloat calculate_dr_dx(cusfloat R, cusfloat dx);
cusfloat calculate_r(
                    cusfloat x,
                    cusfloat y,
                    cusfloat xi,
                    cusfloat eta
                    );

cusfloat wave_term_inf_depth(
                            cusfloat x_ndim,
                            cusfloat y_ndim,
                            IntegralsDb &idb
                            );
std::tuple<cusfloat, cusfloat> wave_term_inf_depth_dhoriz(
                                                            cusfloat x,
                                                            cusfloat y,
                                                            cusfloat z,
                                                            cusfloat xi,
                                                            cusfloat eta,
                                                            cusfloat zeta,
                                                            cusfloat nu,
                                                            IntegralsDb &idb
                                                            );
cusfloat wave_term_inf_depth_dx(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dxndim(
                                cusfloat x_ndim,
                                cusfloat y_ndim,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dy(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dyndim(
                                    cusfloat x_ndim,
                                    cusfloat y_ndim,
                                    IntegralsDb &idb
                                    );
cusfloat wave_term_inf_depth_dz(
                                cusfloat x_ndim,
                                cusfloat z,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );
cusfloat wave_term_inf_depth_dz(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat nu,
                                IntegralsDb &idb
                                );


#endif