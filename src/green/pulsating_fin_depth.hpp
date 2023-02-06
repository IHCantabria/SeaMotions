
#ifndef __pulsating_hpp
#define __pulsating_hpp

// Include general usage libraries
#include <tuple>

// Include local modules
#include "../config.hpp"
#include "../waves.hpp"

// Include namespaces objects into the local namespace
using namespace std;


cuscomplex john_series(
                        cusfloat R,
                        cusfloat z, 
                        cusfloat zeta, 
                        cusfloat h,
                        WaveDispersionData &wave_data
                        );

cuscomplex john_series(
                        cusfloat x,
                        cusfloat y,
                        cusfloat z,
                        cusfloat xi,
                        cusfloat eta,
                        cusfloat zeta,
                        cusfloat h,
                        WaveDispersionData &wave_data
                        );

tuple<cuscomplex, cuscomplex> john_series_dhoriz(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionData &wave_data
                                                );

cuscomplex john_series_dr(
                            cusfloat R,
                            cusfloat z,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            );

cuscomplex john_series_dx(
                            cusfloat x,
                            cusfloat y,
                            cusfloat z,
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            );

cuscomplex john_series_dy(
                            cusfloat x,
                            cusfloat y,
                            cusfloat z,
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            );

cuscomplex john_series_dz(
                            cusfloat R,
                            cusfloat z,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            );

                        
cuscomplex john_series_dz(
                            cusfloat x,
                            cusfloat y,
                            cusfloat z,
                            cusfloat xi,
                            cusfloat eta,
                            cusfloat zeta,
                            cusfloat h,
                            WaveDispersionData &wave_data
                            );

#endif