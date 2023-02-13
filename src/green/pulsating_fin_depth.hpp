
#ifndef __pulsating_fin_depth
#define __pulsating_fin_depth

// Include general usage libraries
#include <tuple>

// Include local modules
#include "../config.hpp"
#include "./integrals_db.hpp"
#include "../waves.hpp"

// Include namespaces objects into the local namespace
using namespace std;


cuscomplex G_integral(
                        cusfloat R,
                        cusfloat z,
                        cusfloat zeta,
                        cusfloat h,
                        WaveDispersionData &wave_data,
                        IntegralsDb &idb
                    );


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


cuscomplex wave_term_fin_depth(
                                cusfloat x,
                                cusfloat y,
                                cusfloat z,
                                cusfloat xi,
                                cusfloat eta,
                                cusfloat zeta,
                                cusfloat h,
                                WaveDispersionData &wave_data,
                                IntegralsDb &idb
                                );

cuscomplex wave_term_fin_depth_integral(
                                        cusfloat R,
                                        cusfloat z,
                                        cusfloat zeta,
                                        cusfloat h,
                                        WaveDispersionData &wave_data,
                                        IntegralsDb &idb
                                        );

#endif