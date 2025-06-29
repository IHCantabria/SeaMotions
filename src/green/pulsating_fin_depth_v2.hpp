
#pragma once

// Include general usage libraries
#include <tuple>

// Include local modules
#include "../config.hpp"
#include "./integrals_db.hpp"
#include "../waves/wave_dispersion_fo.hpp"

// Include namespaces objects into the local namespace
using namespace std;

// Define type alias
using tuplecc = tuple<cuscomplex,cuscomplex>;

// Define module constants
#define     GREEN_ZEROTH_DR     1e-4
#define     GREEN_DR_EPS        1e-6

// Declare module functions
cuscomplex  G_integral(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                    );


cuscomplex  G_integral(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                    );


cuscomplex  G_integral_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                        );


cuscomplex  G_integral_steady_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h
                                );


cuscomplex  G_integral_wave_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                                );
                    

cuscomplex  G_integral_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                        );


cuscomplex  G_integral_steady_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h
                                );


cuscomplex  G_integral_wave_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                                );


cuscomplex  G_integral_steady(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h
                            );


cuscomplex  G_integral_wave(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                            );


cuscomplex  john_series(
                                                cusfloat R,
                                                cusfloat z, 
                                                cusfloat zeta, 
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                        );


cuscomplex  john_series(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                        );


tuplecc     john_series_dhoriz(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                                );


cuscomplex  john_series_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                            );


cuscomplex  john_series_dx(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                            );


cuscomplex  john_series_dy(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                            );


cuscomplex  john_series_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                            );

                        
cuscomplex  john_series_dz(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data
                            );


cuscomplex  wave_term_fin_depth(
                                                cusfloat x,
                                                cusfloat y,
                                                cusfloat z,
                                                cusfloat xi,
                                                cusfloat eta,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                                );

cuscomplex  wave_term_fin_depth_integral(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                                        );


cuscomplex  wave_term_fin_depth_integral_dr(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                                            );


cuscomplex  wave_term_fin_depth_integral_dz(
                                                cusfloat R,
                                                cusfloat z,
                                                cusfloat zeta,
                                                cusfloat h,
                                                WaveDispersionFO &wave_data,
                                                IntegralsDb &idb
                                            );
