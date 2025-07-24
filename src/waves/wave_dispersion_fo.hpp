
#pragma once

// Include local modules
#include "../config.hpp"


// Define modules classes
/**
 * @brief WaveDisperionData container
 * 
 * This class contains information relative to the wave dispersion at some predefined
 * depth and an angular frequency w0. This class also contains information derived from 
 * the dispersion relation as the John series coefficients.
 * 
 */
template<int NumKn>
struct WaveDispersionFO
{
private:
    // Define private class methods
    void    _calculate_dispersion(
                                    void
                                );

    void    _calculate_john_terms(
                                    void
                                );
public:
    // Define class attributes
    cusfloat        grav_acc    = 0.0;      /** Gravitational acceleration*/
    bool            is_john     = false;    /** Flag to track if John's coefficients are calculated*/
    cusfloat        k0          = 0.0;      /** Real wave number*/
    cusfloat        k0nu        = 0.0;      /** Real wave number John's coefficient*/
    cusfloat        kn[NumKn];              /** Imaginary wave numbers*/
    cusfloat        knnu[NumKn];            /** Imaginary wave number John's coefficient*/
    cusfloat        nu          = 0.0;      /** Real wave number at infinite depth: w^2/g*/
    const int       num_kn      = NumKn;    /** Number of imaginary wave number requested*/
    cusfloat        w0          = 0.0;      /** Wave angular frequency*/
    cusfloat        water_depth = 0.0;      /** Working water depth*/

    // Check that the number if imaginary input frequencies is higher than 0
    static_assert( NumKn > 0 );

    // Define class constructors and destructor
    WaveDispersionFO(
                                    cusfloat    w0_i, 
                                    cusfloat    water_depth_i, 
                                    cusfloat    grav_acc_i
                    );

    // Define class public methods
    void    print(
                                    void
                );

    
    void    update(  
                                    cusfloat    w0_i, 
                                    cusfloat    water_depth_i, 
                                    cusfloat    grav_acc_i
                    );


    void    update_full(  
                                    cusfloat    w0_i, 
                                    cusfloat    water_depth_i, 
                                    cusfloat    grav_acc_i
                        );
    
};

// Include template definitions
#include "wave_dispersion_fo.txx"

// Create usefull alias to avoid using template expression
using WaveDispersionFONK = WaveDispersionFO<NUM_KN>;