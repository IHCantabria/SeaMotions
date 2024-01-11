
#ifndef __wave_dispersion_fo_hpp
#define __wave_dispersion_fo_hpp

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
struct WaveDispersionFO
{
public:
    // Define class attributes
    cusfloat    grav_acc    = 0.0;      /** Gravitational acceleration*/
    bool        is_john     = false;    /** Flag to track if John's coefficients are calculated*/
    cusfloat    k0          = 0.0;      /** Real wave number*/
    cusfloat    k0nu        = 0.0;      /** Real wave number John's coefficient*/
    cusfloat*   kn          = nullptr;  /** Imaginary wave numbers*/
    cusfloat*   knnu        = nullptr;  /** Imaginary wave number John's coefficient*/
    cusfloat    nu          = 0.0;      /** Real wave number at infinite depth: w^2/g*/
    int         num_kn      = 0;        /** Number of imaginary wave number requested*/
    cusfloat    w0          = 0.0;      /** Wave angular frequency*/
    cusfloat    water_depth = 0.0;      /** Working water depth*/

    // Define class constructors and destructor
    WaveDispersionFO(
                                    cusfloat    w0_i, 
                                    int         num_kn_i, 
                                    cusfloat    water_depth_i, 
                                    cusfloat    grav_acc_i
                    );

    ~WaveDispersionFO(
                                    void
                    );
                    
    // Define class functions
    void    _calculate_dispersion(
                                    void
                                );

    void    calculate_john_terms(
                                    void
                                );

    void    print(
                                    void
                );
};

#endif