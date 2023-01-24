
#ifndef __waves_hpp
#define __waves_hpp

// Include local modules
#include "config.hpp"
#include "./math/math_tools.hpp"


// Declare module functions
cusfloat w2k(cusfloat w, cusfloat h, cusfloat g);
void w2ki(cusfloat w, cusfloat h, cusfloat g, int n, cusfloat* kn);


// Define modules classes
/**
 * @brief WaveDisperionData container
 * 
 * This class contains information relative to the wave dispersion at some predefined
 * depth and an angular frequency w0. This class also contains information derived from 
 * the dispersion relation as the John series coefficients.
 * 
 */
struct WaveDispersionData
{
    cusfloat grav_acc = 0.0; /** Gravitational acceleration*/
    bool is_john = false; /** Flag to track if John's coefficients are calculated*/
    cusfloat k0 = 0.0; /** Real wave number*/
    cusfloat k0nu = 0.0; /** Real wave number John's coefficient*/
    cusfloat* kn = nullptr; /** Imaginary wave numbers*/
    cusfloat* knnu = nullptr; /** Imaginary wave number John's coefficient*/
    cusfloat nu = 0.0; /** Real wave number at infinite depth: w^2/g*/
    int num_kn = 0; /** Number of imaginary wave number requested*/
    cusfloat w0 = 0.0; /** Wave angular frequency*/
    cusfloat water_depth = 0.0; /** Working water depth*/

    WaveDispersionData(cusfloat w0_i, int num_kn_i, cusfloat water_depth_i, cusfloat grav_acc_i)
    /**
     * @brief WaveDispersionData constructor
     * 
     * \param w0_i Wave Angular frequency
     * \param num_kn_i Number of imaginary wave numbers requested
     * \param water_depth Working water depth
     * \param grav_acc_i Gravitational acceleration
     */
    {
        // Save input arguments
        this->grav_acc = grav_acc_i;
        this->num_kn = num_kn_i;
        this->w0 = w0_i;
        this->water_depth = water_depth_i;

        // Calculate derivative properties
        this->nu = pow2s(this->w0)/this->grav_acc;

        // Calculate dispersion relation terms
        this->_calculate_dispersion();
        
    }

    ~WaveDispersionData()
    {
        if (this->num_kn > 0)
        {
            delete [] this->kn;
            delete [] this->knnu;
        }
    }

    void _calculate_dispersion(void)
    /**
     * @brief Calcualte real and imaginary wave numbers. Hiden method.
     * 
     * \return
     */
    {
        // Calculate real wave number
        this->k0 = w2k(this->w0, this->water_depth, this->grav_acc);

        // Calculate imag wave number
        if (this->num_kn > 0)
        {
            this->kn = new cusfloat [this->num_kn];
            w2ki(this->w0, this->water_depth, this->grav_acc, this->num_kn, this->kn);
        }

    }

    void calculate_john_terms(void)
    /**
     * @brief Calculate John series parameters
     * 
     * The John series to calculate the finite water depth green function has
     * some paramaters that are only dependent on the wave dispersion relation.
     * Therefore, in order to save time it is better to have them precomputed.
     * 
     * \return
     */
    {
        // Rename object properties to fit better on the formulas
        cusfloat k0 = this->k0;
        cusfloat nu = this->nu;
        cusfloat h = this->water_depth;

        // Calculate real wave number coefficient
        this->k0nu = (
                    2*PI*pow2s(k0)
                    /
                    (((pow2s(k0)-pow2s(nu))*h+nu)*pow2s(1+exp(-2*k0*h)))
                    );
        
        // Calculate imaginary wave number coefficient
        if (this->num_kn > 0)
        {
            cusfloat kni = 0.0;
            this->knnu = new cusfloat [this->num_kn];
            cusfloat knnu_i = 0.0;
            for (int i=0; i<this->num_kn; i++)
            {
                kni = kn[i];
                knnu_i = pow2s(kni)+pow2s(nu);
                this->knnu[i] = knnu_i/(knnu_i*h-nu);
            }
        }

        // Activate flag to calculated mode
        this->is_john = true;

    }

    void print(void)
    {
        std::cout << "Wave dispersion properties for:" << std::endl;
        std::cout << "------> w0: " << this->w0 << std::endl;
        std::cout << "------> water_depth: " << this->water_depth << std::endl;
        std::cout << "------> grav_acc: " << this->grav_acc << std::endl;
        std::cout << " - Real wave Nmber (k0): " << this->k0 << std::endl;
        std::cout << " - Imaginary wave Nmber (kn): " << this->k0 << std::endl;
        for (int i=0; i<this->num_kn; i++)
        {
            std::cout << " ---> kn[" << i << "]: " << this->kn[i] << std::endl;
        }
    }
};

#endif