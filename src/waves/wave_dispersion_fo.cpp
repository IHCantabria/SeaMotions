
// Include general usage libraries
#include <iostream>

// Include local modules
#include "waves_common.hpp"
#include "wave_dispersion_fo.hpp"


void    WaveDispersionFO::_calculate_dispersion(void)
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


void    WaveDispersionFO::calculate_john_terms(void)
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


void    WaveDispersionFO::print(void)
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


WaveDispersionFO::WaveDispersionFO(
                                        cusfloat    w0_i, 
                                        int         num_kn_i, 
                                        cusfloat    water_depth_i, 
                                        cusfloat    grav_acc_i
                                    )
/**
 * @brief WaveDispersionFO constructor
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


WaveDispersionFO::~WaveDispersionFO( 
                                        void
                                    )
{
    if (this->num_kn > 0)
    {
        delete [] this->kn;
        delete [] this->knnu;
    }
}