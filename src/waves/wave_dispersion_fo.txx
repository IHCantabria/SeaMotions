
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// Include general usage libraries
#include <iostream>

// Include local modules
#include "wave_dispersion_fo.hpp"

#include "../math/math_tools.hpp"
#include "wave_dispersion_base_fo.hpp"


template<int NumKn>
void    WaveDispersionFO<NumKn>::_calculate_dispersion( void )
/**
 * @brief Calcualte real and imaginary wave numbers. Hiden method.
 * 
 * \return
 */
{
    // Calculate real wave number
    this->k0 = w2k(this->w0, this->water_depth, this->grav_acc);

    // Calculate imag wave number
    w2ki(this->w0, this->water_depth, this->grav_acc, this->num_kn, this->kn);

}


template<int NumKn>
void    WaveDispersionFO<NumKn>::_calculate_john_terms( void )
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
                    2 * PI * pow2s( k0 )
                    /
                    (
                        ( ( pow2s( k0 ) - pow2s( nu ) ) * h + nu )
                        *
                        pow2s( 1.0 + exp(-2.0 * k0 * h ) )
                    )
                );
    
    // Calculate imaginary wave number coefficient
    cusfloat knnu_i = 0.0;
    for (int i=0; i<NumKn; i++)
    {
        knnu_i          = pow2s( this->kn[i] ) + pow2s( nu );
        this->knnu[i]   = knnu_i/(knnu_i*h-nu);
    }

    // Activate flag to calculated mode
    this->is_john = true;

}


template<int NumKn>
void    WaveDispersionFO<NumKn>::print(void)
{
    std::cout << "Wave dispersion properties for:" << std::endl;
    std::cout << "------> w0: " << this->w0 << std::endl;
    std::cout << "------> water_depth: " << this->water_depth << std::endl;
    std::cout << "------> grav_acc: " << this->grav_acc << std::endl;
    std::cout << " - Real wave Nmber (k0): " << this->k0 << std::endl;
    std::cout << " - Imaginary wave Nmber (kn): " << this->k0 << std::endl;
    for (int i=0; i<NumKn; i++)
    {
        std::cout << " ---> kn[" << i << "]: " << this->kn[i] << std::endl;
    }
}


template<int NumKn>
void    WaveDispersionFO<NumKn>::update(
                                            cusfloat    w0_i, 
                                            cusfloat    water_depth_i, 
                                            cusfloat    grav_acc_i
                                        )
{
    // Save input arguments
    this->grav_acc      = grav_acc_i;
    this->w0            = w0_i;
    this->water_depth   = water_depth_i;

    // Calculate derivative properties
    this->nu = pow2s(this->w0)/this->grav_acc;

    // Calculate dispersion relation terms
    this->_calculate_dispersion();
}


template<int NumKn>
void    WaveDispersionFO<NumKn>::update_full(
                                                cusfloat    w0_i, 
                                                cusfloat    water_depth_i, 
                                                cusfloat    grav_acc_i
                                            )
{
    // Update real wave dispersion parameters
    this->update( w0_i, water_depth_i, grav_acc_i );

    // Update imaginary wave dispersion parameters
    this->_calculate_john_terms( );
}



template<int NumKn>
WaveDispersionFO<NumKn>::WaveDispersionFO(
                                            cusfloat    w0_i, 
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
    this->update( w0_i, water_depth_i, grav_acc_i );
}