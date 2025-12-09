
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

#pragma once

// Include general usage libraries
#include <string>

// Include local modules
#include "../config.hpp"


struct LoadCondition
{
private:
    /* Declare class private attributes */
    std::string     _fopath;           // Folder path where the load condition files are located
    std::string     _finame;           // File name of the load condition file

    /* Declare class private methods */
    /**
     * @brief   Method to read load condition file
     * 
     * @param   fopath_in       Folder path where the load condition files are located
     * @param   finame_in       File name of the load condition file
     * 
     */
    void   _read_load_condition_file( 
                                            std::string    fopath_in,
                                            std::string    finame_in
                                        );

public:
    /* Declare class public attributes */
    cusfloat        cog[3]          = { 0.0, 0.0, 0.0 };                // Load condition center of gravity [m]
    cusfloat        draft           = 0.0;                              // Load condition draft [m]
    cusfloat        inertia[6]      = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // Load condition inertia tensor [kg*m^2]
    bool            interia_by_rad  = false;                            // Switch to indicate if inertia is defined by radii of gyration
    cusfloat        mass            = 0.0;                              // Load condition mass [kg]
    cusfloat        rad_inertia[3]  = { 0.0, 0.0, 0.0 };                // Load condition radii of gyration [m]
    bool            use_mass        = false;                            // Switch to indicate if mass is used to define the loading condition or it is calculated from the draft

    /* Define class constructor */
    LoadCondition( 
                    std::string     fopath_in,
                    std::string     finame_in
                );

};