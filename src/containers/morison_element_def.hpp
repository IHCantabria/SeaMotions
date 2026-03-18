
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
#include <array>
#include <fstream>
#include <string>

// Include local modules
#include "../config.hpp"


struct MorisonElementDef
{
private:
    /*** Declare class private attributes  ***/

    /*** Declare class private methods  ***/
    void _parse_input_file( 
                                std::ifstream& file
                            );

    void _parse_input_file( 
                                const std::string& fopath,
                                const std::string& finame
                            );
public:
    /*** Declare class public attributes ***/
    cusfloat                area_axial_ca                       = 0.0;                  // Axial reference area for inertial forces calculation
    cusfloat                area_axial_cd                       = 0.0;                  // Axial reference area for drag calculation
    cusfloat                area_transversal_cd                 = 0.0;                  // Transversal reference area for drag calculation
    cusfloat                area_vertical_cd                    = 0.0;                  // Vertical reference area for drag calculation
    std::string             axis_added_mass_coeff_file          ;                       // Axial added-mass coefficient file
    std::string             axis_drag_coeff_file                ;                       // Axial drag coefficient file
    cusfloat                azimuth_angle                       = 0.0;                  // Azimuth angle (rad)
    cusfloat                cog[3]                              = { 0.0, 0.0, 0.0 };    // Center of gravity of the element (x, y, z)
    cusfloat                declination_angle                   = 0.0;                  // Declination angle (rad)
    std::size_t             divisions_np                        = 0;                    // Number of axial divisions
    cusfloat                kc_length                           = 0.0;                  // Characteristic length used for the Keulegan-Carpenter number calculation
    cusfloat                length                              = 0.0;                  // Element length
    cusfloat                rotation_angle                      = 0.0;                  // Rotation around element axis (rad)
    std::array<cusfloat, 3> start_pos                           = { 0.0, 0.0, 0.0 };    // Element start position (x, y, z)
    std::string             transversal_added_mass_coeff_file   ;                       // Transversal added-mass coefficient file
    std::string             transversal_drag_coeff_file         ;                       // Transversal drag coefficient file
    std::string             vertical_added_mass_coeff_file      ;                       // Vertical added-mass coefficient file
    std::string             vertical_drag_coeff_file            ;                       // Vertical drag coefficient file

    /*** Declare class constructor ***/
    MorisonElementDef( ) = default;

    MorisonElementDef( 
                            std::ifstream& file,
                            cusfloat*      cog
                        );

    MorisonElementDef(
                            const std::string& fopath,
                            const std::string& finame,
                            cusfloat*          cog
                        );

    /*** Declare class copy and move constructors and assignment operators ***/
    MorisonElementDef( const MorisonElementDef& )            = delete;
    MorisonElementDef& operator=( const MorisonElementDef& ) = delete;

    MorisonElementDef( MorisonElementDef&& )                 = default;
    MorisonElementDef& operator=( MorisonElementDef&& )      = default;

};