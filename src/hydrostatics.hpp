
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

// Include local modules
#include "./containers/mpi_config.hpp"
#include "./inout/input.hpp"
#include "./mesh/mesh.hpp"


struct Hydrostatics
{
private:
    // Declare class methods
    void _calculate( 
                        Mesh*       mesh,
                        MpiConfig*  mpi_config
                    );

public:
    // Declare class attributes
    cusfloat    bmx            = 0.0;
    cusfloat    bmy            = 0.0;
    cusfloat    cob[3]         = { 0.0, 0.0, 0.0 };
    cusfloat    cog[3]         = { 0.0, 0.0, 0.0 };
    cusfloat    displacement   = 0.0;
    cusfloat    grav_acc       = 0.0;
    cusfloat    gmx            = 0.0;
    cusfloat    gmy            = 0.0;
    cusfloat    hydstiffmat[36];
    cusfloat    kb             = 0.0;
    cusfloat    kg             = 0.0;
    cusfloat    kmx            = 0.0;
    cusfloat    kmy            = 0.0;
    cusfloat    mass           = 0.0;
    cusfloat    rad_inertia[3] = { 0.0, 0.0, 0.0 };
    cusfloat    rho_water      = 0.0;
    cusfloat    volume         = 0.0;
    cusfloat    wl_area        = 0.0;
    cusfloat    wl_area_cog[3] = { 0.0, 0.0, 0.0 };
    cusfloat    wl_area_ixx    = 0.0;
    cusfloat    wl_area_ixy    = 0.0;
    cusfloat    wl_area_iyy    = 0.0;
    cusfloat    wl_area_mx     = 0.0;
    cusfloat    wl_area_my     = 0.0;

    // Declare class constructors and destructor
    Hydrostatics( 
                    Mesh*       mesh,
                    cusfloat    rhow_in,
                    cusfloat    grav_acc_in,
                    cusfloat    mass_in,
                    cusfloat*   cog_in,
                    cusfloat*   rad_inertia_in,
                    MpiConfig*  mpi_config
                );

    // Declare class methods
    void print( void );

};
