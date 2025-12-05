
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
#include "../mesh/mesh.hpp"


struct BodyDef
{
public:
    // Define class attributes
    cusfloat        cog[3]          = { 0.0, 0.0, 0.0 };
    cusfloat        inertia[6]      = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bool            interia_by_rad  = false;
    bool            is_fix          = false;
    bool            is_mesh         = false;
    bool            is_mesh_fs_qtf  = false;
    int             lid_type        = 0;
    cusfloat        mass            = 0.0;
    Mesh*           mesh            = nullptr;
    Mesh*           mesh_fs_qtf     = nullptr;
    std::string     mesh_finame     = "";
    std::string     mesh_body_name  = "";
    cusfloat        rad_inertia[3]  = { 0.0, 0.0, 0.0 };

    // Define class constructor and destructor
    ~BodyDef( void );

    // Define class methods
    void print( void );

};