
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
#include "morison_element_def.hpp"


struct BodyDef
{
public:
    // Define class attributes
    cusfloat                        cog[3]                  = { 0.0, 0.0, 0.0 };
    bool                            dof_restrictions[6]     = { false, false, false, false, false, false };  ///< Per-DOF kinematic restrictions (false=free, true=fixed): surge, sway, heave, roll, pitch, yaw
    cusfloat                        ext_lid_damp_f          = 0.0;
    cusfloat                        inertia[6]              = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    bool                            interia_by_rad          = false;
    bool                            is_fix                  = false;
    bool                            is_mesh                 = false;
    bool                            is_mesh_int_lid         = false;
    bool                            is_mesh_fs_qtf          = false;
    bool                            is_mesh_total           = false;
    int                             lid_type                = 0;
    cusfloat                        mass                    = 0.0;
    Mesh*                           mesh                    = nullptr;
    std::string                     mesh_body_name          = "";
    Mesh*                           mesh_ext_lid            = nullptr;
    std::string                     mesh_finame             = "";
    Mesh*                           mesh_fs_qtf             = nullptr;
    Mesh*                           mesh_int_lid            = nullptr;
    int                             mesh_items_np           = 0;
    Mesh*                           mesh_total              = nullptr;
    std::vector<MorisonElementDef>  morison_elements        ;
    std::vector<std::string>        morison_elements_names  ;
    std::size_t                     morison_elements_np     = 0;
    cusfloat                        rad_inertia[3]          = { 0.0, 0.0, 0.0 };
    bool                            use_ext_lid             = false;
    bool                            use_morison_elements    = false;

    // Define class constructor and destructor
    ~BodyDef( void );

    // Define class methods
    void print( void );

};