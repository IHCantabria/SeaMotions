
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


struct FieldPointsDef
{
private:
    /*** Declare class private attributes  ***/
    bool            is_mesh         = false;

    /*** Declare class private methods  ***/
    void _load_mesh( 
                                std::string fopath
                    );

    void _parse_input_file( 
                                std::string fopath,
                                std::string finame 
                            );

public:
    /*** Declare class public attributes ***/
    Mesh*           mesh            = nullptr;
    std::string     mesh_body_name  = "";
    std::string     mesh_finame     = "";
    bool            out_components  = false;
    bool            out_potential   = false;
    bool            out_pressure    = false;
    bool            out_velocity    = false;

    // Define class constructor and destructor
    FieldPointsDef( 
                    std::string fopath,
                    std::string finame
                );
    
    ~FieldPointsDef( void );

};