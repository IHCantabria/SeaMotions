
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


// Include local modules
#include "body_def.hpp"


BodyDef::~BodyDef( void )
{
    if ( this->is_mesh )
    {
        delete this->mesh;
    }

    if ( this->is_mesh_fs_qtf )
    {
        delete this->mesh_fs_qtf;
    }
}


void BodyDef::print( void )
{
    std::cout << "Mesh file name:   " << this->mesh_finame << std::endl;
    std::cout << "Mass:             " << this->mass << std::endl;
    std::cout << "COG_X:            " << this->cog[0] << std::endl;
    std::cout << "COG_Y:            " << this->cog[1] << std::endl;
    std::cout << "COG_Z:            " << this->cog[2] << std::endl;
    std::cout << "RXX:              " << this->rad_inertia[0] << std::endl; 
    std::cout << "RYY:              " << this->rad_inertia[1] << std::endl; 
    std::cout << "RZZ:              " << this->rad_inertia[2] << std::endl; 
}