
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
#include "../../mesh/rigid_body_mesh.hpp"


template<typename T>
class HydrostaticForcesNLin
{
private:
    /* Define class private attributes */ 
    RigidBodyMesh*  _mesh           = nullptr;  // Storage pointer of the mesh object used to represent the floater external surface
    T               _force_interf   = nullptr;  // Storage pointer of the functor interface used to calculate the force over the panel
    cusfloat        _pos_init[6]    ;           // Storage initial position of the mesh. It will be used to refresh mesh state during execution except for heave
    cusfloat        _weight         = 0.0;      // Weight of the floater to be accounted on external force vector

public:
    /* Define class constructor */
    HydrostaticForcesNLin( 
                            RigidBodyMesh*  mesh_in,
                            cusfloat*       pos_in,
                            T               force_interf_in,
                            cusfloat        weight_in
                        )
    {
        // Storage input arguments
        this->_force_interf = force_interf_in;
        this->_mesh         = mesh_in;
        this->_weight       = weight_in;

        copy_vector( 6, pos_in, this->_pos_init );
    }

    /* Define class overloaded operators */

    /**
     * @brief   Method to calculate hydrostatic forces for a given heave displacement
     * 
     * @param   heave   Heave displacement [m]
     */
    cusfloat operator()   ( 
                            cusfloat heave
                        ) const;
    
};