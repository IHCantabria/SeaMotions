
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
#include "panel_geom.hpp"


struct SourceNode
{
public:
    // Define class attributes
    cusfloat*   normal_vec  = nullptr;
    int         p_order     = 0;
    PanelGeom*  panel       = nullptr;
    int         poly_oder   = 0;
    cusfloat*   position    = nullptr;
    int         q_order     = 0;
    

    // Define class constructors and destructor
    SourceNode(
                    PanelGeom*  panel_in,
                    int         poly_order_in,
                    int         p_order_in,
                    int         q_order_in,
                    cusfloat*   position_in,
                    cusfloat*   normal_vec_in
                );

    // Define class methods

};
