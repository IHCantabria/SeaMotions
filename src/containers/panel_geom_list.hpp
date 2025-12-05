
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


struct PanelGeomList
{
    // Define class attributes
    int         panel_np;
    PanelGeom** panels;

    // Define constructors and destructor
    PanelGeomList( ) = default;

    PanelGeomList(
                    int         panel_np_in,
                    PanelGeom** panels_in
                )
    {
        this->panel_np  = panel_np_in;
        this->panels    = panels_in;
    }

    ~PanelGeomList( void )
    {
        for ( int i=0; i<panel_np; i++ )
        {
            delete [] this->panels[i];
        }
        delete [] this->panels;
    }
};
