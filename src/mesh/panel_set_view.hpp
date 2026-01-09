
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
#include "mesh_group.hpp"
#include "../containers/panel_geom.hpp"


struct PanelSetView 
{
    PanelGeom** panels;         // Array of panel geometries
    int*        panels_cnp;     // Cumulative counts per body (size meshes_np+1)
    int         panels_tnp;     // Total panels
};


static inline PanelSetView make_panel_view(
                                                MeshGroup*  gp, 
                                                bool        use_waterline
                                            ) 
{
    if ( use_waterline ) 
    {
        return { gp->panels_wl, gp->panels_wl_cnp, gp->panels_wl_tnp };
    }
    return { gp->panels, gp->panels_cnp, gp->panels_tnp };
}