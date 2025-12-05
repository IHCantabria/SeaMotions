
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
#include "tools.hpp"


void    refine_element( 
                            PanelGeom*        panel,
                            PanelGeomList*&   panel_list_obj
                        )
{
    // Define the number of new elements
    int new_elems_np = 4;

    // Allocate heap memory the list of panels
    PanelGeom** panel_list = new PanelGeom*[new_elems_np];

    // Define midnodes
    cusfloat xm[4] = { 0.0, 0.0, 0.0 };
    cusfloat ym[4] = { 0.0, 0.0, 0.0 };
    cusfloat zm[4] = { 0.0, 0.0, 0.0 };

    int j = 0;
    for ( int i=0; i<panel->num_nodes; i++ )
    {
        // Define forward index
        j = ( i + 1 ) % panel->num_nodes;

        // Calculate midnode in between element nodes
        xm[i] = ( panel->x[i] + panel->x[j] ) / 2.0;
        ym[i] = ( panel->y[i] + panel->y[j] ) / 2.0;
        zm[i] = ( panel->z[i] + panel->z[j] ) / 2.0;
        
    }

    // Refine element
    if ( panel->num_nodes == 3 )
    {
        // Create new element 0
        PanelGeom* panel_0  = new PanelGeom;

        panel_0->num_nodes  = 3;

        panel_0->x[0]       = panel->x[0];
        panel_0->y[0]       = panel->y[0];
        panel_0->z[0]       = panel->z[0];

        panel_0->x[1]       = xm[0];
        panel_0->y[1]       = ym[0];
        panel_0->z[1]       = zm[0];

        panel_0->x[2]       = xm[2];
        panel_0->y[2]       = ym[2];
        panel_0->z[2]       = zm[2];

        panel_list[0]       = panel_0;

        // Create new element 1
        PanelGeom* panel_1  = new PanelGeom;

        panel_1->num_nodes  = 3;
        
        panel_1->x[0]       = xm[0];
        panel_1->y[0]       = ym[0];
        panel_1->z[0]       = zm[0];

        panel_1->x[1]       = panel->x[1];
        panel_1->y[1]       = panel->y[1];
        panel_1->z[1]       = panel->z[1];

        panel_1->x[2]       = xm[1];
        panel_1->y[2]       = ym[1];
        panel_1->z[2]       = zm[1];

        panel_list[1]       = panel_1;

        // Create new element 2
        PanelGeom* panel_2  = new PanelGeom;

        panel_2->num_nodes  = 3;
        
        panel_2->x[0]       = xm[0];
        panel_2->y[0]       = ym[0];
        panel_2->z[0]       = zm[0];

        panel_2->x[1]       = xm[1];
        panel_2->y[1]       = ym[1];
        panel_2->z[1]       = zm[1];

        panel_2->x[2]       = xm[2];
        panel_2->y[2]       = ym[2];
        panel_2->z[2]       = zm[2];

        panel_list[2]       = panel_2;

        // Create new element 3
        PanelGeom* panel_3  = new PanelGeom;

        panel_3->num_nodes  = 3;
        
        panel_3->x[0]       = xm[2];
        panel_3->y[0]       = ym[2];
        panel_3->z[0]       = zm[2];

        panel_3->x[1]       = xm[1];
        panel_3->y[1]       = ym[1];
        panel_3->z[1]       = zm[1];

        panel_3->x[2]       = panel->x[2];
        panel_3->y[2]       = panel->y[2];
        panel_3->z[2]       = panel->z[2];

        panel_list[3]       = panel_3;

    }
    else if ( panel->num_nodes == 4)
    {
        // Create center point
        cusfloat cpx = 0.0;
        cusfloat cpy = 0.0;
        cusfloat cpz = 0.0;
        for ( int i=0; i<panel->num_nodes; i++ )
        {
            cpx += panel->x[i];
            cpy += panel->y[i];
            cpz += panel->z[i];
        }
        cpx /= panel->num_nodes;
        cpy /= panel->num_nodes;
        cpz /= panel->num_nodes;

        // Create new element 0
        PanelGeom* panel_0  = new PanelGeom;

        panel_0->num_nodes  = 4;

        panel_0->x[0]       = panel->x[0];
        panel_0->y[0]       = panel->y[0];
        panel_0->z[0]       = panel->z[0];

        panel_0->x[1]       = xm[0];
        panel_0->y[1]       = ym[0];
        panel_0->z[1]       = zm[0];

        panel_0->x[2]       = cpx;
        panel_0->y[2]       = cpy;
        panel_0->z[2]       = cpz;

        panel_0->x[3]       = xm[3];
        panel_0->y[3]       = ym[3];
        panel_0->z[3]       = zm[3];

        panel_list[0]       = panel_0;

        // Create new element 1
        PanelGeom* panel_1  = new PanelGeom;

        panel_1->num_nodes  = 4;

        panel_1->x[0]       = xm[0];
        panel_1->y[0]       = ym[0];
        panel_1->z[0]       = zm[0];

        panel_1->x[1]       = panel->x[1];
        panel_1->y[1]       = panel->y[1];
        panel_1->z[1]       = panel->z[1];

        panel_1->x[2]       = xm[1];
        panel_1->y[2]       = ym[1];
        panel_1->z[2]       = zm[1];

        panel_1->x[3]       = cpx;
        panel_1->y[3]       = cpy;
        panel_1->z[3]       = cpz;

        panel_list[1]       = panel_1;

        // Create new element 2
        PanelGeom* panel_2  = new PanelGeom;

        panel_2->num_nodes  = 4;

        panel_2->x[0]       = cpx;
        panel_2->y[0]       = cpy;
        panel_2->z[0]       = cpz;

        panel_2->x[1]       = xm[1];
        panel_2->y[1]       = ym[1];
        panel_2->z[1]       = zm[1];

        panel_2->x[2]       = panel->x[2];
        panel_2->y[2]       = panel->y[2];
        panel_2->z[2]       = panel->z[2];

        panel_2->x[3]       = xm[2];
        panel_2->y[3]       = ym[2];
        panel_2->z[3]       = zm[2];

        panel_list[2]       = panel_2;

        // Create new element 3
        PanelGeom* panel_3  = new PanelGeom;

        panel_3->num_nodes  = 4;

        panel_3->x[0]       = xm[3];
        panel_3->y[0]       = ym[3];
        panel_3->z[0]       = zm[3];

        panel_3->x[1]       = cpx;
        panel_3->y[1]       = cpy;
        panel_3->z[1]       = cpz;

        panel_3->x[2]       = xm[2];
        panel_3->y[2]       = ym[2];
        panel_3->z[2]       = zm[2];

        panel_3->x[3]       = panel->x[3];
        panel_3->y[3]       = panel->y[3];
        panel_3->z[3]       = panel->z[3];

        panel_list[3]       = panel_3;

    }

    // Loop over elements to calculate local properties
    for ( int i=0; i<new_elems_np; i++ )
    {
        panel_list[i]->calculate_properties( panel->body_cog );
    }

    // Create panel geom list
    panel_list_obj  = new PanelGeomList( 
                                            new_elems_np,
                                            panel_list
                                        );

}