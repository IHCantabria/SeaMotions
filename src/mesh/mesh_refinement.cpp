
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
#include "../containers/logger.hpp"
#include "mesh_refinement.hpp"



void    refine_underwater_quadrilateral( 
                                            PanelGeom*                  panel,
                                            std::vector<PanelGeom*>&    fs_panels,
                                            int&                        fs_panels_count,
                                            int&                        last_node_index
                                        )
{
    // Define local variables
    constexpr int       nodes_np    = 4;
    constexpr cusfloat  fs_z        = 0.0;

    // Detect nodes location -> Underwater: 0 | Abovewater: 1
    int nodes_loc[nodes_np] = { 0, 0, 0, 0 };
    int nodes_loc_cum       = 0;
    for ( int i=0; i<nodes_np; i++ )
    {
        if ( panel->z[i] > fs_z )
        {
            nodes_loc[i] = 1;
            nodes_loc_cum++;
        }
    }

    bool triangle_below = false;
    if ( nodes_loc_cum > 2 )
    {
        triangle_below = true;
    }

    // Loop around panel nodes to detect intersections
    int edge_count  = 0;
    int edge_num[2];
    int j           = 0;
    int uw_to_abw   = 0;

    // Get first node from bottom to top
    for ( int i=0; i<nodes_np; i++ )
    {
        // Get forward index
        j = ( i + 1 ) % nodes_np;

        // Check if intersect free surface
        if ( nodes_loc[i] != nodes_loc[j] )
        {
            // Check direction
            if ( 
                    ( nodes_loc[i] == 0 )
                    &&
                    ( nodes_loc[j] == 1 )
                )
            {
                uw_to_abw = i;
            }

        }
    }

    // Get order from first underwater to abovewater
    int nn[ nodes_np ];
    for ( int i=0; i<nodes_np; i++ )
    {
        nn[i] = ( uw_to_abw + i ) % nodes_np;
    }

    // Get consecutive edges
    for ( int i=0; i<nodes_np; i++ )
    {
        // Get forward index
        j = ( i + 1 ) % nodes_np;

        // Check if intersect free surface
        if ( nodes_loc[nn[i]] != nodes_loc[nn[j]] )
        {
            // Check for overpassing the edge count
            if ( edge_count > 1 )
            {
                Logger logger;
                logger.error( "Quadrilateral element has more than two sides intersecting the free surface" );
                throw std::runtime_error( "" );
            }

            // Add new edge
            edge_num[edge_count] = i;
            edge_count++;

        }
    }

    // Check for consecutive sides
    bool is_consecutive = std::abs( edge_num[1] - edge_num[0] ) == 1;
    int  nnf            = 3;
    int  nnb            = 2;

    if ( is_consecutive )
    {
        nnf = 2;
        nnb = 1;
    }

    // Calculate free surface crossings
    cusfloat fsi0[3];
    int      fsi0n  = -1;
    cusfloat fsi1[3];
    int      fsi1n  = -1;
    cusfloat lm     = 0.0;
    int      idx0   = 0;
    int      idx1   = 0;

    idx0    = nn[ ( edge_num[0]     ) % panel->num_nodes ];
    idx1    = nn[ ( edge_num[0] + 1 ) % panel->num_nodes ];
    lm      = -panel->z[ idx0 ] / ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi0[0] = panel->x[ idx0 ] + lm * ( panel->x[ idx1 ] - panel->x[ idx0 ] );
    fsi0[1] = panel->y[ idx0 ] + lm * ( panel->y[ idx1 ] - panel->y[ idx0 ] );
    fsi0[2] = panel->z[ idx0 ] + lm * ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi0n   = last_node_index;
    last_node_index++;

    idx0    = nn[ ( edge_num[1]     ) % panel->num_nodes ];
    idx1    = nn[ ( edge_num[1] + 1 ) % panel->num_nodes ];
    lm      = -panel->z[ idx0 ] / ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi1[0] = panel->x[ idx0 ] + lm * ( panel->x[ idx1 ] - panel->x[ idx0 ] );
    fsi1[1] = panel->y[ idx0 ] + lm * ( panel->y[ idx1 ] - panel->y[ idx0 ] );
    fsi1[2] = panel->z[ idx0 ] + lm * ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi1n   = last_node_index;
    last_node_index++;

    // Compose new nodal arrangements
    cusfloat x_new[ nodes_np ];
    cusfloat y_new[ nodes_np ];
    cusfloat z_new[ nodes_np ];
    int      nodes_pos[ nodes_np ];

    if ( !is_consecutive ) // Crossed edges are not consecutives
    {
        // Storage first node
        x_new[0]        = panel->x[ nn[0] ];
        y_new[0]        = panel->y[ nn[0] ];
        z_new[0]        = panel->z[ nn[0] ];
        nodes_pos[0]    = panel->nodes_pos[ nn[0] ];

        // Storage first free surface intersection
        x_new[1]        = fsi0[0];
        y_new[1]        = fsi0[1];
        z_new[1]        = fsi0[2];
        nodes_pos[1]    = fsi0n;

        // Storage second free surface intersection
        x_new[2]        = fsi1[0];
        y_new[2]        = fsi1[1];
        z_new[2]        = fsi1[2];
        nodes_pos[2]    = fsi1n;

        // Storage fourth node
        x_new[3]        = panel->x[ nn[3] ];
        y_new[3]        = panel->y[ nn[3] ];
        z_new[3]        = panel->z[ nn[3] ];
        nodes_pos[3]    = panel->nodes_pos[ nn[3] ];

        // Add new panel
        fs_panels[ fs_panels_count ]->set_new_properties(
                                                            nodes_np,
                                                            nodes_pos,
                                                            false,
                                                            x_new,
                                                            y_new,
                                                            z_new
                                                        );
        if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
        {
            fs_panels_count++;
        }

    }
    else // Crossed edges are consecutive
    {
        if ( triangle_below )
        {
            // Storage T0 - N0
            x_new[0]        = panel->x[ nn[0] ];
            y_new[0]        = panel->y[ nn[0] ];
            z_new[0]        = panel->z[ nn[0] ];
            nodes_pos[0]    = panel->nodes_pos[ nn[0] ];
    
            // Storage T0 - N1
            x_new[1]        = fsi0[0];
            y_new[1]        = fsi0[1];
            z_new[1]        = fsi0[2];
            nodes_pos[1]    = fsi0n;
    
            // Storage T0 - N2
            x_new[2]        = fsi1[0];
            y_new[2]        = fsi1[1];
            z_new[2]        = fsi1[2];
            nodes_pos[2]    = fsi1n;
    
            // Add new panel for T0
            fs_panels[ fs_panels_count ]->set_new_properties(
                                                                3,
                                                                nodes_pos,
                                                                false,
                                                                x_new,
                                                                y_new,
                                                                z_new
                                                            );
            if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
            {
                fs_panels_count++;
            }
        }
        else
        {
            /* Calculate first triange */
    
            // Storage T0 - N0
            x_new[0]        = panel->x[ nn[0] ];
            y_new[0]        = panel->y[ nn[0] ];
            z_new[0]        = panel->z[ nn[0] ];
            nodes_pos[0]    = panel->nodes_pos[ nn[0] ];
    
            // Storage T0 - N1
            x_new[1]        = fsi0[0];
            y_new[1]        = fsi0[1];
            z_new[1]        = fsi0[2];
            nodes_pos[1]    = fsi0n;
    
            // Storage T0 - N2
            x_new[2]        = panel->x[ nn[3] ];
            y_new[2]        = panel->y[ nn[3] ];
            z_new[2]        = panel->z[ nn[3] ];
            nodes_pos[2]    = panel->nodes_pos[ nn[3] ];
    
            // Add new panel for T0
            fs_panels[ fs_panels_count ]->set_new_properties(
                                                                3,
                                                                nodes_pos,
                                                                false,
                                                                x_new,
                                                                y_new,
                                                                z_new
                                                            );
            if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
            {
                fs_panels_count++;
            }
    
            /* Calculate second triange */
    
            // Storage T1 - N0
            x_new[0]        = fsi0[0];
            y_new[0]        = fsi0[1];
            z_new[0]        = fsi0[2];
            nodes_pos[0]    = fsi0n;

    
            // Storage T1 - N1
            x_new[1]        = fsi1[0];
            y_new[1]        = fsi1[1];
            z_new[1]        = fsi1[2];
            nodes_pos[1]    = fsi1n;
    
            // Storage T1 - N2
            x_new[2]        = panel->x[ nn[3] ];
            y_new[2]        = panel->y[ nn[3] ];
            z_new[2]        = panel->z[ nn[3] ];
            nodes_pos[2]    = panel->nodes_pos[ nn[3] ];
    
            // Add new panel for T1
            fs_panels[ fs_panels_count ]->set_new_properties(
                                                                3,
                                                                nodes_pos,
                                                                false,
                                                                x_new,
                                                                y_new,
                                                                z_new
                                                            );
            if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
            {
                fs_panels_count++;
            }
    
            /* Calculate thrid triange */
    
            // Storage T2 - N0
            x_new[0]        = fsi1[0];
            y_new[0]        = fsi1[1];
            z_new[0]        = fsi1[2];
            nodes_pos[0]    = fsi1n;
    
            // Storage T2 - N1
            x_new[1]        = panel->x[ nn[2] ];
            y_new[1]        = panel->y[ nn[2] ];
            z_new[1]        = panel->z[ nn[2] ];
            nodes_pos[1]    = panel->nodes_pos[ nn[2] ];
    
            // Storage T2 - N2
            x_new[2]        = panel->x[ nn[3] ];
            y_new[2]        = panel->y[ nn[3] ];
            z_new[2]        = panel->z[ nn[3] ];
            nodes_pos[2]    = panel->nodes_pos[ nn[3] ];
    
            // Add new panel for T2
            fs_panels[ fs_panels_count ]->set_new_properties(
                                                                3,
                                                                nodes_pos,
                                                                false,
                                                                x_new,
                                                                y_new,
                                                                z_new
                                                            );
            if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
            {
                fs_panels_count++;
            }

        }

    }

}


void    refine_underwater_triangle( 
                                            PanelGeom*                  panel,
                                            std::vector<PanelGeom*>&    fs_panels,
                                            int&                        fs_panels_count,
                                            int&                        last_node_index
                                    )
{
    // Define local variables
    constexpr int       nodes_np    = 3;
    constexpr cusfloat  fs_z        = 0.0;

    // Detect nodes location -> Underwater: 0 | Abovewater: 1
    int nodes_loc[nodes_np] = { 0, 0, 0 };
    int nodes_loc_cum       = 0;
    for ( int i=0; i<nodes_np; i++ )
    {
        if ( panel->z[i] > fs_z )
        {
            nodes_loc[i] = 1;
            nodes_loc_cum++;
        }
    }

    // Loop around panel nodes to detect intersections
    int edge_count  = 0;
    int edge_num[2];
    int j           = 0;
    int uw_to_abw   = 0;

    // Get first node from bottom to top
    for ( int i=0; i<nodes_np; i++ )
    {
        // Get forward index
        j = ( i + 1 ) % nodes_np;

        // Check if intersect free surface
        if ( nodes_loc[i] != nodes_loc[j] )
        {
            // Check direction
            if ( 
                    ( nodes_loc[i] == 0 )
                    &&
                    ( nodes_loc[j] == 1 )
                )
            {
                uw_to_abw = i;
            }

        }
    }

    // Get order from first underwater to abovewater
    int nn[ nodes_np ];
    for ( int i=0; i<nodes_np; i++ )
    {
        nn[i] = ( uw_to_abw + i ) % nodes_np;
    }

    // Get consecutive edges
    for ( int i=0; i<nodes_np; i++ )
    {
        // Get forward index
        j = ( i + 1 ) % nodes_np;

        // Check if intersect free surface
        if ( nodes_loc[nn[i]] != nodes_loc[nn[j]] )
        {
            // Check for overpassing the edge count
            if ( edge_count > 1 )
            {
                Logger logger;
                logger.error( "Quadrilateral element has more than two sides intersecting the free surface" );
                throw std::runtime_error( "" );
            }

            // Add new edge
            edge_num[edge_count] = i;
            edge_count++;

        }
    }

    // Check for consecutive sides
    bool is_consecutive = nodes_loc_cum == 1;
    int  nnf            = 0;
    int  nnb            = 2;

    if ( is_consecutive )
    {
        nnf = 2;
        nnb = 1;
    }

    // Calculate free surface crossings
    cusfloat fsi0[3];
    int      fsi0n  = -1;
    cusfloat fsi1[3];
    int      fsi1n  = -1;
    cusfloat lm     = 0.0;
    int      idx0   = 0;
    int      idx1   = 0;

    idx0    = nn[ ( edge_num[0]     ) % panel->num_nodes ];
    idx1    = nn[ ( edge_num[0] + 1 ) % panel->num_nodes ];
    lm      = -panel->z[ idx0 ] / ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi0[0] = panel->x[ idx0 ] + lm * ( panel->x[ idx1 ] - panel->x[ idx0 ] );
    fsi0[1] = panel->y[ idx0 ] + lm * ( panel->y[ idx1 ] - panel->y[ idx0 ] );
    fsi0[2] = panel->z[ idx0 ] + lm * ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi0n   = last_node_index;
    last_node_index++;

    idx0    = nn[ ( edge_num[1]     ) % panel->num_nodes ];
    idx1    = nn[ ( edge_num[1] + 1 ) % panel->num_nodes ];
    lm      = -panel->z[ idx0 ] / ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi1[0] = panel->x[ idx0 ] + lm * ( panel->x[ idx1 ] - panel->x[ idx0 ] );
    fsi1[1] = panel->y[ idx0 ] + lm * ( panel->y[ idx1 ] - panel->y[ idx0 ] );
    fsi1[2] = panel->z[ idx0 ] + lm * ( panel->z[ idx1 ] - panel->z[ idx0 ] );
    fsi1n   = last_node_index;
    last_node_index++;

    // Allocate space for new nodes
    cusfloat x_new[ nodes_np ];
    cusfloat y_new[ nodes_np ];
    cusfloat z_new[ nodes_np ];
    int      nodes_pos[ nodes_np ];

    // Re-mesh according to the case
    if ( is_consecutive )
    {
        // Calculate mid-point in between the first
        // and the last nodes
        cusfloat midp[3];
        int      midpn = 0;
        midp[0] = ( panel->x[ nn[0] ] + panel->x[ nn[2] ] ) / 2.0;
        midp[1] = ( panel->y[ nn[0] ] + panel->y[ nn[2] ] ) / 2.0;
        midp[2] = ( panel->z[ nn[0] ] + panel->z[ nn[2] ] ) / 2.0;
        midpn   = last_node_index;
        last_node_index++;
    
        /* Calculate first triange */
    
        // Storage T0 - N0
        x_new[0]        = panel->x[ nn[0] ];
        y_new[0]        = panel->y[ nn[0] ];
        z_new[0]        = panel->z[ nn[0] ];
        nodes_pos[0]    = panel->nodes_pos[ nn[0] ];
    
        // Storage T0 - N1
        x_new[1]        = fsi0[0];
        y_new[1]        = fsi0[1];
        z_new[1]        = fsi0[2];
        nodes_pos[1]    = fsi0n;
    
    
        // Storage T0 - N2
        x_new[2]        = midp[0];
        y_new[2]        = midp[1];
        z_new[2]        = midp[2];
        nodes_pos[2]    = midpn;
    
        // Add new panel for T0
        fs_panels[ fs_panels_count ]->set_new_properties(
                                                            3,
                                                            nodes_pos,
                                                            false,
                                                            x_new,
                                                            y_new,
                                                            z_new
                                                        );
        if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
        {
            fs_panels_count++;
        }
    
        /* Calculate second triange */
    
        // Storage T1 - N0
        x_new[0]        = fsi0[0];
        y_new[0]        = fsi0[1];
        z_new[0]        = fsi0[2];
        nodes_pos[0]    = fsi0n;
    
        // Storage T1 - N1
        x_new[1]        = fsi1[0];
        y_new[1]        = fsi1[1];
        z_new[1]        = fsi1[2];
        nodes_pos[1]    = fsi1n;
    
        // Storage T1 - N2
        x_new[2]        = midp[0];
        y_new[2]        = midp[1];
        z_new[2]        = midp[2];
        nodes_pos[2]    = midpn;
    
        // Add new panel for T1
        fs_panels[ fs_panels_count ]->set_new_properties(
                                                            3,
                                                            nodes_pos,
                                                            false,
                                                            x_new,
                                                            y_new,
                                                            z_new
                                                        );
        if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
        {
            fs_panels_count++;
        }
    
        /* Calculate third triange */
    
        // Storage T2 - N0
        x_new[0]        = fsi1[0];
        y_new[0]        = fsi1[1];
        z_new[0]        = fsi1[2];
        nodes_pos[0]    = fsi1n;
    
        // Storage T2 - N1
        x_new[1]        = panel->x[ nn[2] ];
        y_new[1]        = panel->y[ nn[2] ];
        z_new[1]        = panel->z[ nn[2] ];
        nodes_pos[1]    = panel->nodes_pos[ nn[2] ];
    
        // Storage T2 - N2
        x_new[2]        = midp[0];
        y_new[2]        = midp[1];
        z_new[2]        = midp[2];
        nodes_pos[2]    = midpn;
    
        // Add new panel for T2
        fs_panels[ fs_panels_count ]->set_new_properties(
                                                            3,
                                                            nodes_pos,
                                                            false,
                                                            x_new,
                                                            y_new,
                                                            z_new
                                                        );
        if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
        {
            fs_panels_count++;
        }

    }
    else
    {
        // Storage T1 - N0
        x_new[0]        = fsi0[0];
        y_new[0]        = fsi0[1];
        z_new[0]        = fsi0[2];
        nodes_pos[0]    = fsi0n;
    
        // Storage T1 - N1
        x_new[1]        = fsi1[0];
        y_new[1]        = fsi1[1];
        z_new[1]        = fsi1[2];
        nodes_pos[1]    = fsi1n;
    
        // Storage T1 - N2
        x_new[2]        = panel->x[ nn[0] ];
        y_new[2]        = panel->y[ nn[0] ];
        z_new[2]        = panel->z[ nn[0] ];
        nodes_pos[2]    = panel->nodes_pos[ nn[0] ];
    
        // Add new panel for T1
        fs_panels[ fs_panels_count ]->set_new_properties(
                                                            3,
                                                            nodes_pos,
                                                            false,
                                                            x_new,
                                                            y_new,
                                                            z_new
                                                        );
        if ( fs_panels[ fs_panels_count ]->area > MIN_PANEL_AREA )
        {
            fs_panels_count++;
        }
    }

}