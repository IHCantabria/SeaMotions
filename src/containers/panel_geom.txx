
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
#include "panel_geom.hpp"
#include "../math/gauss_t.hpp"
#include "../math/topology.hpp"


template<int NGP>
void PanelGeom::calculate_integration_properties( void )
{
    cusfloat    gp_global[3]     = { 0.0, 0.0, 0.0 };
    for ( int i=0; i<NGP*NGP; i++ )
    {
        // Get local panel points in global coordinates
        this->local_to_global( 
                                    GaussPointsT<2, NGP>::roots_x[i],
                                    GaussPointsT<2, NGP>::roots_y[i],
                                    gp_global
                                );
        
        this->gauss_points_global_x[i]  = gp_global[0];
        this->gauss_points_global_y[i]  = gp_global[1];
        this->gauss_points_global_z[i]  = gp_global[2];

        // Calculate jacobi determinant for each Gauss point
        this->jac_det_gauss_points[i]   = jacobi_det_2d( 
                                                            this->num_nodes,
                                                            this->xl,
                                                            this->yl,
                                                            GaussPointsT<2, NGP>::roots_x[i],
                                                            GaussPointsT<2, NGP>::roots_y[i]
                                                        );

    }
}


template<int NGP>
void PanelGeom::check_underwater( 
                                            void
                                )
{
    // Define local variables
    constexpr cusfloat fs_z = 0.0;

    // Check for points above water surface ( Z = 0 )
    int above_np = 0;
    for ( int i=0; i<this->num_nodes; i++ )
    {
        if ( this->z[i] > fs_z )
        {
            above_np++;
        }
    }

    // Check how many nodes are lying on the FS
    int count_fs = 0;
    for ( int i=0; i<this->num_nodes; i++ )
    {
        if ( std::abs( this->z[i] - fs_z ) < FS_SEL_THR )
        {
            count_fs++;
        }
    }

    // Loop over gauss points to make z=0 
    // in order to perfom hydrostatic 
    // pressure integration correctly
    for ( std::size_t i=0; i<NGP*NGP; i++ )
    {
        if ( this->gauss_points_global_z[i] > fs_z )
        {
            this->gauss_points_global_z[i]  = fs_z;
        }
    }

    // Check zone location of the panel
    // UW: -1 -> underwater
    // UW:  0 -> free surface
    // UW:  1 -> abovewater
    if ( above_np == 0 )
    {
        this->location_zone = -1;
    }
    else if ( above_np == this->num_nodes )
    {
        this->location_zone = 1;
    }
    else
    {
        // Check if the free surface cuts along an edge 
        // of the element
        if ( count_fs == 2 )
        {
            if ( above_np > 0 )
            {
                this->location_zone = 1;
            }
            else
            {
                this->location_zone = -1;
            }
        }
        else
        {
            this->location_zone = 0;
        }
    }

}