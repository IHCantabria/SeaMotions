
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

// Include general usage libraries
#include <filesystem>
#include "mpi.h"

// Include local modules
#include "../inout/vtu.hpp"
#include "../math/euler_transforms.hpp"
#include "mesh_operations.hpp"
#include "mesh_refinement.hpp"
#include "stability_mesh.hpp"


void StabilityMesh::check_underwater_panels(
                                                void
                                            )
{
    // Loop over all panel to check the location zone
    for ( int i=0; i<this->elems_np; i++ )
    {
        this->panels[i]->check_underwater<NUM_GP>( );
    }

    // Loop over free surface panels and divide them to have the
    // underwater ones
    this->fs_panels_np = 0;
    int last_node_index = this->nodes_np;
    for ( int i=0; i<this->elems_np; i++ )
    {
        if ( this->panels[i]->location_zone == 0 )
        {
            if ( this->panels[i]->num_nodes == 3 )
            {
                refine_underwater_triangle( 
                                                this->panels[i], 
                                                this->fs_panels, 
                                                this->fs_panels_np, 
                                                last_node_index 
                                            );
            }
            else if ( this->panels[i]->num_nodes == 4 )
            {
                refine_underwater_quadrilateral( 
                                                    this->panels[i], 
                                                    this->fs_panels, 
                                                    this->fs_panels_np, 
                                                    last_node_index 
                                                );
            }
        }
    }

    this->fs_nodes_np = last_node_index - this->nodes_np;

}


int     StabilityMesh::get_elems_np(
                                        void
                                    ) const
{
    return this->elems_np + this->fs_panels_np;
}                        


void StabilityMesh::_initialize( 
                                    void 
                                )
{
    // Allocate space for backup nodes
    this->_x_backup  = generate_empty_vector<cusfloat>( this->nodes_np );
    this->_y_backup  = generate_empty_vector<cusfloat>( this->nodes_np );
    this->_z_backup  = generate_empty_vector<cusfloat>( this->nodes_np );

    // Copy data from nodes to dynamic nodes
    for ( int i=0; i<this->nodes_np; i++ )
    {
        this->_x_backup[i] = this->x[i];
        this->_y_backup[i] = this->y[i];
        this->_z_backup[i] = this->z[i];
    }

    // Resize fs_panels vector
    this->fs_panels.resize( this->elems_np );
    for ( int i=0; i<this->elems_np; i++ )
    {
        this->fs_panels[i] = new PanelGeom;
    }

}


StabilityMesh::StabilityMesh(  
                                std::string         file_path,
                                std::string         body_name,
                                cusfloat*           cog_in,
                                bool                is_fix,
                                int                 panel_type,
                                cusfloat            draft_in
                            ): Mesh(
                                        file_path,
                                        body_name,
                                        cog_in,
                                        is_fix,
                                        panel_type,
                                        false
                                    )
{
    // Storage input arguments
    this->_cog_backup[0]    = cog_in[0];
    this->_cog_backup[1]    = cog_in[1];
    this->_cog_backup[2]    = cog_in[2];
    this->draft             = draft_in;

    // Initialize class
    this->_initialize( );

    // // Check for underwater panels
    // this->check_underwater_panels( );

}


StabilityMesh::~StabilityMesh(
                                void
                            )
{
    // Free memory for dynamic nodes
    mkl_free( this->_x_backup );
    mkl_free( this->_y_backup );
    mkl_free( this->_z_backup );

    // Free memory for free surface panels
    for ( int i=0; i<this->elems_np; i++ )
    {
        delete this->fs_panels[i];
    }
}


void    StabilityMesh::move(
                                                    cusfloat            dx,
                                                    cusfloat            dy,
                                                    cusfloat            dz,
                                                    cusfloat            drx,
                                                    cusfloat            dry,
                                                    cusfloat            drz
                            )
{
    /* 1. Move and rotate to the current position */

    // 1.1 Refresh information from back-up nodes
    copy_vector( this->nodes_np, this->_x_backup, this->x );
    copy_vector( this->nodes_np, this->_y_backup, this->y );
    copy_vector( this->nodes_np, this->_z_backup, this->z );

    // 1.2 Perform rotations about the C.O.G
    mesh_points_rotation_refp( 
                                    this->nodes_np,
                                    this->x,
                                    this->y,
                                    this->z,
                                    this->_cog_backup[0],
                                    this->_cog_backup[1],
                                    this->_cog_backup[2],
                                    drx,
                                    dry,
                                    drz
                                );

    // 1.3 Translate to the current position
    mesh_points_translation( 
                                    this->nodes_np,
                                    this->x,
                                    this->y,
                                    this->z,
                                    dx,
                                    dy,
                                    dz
                            );

    // 1.4 Move the current C.O.G position
    this->cog[0] = this->_cog_backup[0] + dx;
    this->cog[1] = this->_cog_backup[1] + dy;
    this->cog[2] = this->_cog_backup[2] + dz;

    /* 2. Update mesh properties */
    this->_update_panels_properties( this->cog );
    
    /* 3. Reset free surface panels intersections */
    this->_reset_fs_intersect( );

}


PanelGeom* StabilityMesh::get_panel( 
                                                    const int   idx
                                    ) const
{
    PanelGeom*  p = nullptr;
    if ( idx > ( this->elems_np - 1 ) )
    {
        p = this->fs_panels[ idx - this->elems_np ];
    }
    else
    {
        p = this->panels[idx];
    }

    return p;
}


void    StabilityMesh::_reset_fs_intersect(
                                                    void
                                            )
{
    this->fs_nodes_np   = 0;
    this->fs_panels_np  = 0;
}


void    StabilityMesh::write_underwater_panels(
                                                    std::string fopath,
                                                    std::string finame
                                                )
{
    // Calculate number of underwater panels
    int past_uw_elems  = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        if ( this->panels[i]->location_zone < 0 )
        {
            past_uw_elems++;
        }
    }

    std::cout << "past_uw_elems: " << past_uw_elems << std::endl;

    // Allocate resources to storage local variables
    int         total_nodes_np      = this->nodes_np + this->fs_nodes_np;
    int         total_elems_np      = past_uw_elems + this->fs_panels_np;
    cusfloat*   x_total             = generate_empty_vector<cusfloat>( total_nodes_np );
    cusfloat*   y_total             = generate_empty_vector<cusfloat>( total_nodes_np );
    cusfloat*   z_total             = generate_empty_vector<cusfloat>( total_nodes_np );
    int*        total_elems         = generate_empty_vector<int>( this->enrl * total_elems_np );
    int*        total_elems_type    = generate_empty_vector<int>( total_elems_np );

    // Copy mesh nodes
    for ( int i=0; i<this->nodes_np; i++ )
    {
        x_total[i] = this->x[i];
        y_total[i] = this->y[i];
        z_total[i] = this->z[i];
    }

    // Add nodes from new free surface panels
    PanelGeom* fspi     = nullptr;
    int        node_num = 0;
    for ( int i=0; i<this->fs_panels_np; i++ )
    {
        fspi = this->fs_panels[i];
        for ( int j=0; j<this->fs_panels[i]->num_nodes; j++ )
        {
            node_num            = fspi->nodes_pos[ j ];
            if ( node_num > this->nodes_np-1 )
            {
                x_total[ node_num ] = fspi->x[ j ];
                y_total[ node_num ] = fspi->y[ j ];
                z_total[ node_num ] = fspi->z[ j ];
            }
        }
    }

    // Copy mesh elements
    PanelGeom* paneli;
    int offset  = 0;
    int count   = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        paneli = this->panels[i];
        if ( paneli->location_zone < 0 )
        {
            offset = count * this->enrl;
            total_elems[ offset + 0 ] = paneli->num_nodes;
            for ( int j=0; j<paneli->num_nodes; j++ )
            {
                total_elems[ offset + 1 + j ] = paneli->nodes_pos[j];
            }
            count++;
        }
    }

    std::cout << "Panel 2:" << std::endl;
    int id = 2;
    for ( int i=0; i<this->panels[id]->num_nodes; i++ )
    {
        std::cout << " -> " << this->panels[id]->x[i] << " " << this->panels[id]->y[i] << " " << this->panels[id]->z[i] << std::endl;
    }

    // Add extra elements coming from the free surface refinement
    for ( int i=0; i<this->fs_panels_np; i++ )
    {
        offset                  = ( i + past_uw_elems ) * this->enrl;
        total_elems[offset+0]   = this->fs_panels[i]->num_nodes;

        for ( int j=0; j<this->fs_panels[i]->num_nodes; j++ )
        {
            total_elems[offset+1+j] = this->fs_panels[i]->nodes_pos[j];
        }

    }

    // Copy elements type
    count = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        paneli = this->panels[i];
        if ( paneli->location_zone < 0 )
        {
            total_elems_type[ count ] = this->elems_type[i];
            count++;
        }
    }

    // Add new panels type
    for ( int i=0; i<this->fs_panels_np; i++ )
    {
        offset = i + past_uw_elems;

        if ( this->fs_panels[i]->num_nodes == 3 )
        {
            total_elems_type[ offset ] = 5;
        }
        else if ( this->fs_panels[i]->num_nodes == 4 )
        {
            total_elems_type[ offset ] = 9;
        }

    }

    // Get current process ID
    int proc_rank = 0;
    MPI_Comm_rank(
                    MPI_COMM_WORLD,
                    &proc_rank
                );
    if ( proc_rank == MPI_ROOT_PROC_ID )
    {
        // Compose file path
        std::string finame_ext = finame + ".vtu";
        std::filesystem::path dir ( fopath );
        std::filesystem::path file ( finame_ext );
        std::filesystem::path filepath = dir / file;
        
        // Write mesh
        write_vtu_binary_appended( 
                                    filepath.string( ),
                                    total_nodes_np,
                                    x_total,
                                    y_total,
                                    z_total,
                                    total_elems_np,
                                    this->enrl,
                                    total_elems,
                                    total_elems_type
                                );
    }

    // Free local heap memory allocations
    mkl_free( x_total );
    mkl_free( y_total );
    mkl_free( z_total );
    mkl_free( total_elems );
    mkl_free( total_elems_type );

}