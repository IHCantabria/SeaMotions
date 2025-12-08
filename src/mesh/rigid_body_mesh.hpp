
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
#include "mesh.hpp"


class RigidBodyMesh: public Mesh
{
private:
    /* Define private attributes */
    int                         _auto_flush_seed    = 0;            // Seed number used to track the number of times that a mesh have been flushed automatically
    cusfloat                    _cog_backup[3]      ;               // Position of the cog w.r.t to the keel. Back-up value used to refresh cog value after a movement.
    cusfloat*                   _x_backup           = nullptr;      // X node positions back-up. Used to recalculate data when performing affine transformations
    cusfloat*                   _y_backup           = nullptr;      // Y node positions back-up. Used to recalculate data when performing affine transformations
    cusfloat*                   _z_backup           = nullptr;      // Z node positions back-up. Used to recalculate data when performing affine transformations

    /* Define private methods */


    /**
     * @brief   This allows to flush current mesh in a .vtu
     * 
     * This method is used for debugging purposes to track how the mesh evolves during
     * iterations or in time evolution.
     * 
     */
    void    _auto_flush( 
                            void
                        );

    /**
     * @brief   Initialization method to have back-up copies of the input mesh.
     */
    void    _initialize( 
                                    void 
                        );
    
    /**
     * @brief Reset panel free surface intersections
     */
    void    _reset_fs_intersect( 
                                    void
                                );

public:
    /* Define public attributes */
    cusfloat                    cog[3];                         // Position of the cog w.r.t to the keel. Position w.r.t reference system after a movement.
    cusfloat                    draft           = 0.0;          // Draft of the floater w.r.t to the keel
    int                         fs_nodes_np     = 0;            // Number of nodes added by the free surface refinement process
    std::vector<PanelGeom*>     fs_panels;                      // New panels generated to refine the mesh around free surface and to cut it cleanly
    int                         fs_panels_np = 0;               // Number of panels aded by the free surface refinement process

    /* Define class constructor */
    RigidBodyMesh( ) = default;

    RigidBodyMesh( 
                                        std::string         file_path,
                                        std::string         body_name,
                                        cusfloat*           cog_in,
                                        bool                is_fix,
                                        int                 panel_type,
                                        cusfloat            draft_in
                    );

    ~RigidBodyMesh( );

    /* Declare class methods */
    void        check_underwater_panels( 
                                                void 
                                    );

    int         get_elems_np(           
                                                void
                            ) const override;

    PanelGeom*  get_panel(           
                                                const int idx
                            ) const override;

    void        move(
                                                cusfloat            x,
                                                cusfloat            y,
                                                cusfloat            z,
                                                cusfloat            rx,
                                                cusfloat            ry,
                                                cusfloat            rz
                    );

    void        write_underwater_panels(
                                                std::string         fopath,
                                                std::string         finame
                                        );
    
};