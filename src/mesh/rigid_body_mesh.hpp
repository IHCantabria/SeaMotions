
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
    std::string                 _auto_flush_fopath  ;               // Folder path where to storage all the mesh seed outputed by the auto flush functionality
    bool                        _is_auto_flush      = false;        // Switch to state if the auto_flush_mesh system is working for the current instance
    cusfloat                    _cog_backup[3]      ;               // Position of the cog w.r.t to the keel. Back-up value used to refresh cog value after a movement.
    cusfloat*                   _x_backup           = nullptr;      // X node positions back-up. Used to recalculate data when performing affine transformations
    cusfloat*                   _y_backup           = nullptr;      // Y node positions back-up. Used to recalculate data when performing affine transformations
    cusfloat*                   _z_backup           = nullptr;      // Z node positions back-up. Used to recalculate data when performing affine transformations
    bool                        _underwater_checked = false;        // True once check_underwater_panels() has been called at least once
    std::vector<int>            _uw_panel_indices   ;               // Indices (into panels[]) of fully-underwater original panels, cached by check_underwater_panels()

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
     * @brief   Delegate method for class construction
     * 
     * This method allows to have a single implementation for class construction 
     * while having multiple constructor interfaces.
     */
    void    _build(
                            cusfloat*   cog_in,
                            cusfloat    draft
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
    int                         fs_panels_np    = 0;            // Number of panels aded by the free surface refinement process
    int                         uw_elems_np     = 0;            // Number of fully-underwater original panels (updated by check_underwater_panels())

    /* Define class constructor */
    RigidBodyMesh( ) = default;

    RigidBodyMesh( 
                                        std::string         file_path,
                                        std::string         body_name,
                                        cusfloat*           cog_in,
                                        bool                is_fix,
                                        PanelTypeE          panel_type,
                                        cusfloat            draft_in
                    );

    RigidBodyMesh( 
                                        std::string         file_path,
                                        std::string         body_name,
                                        cusfloat*           cog_in,
                                        bool                is_fix,
                                        PanelTypeE          panel_type,
                                        cusfloat            draft_in,
                                        std::string         auto_flush_fopath,
                                        bool                auto_flush = false
                    );

    /**
     * @brief Construct from an existing Mesh (e.g. a pre-loaded mesh_total).
     *
     * Copies all node/panel data from @p src_mesh via the Mesh copy constructor
     * and initialises the rigid-body motion back-up buffers.
     *
     * @param src_mesh  Source mesh to copy geometry from.
     * @param cog_in    Centre of gravity (3 components) used as the rotation pivot.
     * @param draft_in  Nominal draft of the body [m] (used for reference only).
     */
    RigidBodyMesh(
                                        const Mesh&         src_mesh,
                                        cusfloat*           cog_in,
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