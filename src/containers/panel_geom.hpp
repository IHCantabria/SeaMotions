
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
#include "../config.hpp"
#include "../tools/memory_report.hpp"



struct PanelGeom
{
private:
    // Define prive class attributes
    bool                    _is_source_nodes        = false;
    cusfloat*               _source_normal_vec      = nullptr;
    cusfloat*               _source_positions       = nullptr;

    // Define private class methods
    void    _read_node_list( 
                                    int         npe,
                                    cusfloat*   x_in,
                                    cusfloat*   y_in,
                                    cusfloat*   z_in
                            );

    void    _read_node_list_mask(  
                                    int         npe,
                                    int*        nodes_pos_in,
                                    bool        use_node_mask,
                                    cusfloat*   x_in,
                                    cusfloat*   y_in,
                                    cusfloat*   z_in
                                );

    void    _copy_trivial_from(
                                    const PanelGeom& other
                                );

    void    _release_source_nodes(
                                    void
                                );

public:
    // Define class attributes
    cusfloat                        area                    = 0.0;
    cusfloat                        body_cog[3]             = { 0.0, 0.0, 0.0 };
    cusfloat                        center[3]               = {0.0, 0.0, 0.0};
    cusfloat                        center_wl[3]            = {0.0, 0.0, 0.0};
    cusfloat                        free_surface_log_int    = 0.0;
    cusfloat                        ext_lid_damp_f          = 0.0;
    cusfloat                        global_to_local_mat[9];
    cusfloat                        gauss_points_global_x[NUM_GP2];
    cusfloat                        gauss_points_global_y[NUM_GP2];
    cusfloat                        gauss_points_global_z[NUM_GP2];
    static constexpr std::size_t    gauss_points_np         = NUM_GP;
    cusfloat                        jac_det_gauss_points[NUM_GP2];
    cusfloat                        is_move_f               = 0.0;
    bool                            is_wl_boundary          = false;
    cusfloat                        len_wl                  = 0;
    cusfloat                        length                  = 0.0;
    cusfloat                        local_to_global_mat[9];
    int                             location_zone           = -999;
    static constexpr int            MAX_PANEL_NODES         = 4;
    cusfloat                        moments_fo[3]           = { 0.0, 0.0, 0.0 };
    cusfloat                        moments_so[3]           = { 0.0, 0.0, 0.0 };
    int                             nodes_pos[MAX_PANEL_NODES];
    cusfloat                        normal_vec[6]           = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    cusfloat                        normal_vec_wl[6]        = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    int                             num_nodes               = 0;
    // Stable identity used by the time-domain Duhamel σ history so that
    // _sigma_hist can be keyed by the underlying physical face rather than
    // by the current BEM-panel index (which is reassigned every time
    // RigidBodyMesh::check_underwater_panels() reclassifies the wetted
    // surface).
    //   parent_body_id  — owning rigid body (set by MeshGroup at assembly).
    //   parent_face_id  — body-local index into RigidBodyMesh::panels[],
    //                     i.e. the original-face id from the input mesh.
    //                     FS-refined sub-panels inherit their parent face id.
    int                             parent_body_id          = -1;
    int                             parent_face_id          = -1;
    cusfloat                        sysref_centre[3]        = { 0.0, 0.0, 0.0 };
    PanelTypeE                      type                    = PanelTypeE::NONE;
    int                             wl_nodes[2]             = { 0, 0 };
    cusfloat                        x[MAX_PANEL_NODES];
    cusfloat                        x_wl[2]                 = { 0.0, 0.0 };
    cusfloat                        xl[MAX_PANEL_NODES];
    cusfloat                        xlc[MAX_PANEL_NODES];
    cusfloat                        y[MAX_PANEL_NODES];
    cusfloat                        y_wl[2]                 = { 0.0, 0.0 };
    cusfloat                        yl[MAX_PANEL_NODES];
    cusfloat                        ylc[MAX_PANEL_NODES];
    cusfloat                        z[MAX_PANEL_NODES];
    cusfloat                        z_wl[2]                 = { 0.0, 0.0 };
    cusfloat                        zl[MAX_PANEL_NODES];
    cusfloat                        zlc[MAX_PANEL_NODES];
    cusfloat                        _x2d[3];
    cusfloat                        _global_pos[3];

    // Define class constructors and destructor
    PanelGeom( ) = default;

    PanelGeom(
                                        int         npe,
                                        cusfloat*   x_in,
                                        cusfloat*   y_in,
                                        cusfloat*   z_in,
                                        int         is_move_f_in,
                                        PanelTypeE  type_in,
                                        cusfloat*   cog,
                                        bool        force_auto_type=true
                );

    PanelGeom(
                                        int         npe,
                                        int*        nodes_pos_in,
                                        bool        use_node_mask,
                                        cusfloat*   x_in,
                                        cusfloat*   y_in,
                                        cusfloat*   z_in,
                                        int         is_move_f_in,
                                        PanelTypeE  type_in,
                                        cusfloat*   cog,
                                        bool        force_auto_type=true
                );
    
    PanelGeom(
                                        const PanelGeom& other
                );

    PanelGeom(
                                        PanelGeom&&     other
                ) noexcept;

    PanelGeom& operator=(
                                        const PanelGeom& other
                        );

    PanelGeom& operator=(
                                        PanelGeom&&     other
                        ) noexcept;
    
    ~PanelGeom( 
                                        void 
                );

    std::size_t memory_bytes( void ) const;

    void append_memory_report(
                                    std::vector<MemoryReportEntry>& entries,
                                    const std::string& prefix
                                ) const;

    // Add method to calculate the geometric propertiess
    void    calcualte_free_surface_singularity( 
                                                    void 
                                                );

    template<int NGP>
    void    calculate_integration_properties( 
                                                    void 
                                            );

    void    calculate_properties(       
                                                    cusfloat*   cog
                                );

    template<int NGP>
    void    check_underwater( 
                                                    void
                            );

    void    calculate_source_nodes(
                                                    int         poly_order,
                                                    cusfloat*   cog
                                    );

    void    get_node_position( 
                                                    int         num_node, 
                                                    cusfloat*   node_pos
                               );

    void    get_panel_xy_proj( 
                                                    PanelGeom*  new_panel 
                            );

    void    get_source_nodes_data(
                                                    cusfloat*&  position,
                                                    cusfloat*&  normals_vec
                                );

    bool    is_john( 
                                                    cusfloat*   position,
                                                    cusfloat    water_depth
                                );

    void    local_coords_from_z_proj(
                                                    cusfloat    x,
                                                    cusfloat    y,
                                                    cusfloat&   xi,
                                                    cusfloat&   eta
                                    );

    void    get_node_local_position(   
                                                    int          num_node, 
                                                    cusfloat*    node_pos
                                   );

    void    get_node_local_position_c(   
                                                    int          num_node, 
                                                    cusfloat*    node_pos
                                    );

    void    initialize(
                                                    cusfloat*   cog,
                                                    bool        force_auto_type
                        );

    int     is_inside(                  
                                                    cusfloat*   field_point
                        );

    void    local_to_global( 
                                                    cusfloat    xi, 
                                                    cusfloat    eta,
                                                    cusfloat*   global_pos
                               );

    void    set_ext_lid_damp_f(
                                                    cusfloat    damp_f
                            );

    void    set_new_properties(
                                                    int         npe,
                                                    cusfloat*   x_in,
                                                    cusfloat*   y_in,
                                                    cusfloat*   z_in
                                );

    void    set_new_properties(
                                                    int         npe,
                                                    cusfloat*   x_in,
                                                    cusfloat*   y_in,
                                                    cusfloat*   z_in,
                                                    cusfloat*   cog_in
                                );

    void    set_new_properties(
                                                    int         npe,
                                                    int*        nodes_pos_in,
                                                    bool        use_node_mask,
                                                    cusfloat*   x_in,
                                                    cusfloat*   y_in,
                                                    cusfloat*   z_in
                                );

    void    set_new_properties(
                                                    int         npe,
                                                    int*        nodes_pos_in,
                                                    bool        use_node_mask,
                                                    cusfloat*   x_in,
                                                    cusfloat*   y_in,
                                                    cusfloat*   z_in,
                                                    cusfloat*   cog_in
                                );

    void    write( 
                                                    std::string finame
                );

    friend std::ostream& operator<< ( std::ostream& os, PanelGeom& panel );
};  

// Include template definitions
#include "panel_geom.txx"