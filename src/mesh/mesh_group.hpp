
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
#include "../containers/panel_geom.hpp"
#include "../containers/source_node.hpp"
#include "mesh.hpp"
#include "../tools/memory_report.hpp"


struct MeshGroup
{
private:
    bool            _is_wl_points           = false;
    bool            _owns_source_nodes      = false;    ///< true when source_nodes[] entries were allocated by define_source_nodes()
    
public:
    // Define class attributes
    bool            is_panels_mirror    = false;
    Mesh**          meshes              = nullptr;
    int             meshes_np           = 0;
    PanelGeom**     panels              = nullptr;
    PanelGeom**     panels_mirror       = nullptr;
    int*            panels_raddif_np    = nullptr;
    int*            panels_raddif_cnp   = nullptr;
    int             panels_raddif_tnp   = 0;
    int*            panels_np           = nullptr;
    int*            panels_cnp          = nullptr;
    int             panels_tnp          = 0;
    PanelGeom**     panels_wl           = nullptr;
    int*            panels_wl_np        = nullptr;
    int*            panels_wl_cnp       = nullptr;
    int             panels_wl_tnp       = 0;
    SourceNode**    source_nodes        = nullptr;
    int*            source_nodes_np     = nullptr;
    int*            source_nodes_cnp    = nullptr;
    int             source_nodes_tnp    = 0;

    // Define class constructor and destructor
    MeshGroup(
                    Mesh**  meshes,
                    int     mesh_np,
                    bool    is_wl_points
                );

    ~MeshGroup(
                    void
                );

    void append_memory_report(
                                    std::vector<MemoryReportEntry>& entries,
                                    const std::string& prefix
                                ) const;

    // Define class methods
    void    define_mirror_panels(
                                    void
                                );

    /**
     * @brief Build source nodes from the MeshGroup's own panels[] array.
     *
     * Must be called after the MeshGroup is fully built (and after
     * define_mirror_panels() if mirror panels are needed).  Only
     * poly_order == 0 (one source node per panel at the panel centre)
     * is currently supported.  The resulting source_nodes[] array is
     * 1-to-1 with panels[]: source_nodes[j] corresponds to panels[j].
     *
     * The MeshGroup takes ownership of the created SourceNode objects
     * and frees them on destruction or when this method is called again.
     */
    void    define_source_nodes(
                                    int     poly_order
                                );

    int     get_body_id(  
                                    int panel_id
                                );
    
};
