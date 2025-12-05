
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

// Include general usage libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/panel_geom.hpp"
#include "../../src/containers/source_node.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/tools.hpp"

// Define panels type



struct Mesh
{
protected:
    // Define class attributes
    cusfloat    _fs_centre_x        = 0.0;
    cusfloat    _fs_centre_y        = 0.0;
    cusfloat    _fs_radius          = 0.0;
    bool        _is_bouding_box     = false;
    bool        _is_fs_centre       = false;
    bool        _is_fs_radius       = false;
    cusfloat    _is_move_f          = 0.0;
    bool        _is_source_nodes    = false;
    int         _valid_elem_type[2] = { 3, 4 };
    int         _valid_elem_type_np = 2;

    // Define class methods
    void    _calculate_bounding_box(
                                        void
                                    );

    void    _calculate_fs_centre(
                                        void
                                );
    
    void    _create_panels(
                                        int                 panel_type,
                                        cusfloat*           cog
                            );
    
    bool    _is_valid_type( 
                                        int                 elem_type 
                           );
    
    void    _joint_meshes(
                                        std::vector<Mesh*>  meshes
                        );

    void    _load_gmsh_mesh( 
                                        std::string         file_path,
                                        std::string         body_name
                            );

    void    _load_poly_mesh( 
                                        std::string         file_path,
                                        std::string         body_name
                           );

    void    _load_simply_mesh( 
                                        std::string         file_path,
                                        std::string         body_name
                            );

    void    _update_panels_properties(
                                        cusfloat*           cog
                                    );

public:
    // Define class attributes
    int             bodies_np       = 1;
    int*            elems           = nullptr;
    int*            elems_type      = nullptr;
    int             elems_np        = 0;
    int             enrl            = 0;
    int             mnpe            = 0;
    std::string     name            = "NoName";
    int             nodes_np        = 0;
    PanelGeom**     panels          = nullptr;
    PanelGeom**     panels_wl       = nullptr;
    int             panels_wl_np    = 0;
    int*            panels_type     = nullptr;
    SourceNode**    source_nodes    = nullptr;
    int             source_nodes_np = 0;
    cusfloat*       x               = nullptr;
    cusfloat        x_max           = 0.0;
    cusfloat        x_min           = 0.0;
    cusfloat*       y               = nullptr;
    cusfloat        y_max           = 0.0;
    cusfloat        y_min           = 0.0;
    cusfloat*       z               = nullptr;
    cusfloat        z_max           = 0.0;
    cusfloat        z_min           = 0.0;

    // Define class constructor and destructor
            Mesh( ) = default;

            Mesh( 
                                                std::string         file_path,
                                                std::string         body_name,
                                                cusfloat*           cog,
                                                bool                is_fix,
                                                int                 panel_type
                );

            Mesh(
                                                std::vector<Mesh*>  meshes,
                                                cusfloat*           cog,
                                                bool                is_fix

                );

    virtual ~Mesh( 
                                                void 
                );

    // Define class methods
            void        calculate_fs_radius(
                                                void
                                            );

            void        define_source_nodes(
                                                int                 poly_order,
                                                cusfloat*           cog
                                        );

            void        detect_pc_points(
                                                cusfloat            wl_det_prec
                                        );

            void        detect_wl_points(
                                                cusfloat           wl_det_prec
                                        );

            cusfloat    get_fs_radius(
                                                void
                                        );
    
            void        get_elem_nodes( 
                                                int                 elem_num, 
                                                int&                npe, 
                                                cusfloat*           xn, 
                                                cusfloat*           yn,
                                                cusfloat*           zn
                                    );

    virtual int        get_elems_np(           
                                                void
                                    ) const;

    virtual PanelGeom*  get_panel(
                                                const int           idx
                                    ) const;

            void        set_all_panels_type(
                                                int                 panel_type
                                        );


            void        set_elements_type(
                                                void
                                        );

    
            void        write( 
                                                std::string         fipath
                                );
    
};
