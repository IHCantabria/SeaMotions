
#ifndef __mesh_hpp
#define __mesh_hpp

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
#define DIFFRAC_PANEL_CODE  0
#define LID_PANEL_CODE      1


struct Mesh
{
private:
    // Define class attributes
    cusfloat    _fs_centre_x        = 0.0;
    cusfloat    _fs_centre_y        = 0.0;
    cusfloat    _fs_radius          = 0.0;
    bool        _is_bouding_box     = false;
    bool        _is_fs_centre       = false;
    bool        _is_fs_radius       = false;
    bool        _is_source_nodes    = false;
    int         valid_elem_type[2]  = { 3, 4 };
    int         valid_elem_type_np  = 2;

    // Define class methods
    void    _calculate_bounding_box(
                                        void
                                    );

void        _calculate_fs_centre(
                                        void
                                );
    
    void    _create_panels(
                                        cusfloat*           cog
                            );
    
    bool    _is_valid_type( 
                                        int                 elem_type 
                           );
    
    void    _joint_meshes(
                                        std::vector<Mesh*>  meshes
                        );

    void    _load_poly_mesh( 
                                        std::string         file_path,
                                        std::string         body_name
                           );

public:
    // Define class attributes
    int             bodies_np       = 1;
    int*            elems           = nullptr;
    int             elems_np        = 0;
    int             enrl            = 0;
    int             mnpe            = 0;
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
                                int                 panel_type
        );

    Mesh(
                                std::vector<Mesh*>  meshes,
                                cusfloat*           cog

        );

    ~Mesh( 
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

    void        set_all_panels_type(
                                        int                 panel_type
                                   );
    
};

#endif