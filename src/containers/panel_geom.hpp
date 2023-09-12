
#ifndef __panel_geom_hpp
#define __panel_geom_hpp

// Include local modules
#include "../config.hpp"



struct PanelGeom
{
private:
    // Define class attributes
    bool                    _is_source_nodes        = false;
    cusfloat*               _source_normal_vec      = nullptr;
    cusfloat*               _source_positions       = nullptr;

public:
    // Define class attributes
    cusfloat                area                    = 0.0;
    cusfloat                center[3]               = {0.0, 0.0, 0.0};
    cusfloat                global_to_local_mat[9];
    cusfloat                length                  = 0.0;
    cusfloat                local_to_global_mat[9];
    static constexpr int    MAX_PANEL_NODES         = 4;
    cusfloat                normal_vec[3]           = { 0.0, 0.0, 0.0 };
    int                     num_nodes               = 0;
    cusfloat                sysref_centre[3]        = { 0.0, 0.0, 0.0 };
    cusfloat                x[MAX_PANEL_NODES];
    cusfloat                xl[MAX_PANEL_NODES];
    cusfloat                y[MAX_PANEL_NODES];
    cusfloat                yl[MAX_PANEL_NODES];
    cusfloat                z[MAX_PANEL_NODES];
    cusfloat                zl[MAX_PANEL_NODES];

    // Define class constructors and destructor
    ~PanelGeom( 
                                        void 
                );

    // Add method to calculate the geometric propertiess
    void    calculate_properties(       
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

    int     is_inside(                  
                                        cusfloat*   field_point
                        );

    void    local_to_global( 
                                        cusfloat    xi, 
                                        cusfloat    eta,
                                        cusfloat*   global_pos
                               );

    void    write( 
                                        std::string finame
                );
};  

#endif