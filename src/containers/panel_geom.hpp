
#ifndef __panel_geom_hpp
#define __panel_geom_hpp

// Include local modules
#include "../config.hpp"


struct PanelGeom
{
    cusfloat area = 0.0;
    cusfloat center[3] = {0.0, 0.0, 0.0};
    cusfloat global_to_local_mat[9];
    cusfloat length;
    cusfloat local_to_global_mat[9];
    int num_nodes = 0;
    static constexpr int MAX_PANEL_NODES = 4;
    cusfloat x[MAX_PANEL_NODES];
    cusfloat xl[MAX_PANEL_NODES];
    cusfloat y[MAX_PANEL_NODES];
    cusfloat yl[MAX_PANEL_NODES];
    cusfloat z[MAX_PANEL_NODES];
    cusfloat zl[MAX_PANEL_NODES];

    // Add method to calculate the geometric propertiess
    void calculate_properties( void );
    void get_node_position( int num_node, cusfloat* node_pos );
    void get_node_local_position( int num_node, cusfloat* node_pos );
    int is_inside( cusfloat* field_point );
};

#endif