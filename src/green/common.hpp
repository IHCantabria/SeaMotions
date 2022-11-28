
#ifndef __green_common_hpp
#define __green_common_hpp

#include "../config.hpp"
#include "../containers.hpp"


void calculate_distance_node_field(PanelGeom &panel, cusfloat (&field_point_local)[3], cusfloat* node_fieldp_mod, 
    cusfloat* node_fieldp_dx, cusfloat* node_fieldp_dy, cusfloat* node_fieldp_dz);
void calculate_nodes_distance(PanelGeom &panel, cusfloat* delta_xi, cusfloat* delta_eta);


#endif