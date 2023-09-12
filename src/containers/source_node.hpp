
#ifndef __source_node_hpp
#define __source_node_hpp

// Include local modules
#include "panel_geom.hpp"


struct SourceNode
{
public:
    // Define class attributes
    cusfloat*   normal_vec  = nullptr;
    int         p_order     = 0;
    PanelGeom*  panel       = nullptr;
    int         poly_oder   = 0;
    cusfloat*   position    = nullptr;
    int         q_order     = 0;
    

    // Define class constructors and destructor
    SourceNode(
                    PanelGeom*  panel_in,
                    int         poly_order_in,
                    int         p_order_in,
                    int         q_order_in,
                    cusfloat*   position_in,
                    cusfloat*   normal_vec_in
                );

    // Define class methods

};


#endif