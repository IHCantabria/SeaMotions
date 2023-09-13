
// Include local modules
#include "source_node.hpp"


SourceNode::SourceNode( 
                            PanelGeom*  panel_in,
                            int         poly_order_in,
                            int         p_order_in,
                            int         q_order_in,
                            cusfloat*   position_in,
                            cusfloat*   normal_vec_in
                        )
{
    // Storage class attributes
    this->normal_vec    = normal_vec_in;
    this->p_order       = p_order_in;
    this->q_order       = q_order_in;
    this->panel         = panel_in;
    this->poly_oder     = poly_order_in;
    this->position      = position_in;
    
}