
#ifndef __mesh_tools_hpp
#define __mesh_tools_hpp

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/containers/panel_geom_list.hpp"

void    refine_element( 
                            PanelGeom*        panel,
                            PanelGeomList*&   panel_list_obj
                        );

#endif