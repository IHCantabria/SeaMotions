
#ifndef __mesh_gp_hpp
#define __mesh_gp_hpp

// Include local modules
#include "../containers/panel_geom.hpp"
#include "../containers/source_node.hpp"
#include "mesh.hpp"


struct MeshGroup
{
public:
    // Define class attributes
    Mesh**          meshes              = nullptr;
    int             meshes_np           = 0;
    PanelGeom**     panels              = nullptr;
    int*            panels_np           = nullptr;
    int*            panels_cnp          = nullptr;
    int             panels_tnp          = 0;
    SourceNode**    source_nodes        = nullptr;
    int*            source_nodes_np     = nullptr;
    int*            source_nodes_cnp    = nullptr;
    int             source_nodes_tnp    = 0;

    // Define class constructor and destructor
    MeshGroup(
                    Mesh**  meshes,
                    int     mesh_np
                );

    ~MeshGroup(
                    void
                );
    
};

#endif