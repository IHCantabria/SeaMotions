
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
    size_t          meshes_np           = 0;
    PanelGeom**     panels              = nullptr;
    size_t*         panels_np           = nullptr;
    size_t*         panels_cnp          = nullptr;
    size_t          panels_tnp          = 0;
    SourceNode**    source_nodes        = nullptr;
    size_t*         source_nodes_np     = nullptr;
    size_t*         source_nodes_cnp    = nullptr;
    size_t          source_nodes_tnp    = 0;

    // Define class constructor and destructor
    MeshGroup(
                    Mesh**  meshes,
                    size_t  mesh_np
                );

    ~MeshGroup(
                    void
                );
    
};

#endif