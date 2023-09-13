
// Include local modules
#include "mesh_group.hpp"


MeshGroup::MeshGroup(
                        Mesh**  meshes_in,
                        size_t  meshes_np_in
                    )
{
    // Storage necessary input arguments into class
    // attributes
    this->meshes        = meshes_in;
    this->meshes_np     = meshes_np_in;

    // Allocate space for dimensions vectors
    this->panels_np         = new size_t[this->meshes_np];
    this->panels_cnp        = new size_t[this->meshes_np+1];
    this->source_nodes_np   = new size_t[this->meshes_np];
    this->source_nodes_cnp  = new size_t[this->meshes_np+1];
    
    // Loop over meshes to have their dimension
    this->panels_cnp[0]     = 0;
    this->source_nodes_cnp[0]    = 0;
    for ( int i=0; i<this->meshes_np; i++ )
    {
        // Get mesh group panels list dimensionts
        this->panels_np[i]  = this->meshes[i]->elems_np;
        this->panels_cnp[i] = this->panels_cnp[i-1] + this->panels_np[i];

        // Get mesh group source nodes list dimensions
        this->source_nodes_np[i]  = this->meshes[i]->source_nodes_np;
        this->source_nodes_cnp[i] = this->source_nodes_cnp[i-1] + this->source_nodes_np[i];
    }

    this->panels_tnp        = this->panels_cnp[this->meshes_np];
    this->source_nodes_tnp  = this->source_nodes_cnp[this->meshes_np];

    // Allocate space to have a continium list of panels and
    // source nodes
    this->panels        = new PanelGeom*[this->panels_tnp];
    this->source_nodes  = new SourceNode*[this->source_nodes_tnp];
    
    // Loop over meshes to copy all the references into the new vectors
    int start_index = 0;
    for ( int i=0; i<this->meshes_np; i++ )
    {
        // Loop over panels to copy its reference
        start_index = this->panels_cnp[i];
        for ( int j=0; j<this->meshes[i]->elems_np; j++ )
        {
            this->panels[start_index+j] = this->meshes[i]->panels[j];
        }

        // Loop over source nodes to copy its memory address
        start_index = this->source_nodes_cnp[i];
        for ( int j=0; j<this->meshes[i]->source_nodes_np; j++ )
        {
            this->source_nodes[start_index+j] = this->meshes[i]->source_nodes[j];
        }
    }
}


MeshGroup::~MeshGroup(
                        void
                    )
{
    delete [] this->panels;
    delete [] this->panels_np;
    delete [] this->panels_cnp;
    delete [] this->source_nodes;
    delete [] this->source_nodes_np;
    delete [] this->source_nodes_cnp;
}