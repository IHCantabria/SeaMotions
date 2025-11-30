
// Include local modules
#include "mesh_group.hpp"


void    MeshGroup::define_mirror_panels(
                                            void
                                        )
{

    // Allocate heap memory space for the mirror panels
    this->panels_mirror = new PanelGeom*[this->panels_tnp];

    // Loop over panels to define the mirror ones. Only 
    // defined for diffraction panels
    cusfloat    cog[3]  = { 0.0, 0.0, 0.0 };
    PanelGeom*  panel_i = nullptr;
    cusfloat    zm[4]   = { 0.0, 0.0, 0.0, 0.0 };

    for ( int i=0; i<this->panels_tnp; i++ )
    {
        panel_i = this->panels[i];
        // Get the z mirrored values
        for ( int j=0; j<panel_i->num_nodes; j++ )
        {
            zm[j] = -panel_i->z[j];
        }

        // Create new panel
        this->panels_mirror[i] = new PanelGeom(
                                                    panel_i->num_nodes,
                                                    panel_i->x,
                                                    panel_i->y,
                                                    zm,
                                                    panel_i->is_move_f,
                                                    panel_i->type,
                                                    cog,
                                                    false
                                                );
    }
}


        MeshGroup::MeshGroup(
                                            Mesh**  meshes_in,
                                            int     meshes_np_in,
                                            bool    is_wl_points
                            )
{
    // Storage necessary input arguments into class
    // attributes
    this->_is_wl_points = is_wl_points;
    this->meshes        = meshes_in;
    this->meshes_np     = meshes_np_in;

    // Allocate space for dimensions vectors
    this->panels_raddif_np  = new int[this->meshes_np];     clear_vector( this->meshes_np, this->panels_raddif_np );
    this->panels_raddif_cnp = new int[this->meshes_np+1];   clear_vector( this->meshes_np+1, this->panels_raddif_cnp );
    this->panels_np         = new int[this->meshes_np];     clear_vector( this->meshes_np, this->panels_np );
    this->panels_cnp        = new int[this->meshes_np+1];   clear_vector( this->meshes_np+1, this->panels_cnp );
    this->source_nodes_np   = new int[this->meshes_np];     clear_vector( this->meshes_np, this->source_nodes_np );
    this->source_nodes_cnp  = new int[this->meshes_np+1];   clear_vector( this->meshes_np+1, this->source_nodes_cnp );

    if ( this->_is_wl_points )
    {
        this->panels_wl_np      = new int[this->meshes_np];     clear_vector( this->meshes_np, this->panels_wl_np );
        this->panels_wl_cnp     = new int[this->meshes_np+1];   clear_vector( this->meshes_np+1, this->panels_wl_cnp );
    }
    
    // Loop over meshes to have their dimension
    this->panels_cnp[0]         = 0;
    this->panels_raddif_cnp[0]  = 0;
    this->source_nodes_cnp[0]   = 0;
    for ( int i=0; i<this->meshes_np; i++ )
    {
        // Get mesh group panels list dimensions
        this->panels_np[i]          = this->meshes[i]->elems_np;
        this->panels_cnp[i+1]       = this->panels_cnp[i] + this->panels_np[i];

        this->panels_raddif_np[i]   = 0;
        for ( int j=0; j<this->meshes[i]->elems_np; j++ )
        {
            if ( this->meshes[i]->panels[j]->type == DIFFRAC_PANEL_CODE )
            {
                this->panels_raddif_np[i] += 1;
            }
        }
        this->panels_raddif_cnp[i+1] = this->panels_raddif_cnp[i] + this->panels_raddif_np[i];

        // Get mesh group panels wl list dimension
        if ( this->_is_wl_points )
        {
            this->panels_wl_np[i]       = this->meshes[i]->panels_wl_np;
            this->panels_wl_cnp[i+1]    = this->panels_wl_cnp[i] + this->panels_wl_np[i];
        }

        // Get mesh group source nodes list dimensions
        this->source_nodes_np[i]    = this->meshes[i]->source_nodes_np;
        this->source_nodes_cnp[i+1] = this->source_nodes_cnp[i] + this->source_nodes_np[i];
    }

    this->panels_tnp        = this->panels_cnp[this->meshes_np];
    this->panels_raddif_tnp = this->panels_raddif_cnp[this->meshes_np];
    this->source_nodes_tnp  = this->source_nodes_cnp[this->meshes_np];
    
    if ( this->_is_wl_points ) 
    {
        this->panels_wl_tnp     = this->panels_wl_cnp[this->meshes_np];
    }

    // Allocate space to have a continium list of panels and
    // source nodes
    this->panels        = new PanelGeom*[this->panels_tnp];
    this->panels_wl     = new PanelGeom*[this->panels_wl_tnp];
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

        // Loop over panels wl to copy its reference
        if ( this->_is_wl_points )
        {
            start_index = this->panels_wl_cnp[i];
            for ( int j=0; j<this->meshes[i]->panels_wl_np; j++ )
            {
                this->panels_wl[start_index+j] = this->meshes[i]->panels_wl[j];
            }
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

    if ( this->is_panels_mirror ) 
    {
        for ( int i=0; i<this->panels_tnp; i++ )
        {
            delete this->panels_mirror[i];
        }
        delete [] this->panels_mirror;
    }

    if ( this->_is_wl_points )
    {
        delete [] this->panels_wl;
        delete [] this->panels_wl_np;
        delete [] this->panels_wl_cnp;
    }
}