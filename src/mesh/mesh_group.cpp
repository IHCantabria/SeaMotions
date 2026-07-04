
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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


int     MeshGroup::get_body_id(  
                                            int panel_id
                                )
{
    // Loop over meshes to find the one that contains the panel
    for ( int i=0; i<this->meshes_np; i++ )
    {
        if ( panel_id >= this->panels_cnp[i] && panel_id < this->panels_cnp[i+1] )
        {
            return i;
        }
    }

    // If no mesh is found, return an error value
    return static_cast<int>(-1);
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
    this->meshes_np     = meshes_np_in;

    this->meshes = new Mesh*[ this->meshes_np ];
    for ( int i = 0; i < this->meshes_np; i++ )
    {
        this->meshes[i] = meshes_in[i];
    }

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
    // NOTE: get_elems_np() and get_panel() are called via virtual dispatch so that
    // derived classes (e.g. RigidBodyMesh) can expose only their underwater panels
    // (fully-submerged originals + free-surface-refined halves) without changing
    // the MeshGroup interface.
    this->panels_cnp[0]         = 0;
    this->panels_raddif_cnp[0]  = 0;
    this->source_nodes_cnp[0]   = 0;
    for ( int i=0; i<this->meshes_np; i++ )
    {
        // Get mesh group panels list dimensions (virtual: honours RigidBodyMesh)
        this->panels_np[i]          = this->meshes[i]->get_elems_np( );
        this->panels_cnp[i+1]       = this->panels_cnp[i] + this->panels_np[i];

        this->panels_raddif_np[i]   = 0;
        for ( int j=0; j<this->meshes[i]->get_elems_np( ); j++ )
        {
            if ( this->meshes[i]->get_panel( j )->type == PanelTypeE::DIFFRAC )
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
    // Use virtual get_panel() so RigidBodyMesh can supply its filtered
    // underwater panel set (fully-submerged + FS-refined panels).
    int start_index = 0;
    for ( int i=0; i<this->meshes_np; i++ )
    {
        // Loop over panels to copy its reference (virtual dispatch).
        // We also tag every panel with the owning body id so the time-domain
        // Duhamel σ history can key per-physical-face entries by
        // (parent_body_id, parent_face_id) — surviving the remesher's
        // index reshuffles at the waterline.
        start_index = this->panels_cnp[i];
        for ( int j=0; j<this->meshes[i]->get_elems_np( ); j++ )
        {
            PanelGeom* p = this->meshes[i]->get_panel( j );
            p->parent_body_id      = i;
            this->panels[start_index+j] = p;
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


void MeshGroup::define_source_nodes( int poly_order )
{
    // Free previously owned SourceNode objects (if this method was called before)
    if ( this->_owns_source_nodes && this->source_nodes != nullptr )
    {
        for ( int i=0; i<this->source_nodes_tnp; i++ )
        {
            delete this->source_nodes[i];
        }
    }
    delete [] this->source_nodes;
    this->source_nodes      = nullptr;
    this->source_nodes_tnp  = 0;
    this->_owns_source_nodes = false;

    if ( poly_order != 0 )
    {
        std::cerr << "MeshGroup::define_source_nodes: only poly_order=0 is currently supported." << std::endl;
        throw std::runtime_error( "MeshGroup::define_source_nodes: unsupported poly_order" );
    }

    // One source node per panel (poly_order == 0): source_nodes[j] <-> panels[j]
    this->source_nodes_tnp = this->panels_tnp;
    this->source_nodes     = new SourceNode*[this->source_nodes_tnp];

    cusfloat* position = nullptr;
    cusfloat* normals  = nullptr;

    for ( int i=0; i<this->panels_tnp; i++ )
    {
        PanelGeom* panel = this->panels[i];
        panel->calculate_source_nodes( poly_order, nullptr );
        panel->get_source_nodes_data( position, normals );
        this->source_nodes[i] = new SourceNode(
                                                    panel,
                                                    poly_order,
                                                    0,
                                                    0,
                                                    position,
                                                    normals
                                               );
    }

    this->_owns_source_nodes = true;
}


        MeshGroup::~MeshGroup(
                                            void
                            )
{
    delete [] this->meshes;
    delete [] this->panels;
    delete [] this->panels_np;
    delete [] this->panels_cnp;
    if ( this->_owns_source_nodes && this->source_nodes != nullptr )
    {
        for ( int i=0; i<this->source_nodes_tnp; i++ )
        {
            delete this->source_nodes[i];
        }
    }
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


void MeshGroup::append_memory_report(
                                        std::vector<MemoryReportEntry>& entries,
                                        const std::string& prefix
                                    ) const
{
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "panels_ptrs" ),
                        ( this->panels != nullptr )
                            ? static_cast<std::size_t>( this->panels_tnp ) * sizeof( PanelGeom* )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "panels_np" ),
                        ( this->panels_np != nullptr )
                            ? static_cast<std::size_t>( this->meshes_np ) * sizeof( int )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "panels_cnp" ),
                        ( this->panels_cnp != nullptr )
                            ? static_cast<std::size_t>( this->meshes_np + 1 ) * sizeof( int )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "panels_raddif_np" ),
                        ( this->panels_raddif_np != nullptr )
                            ? static_cast<std::size_t>( this->meshes_np ) * sizeof( int )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "panels_raddif_cnp" ),
                        ( this->panels_raddif_cnp != nullptr )
                            ? static_cast<std::size_t>( this->meshes_np + 1 ) * sizeof( int )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "source_nodes_ptrs" ),
                        ( this->source_nodes != nullptr )
                            ? static_cast<std::size_t>( this->source_nodes_tnp ) * sizeof( SourceNode* )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "source_nodes_np" ),
                        ( this->source_nodes_np != nullptr )
                            ? static_cast<std::size_t>( this->meshes_np ) * sizeof( int )
                            : 0
                    );
    add_memory_entry(
                        entries,
                        memory_report_path( prefix, "source_nodes_cnp" ),
                        ( this->source_nodes_cnp != nullptr )
                            ? static_cast<std::size_t>( this->meshes_np + 1 ) * sizeof( int )
                            : 0
                    );

    if ( this->_is_wl_points )
    {
        add_memory_entry(
                            entries,
                            memory_report_path( prefix, "panels_wl_ptrs" ),
                            ( this->panels_wl != nullptr )
                                ? static_cast<std::size_t>( this->panels_wl_tnp ) * sizeof( PanelGeom* )
                                : 0
                        );
        add_memory_entry(
                            entries,
                            memory_report_path( prefix, "panels_wl_np" ),
                            ( this->panels_wl_np != nullptr )
                                ? static_cast<std::size_t>( this->meshes_np ) * sizeof( int )
                                : 0
                        );
        add_memory_entry(
                            entries,
                            memory_report_path( prefix, "panels_wl_cnp" ),
                            ( this->panels_wl_cnp != nullptr )
                                ? static_cast<std::size_t>( this->meshes_np + 1 ) * sizeof( int )
                                : 0
                        );
    }

    if ( this->is_panels_mirror )
    {
        add_memory_entry(
                            entries,
                            memory_report_path( prefix, "panels_mirror_ptrs" ),
                            ( this->panels_mirror != nullptr )
                                ? static_cast<std::size_t>( this->panels_tnp ) * sizeof( PanelGeom* )
                                : 0
                        );

        std::size_t mirror_bytes = 0;
        if ( this->panels_mirror != nullptr )
        {
            for ( int i = 0; i < this->panels_tnp; i++ )
            {
                if ( this->panels_mirror[i] != nullptr )
                {
                    mirror_bytes += this->panels_mirror[i]->memory_bytes( );
                }
            }
        }
        add_memory_entry( entries, memory_report_path( prefix, "panels_mirror" ), mirror_bytes );
    }
}