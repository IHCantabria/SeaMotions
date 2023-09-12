
// Include local modules
#include "../math/shape_functions.hpp"
#include "mesh.hpp"


void Mesh::_calculate_bounding_box(
                                        void
                                    )
{
    // Restart bounding box values
    this->x_max = -1e302;
    this->x_min =  1e302;
    this->y_max = -1e302;
    this->y_min =  1e302;
    this->z_max = -1e302;
    this->z_min =  1e302;

    // Loop over mesh nodes to get the coordinates
    // positions of the mesh
    for ( int i=0; i<this->nodes_np; i++ )
    {
        // Check for X Coordinate
        if ( this->x[i] > this->x_max )
        {
            this->x_max = this->x[i];
        }

        if ( this->x[i] < this->x_min )
        {
            this->x_min = this->x[i];
        }

        // Check for Y Coordinate
        if ( this->y[i] > this->y_max )
        {
            this->y_max = this->y[i];
        }

        if ( this->y[i] < this->y_min )
        {
            this->y_min = this->y[i];
        }

        // Check for Z Coordinate
        if ( this->z[i] > this->z_max )
        {
            this->z_max = this->z[i];
        }

        if ( this->z[i] < this->z_min )
        {
            this->z_min = this->z[i];
        }
    }

}


void Mesh::_create_panels(
                            void
                        )
{
    // Create array to stogate the panels
    this->panels = new PanelGeom* [this->elems_np];

    // Loop over elements to create every single panel
    int node_num    = 0;
    int npe         = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        // Create new panel
        this->panels[i] = new PanelGeom;

        // Get nodes per element
        npe = this->elems[i*this->enrl];

        // Get panel vertexes form element list
        this->panels[i]->num_nodes = npe;
        for ( int j=0; j<npe; j++ )
        {
            node_num                = this->elems[(i*this->enrl)+j+1];
            this->panels[i]->x[j]   = this->x[node_num];
            this->panels[i]->y[j]   = this->y[node_num];
            this->panels[i]->z[j]   = this->z[node_num];
        }

        // Calculate panel properties
        this->panels[i]->calculate_properties( );

    }
}


void Mesh::define_source_nodes(
                                    int         poly_order,
                                    cusfloat*   cog
                                )
{
    // Get nodes per element depending on the element type
    // and the order
    int dofs_quad   = dofs_rectangular_region( poly_order );
    int dofs_tri    = dofs_triangular_region( poly_order );

    // Count total number of source nodes
    int sn_np = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        if ( this->panels[i]->num_nodes == 3 )
        {
            sn_np += dofs_tri;
        }
        else if ( this->panels[i]->num_nodes == 4 )
        {
            sn_np += dofs_quad;
        }
        else
        {
            std::cerr << "Number of nodes: " << this->panels[i]->num_nodes;
            std::cerr << " for panel: " << i << " not available." << std::endl;
            throw std::runtime_error( "" );
        }
    }
    this->source_nodes_np   = sn_np;

    // Allocate space for the source nodes vector
    this->source_nodes      = new SourceNode*[ sn_np ];
    this->_is_source_nodes  = true;
    
    // Loop over elements to create source nodes objects
    int         count       = 0;
    cusfloat*   position    = nullptr;
    int         local_count = 0;
    cusfloat*   normals_vec = nullptr;
    cusfloat*   this_pos    = nullptr;
    for ( int i=0; i<this->elems_np; i++ )
    {
        // Calculate source nodes over panel
        this->panels[i]->calculate_source_nodes(    
                                                    poly_order,
                                                    cog
                                                );

        // Get source nodes position over the panel
        this->panels[i]->get_source_nodes_data( 
                                                    position,
                                                    normals_vec
                                                );

        // Loop over polynomials degree
        local_count = 0;
        if ( this->panels[i]->num_nodes == 3 )
        {
            for ( int pi=0; pi<poly_order; pi++ )
            {
                for ( int qi=0; qi<poly_order-pi; qi++ )
                {
                    this->source_nodes[count] = new SourceNode(
                                                                    this->panels[i],
                                                                    poly_order,
                                                                    pi,
                                                                    qi,
                                                                    &(position[3*local_count]),
                                                                    &(this->panels[i]->normal_vec[6*local_count])
                                                                );
                    count++;
                    local_count++;
                }
            }
        }
        else if ( this->panels[i]->num_nodes == 4 )
        {
            for ( int pi=0; pi<poly_order; pi++ )
            {
                for ( int qi=0; qi<poly_order; qi++ )
                {
                    local_count = pi*poly_order+qi;
                    this->source_nodes[count] = new SourceNode(
                                                                    this->panels[i],
                                                                    poly_order,
                                                                    pi,
                                                                    qi,
                                                                    &(position[3*local_count]),
                                                                    &(this->panels[i]->normal_vec[6*local_count])
                                                                );
                    count++;
                }
            }
        }
    }

}


void Mesh::get_elem_nodes( 
                            int         elem_num, 
                            int&        npe, 
                            cusfloat*   xn, 
                            cusfloat*   yn,
                            cusfloat*   zn
                        )
{
    // Get nodes per element
    int elem_off = elem_num*this->enrl;
    npe = this->elems[elem_off];

    // Loop over nodes
    int node_i = 0;
    for ( int i=0; i<npe; i++ )
    {
        // Get node ith of the element
        node_i  = this->elems[elem_off+(i+1)];

        // Get coordinates of the ith node
        xn[i]   = this->x[node_i];
        yn[i]   = this->y[node_i];
        zn[i]   = this->z[node_i];

    }
}


bool Mesh::_is_valid_type( 
                            int elem_type 
                        )
{
    bool is_valid = false;
    for ( int i=0; i<this->valid_elem_type_np; i++ )
    {
        if ( elem_type == this->valid_elem_type[i] )
        {
            is_valid = true;
            break;
        }
    }

    return is_valid;
}


void Mesh::_load_poly_mesh( 
                            std::string file_path 
                        )
{
    // Define auxiliar variable to help in the file parsing
    int                 aux_int;
    std::string         aux_str;
    int                 elem_count = 0;
    int                 elem_valid_count = 0;
    int                 elem_total = 0;
    int                 elem_id = 0;
    std::istringstream  iss;
    std::string         line;
    int                 mnpe = 0;
    int                 node_count = 0;
    int                 node_id = 0;
    int                 npe = 0;
    int                 section_id = 0;

    // Open file unit
    std::ifstream infile( file_path );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Skip header lines
    for ( int i=0; i<2; i++ )
    {
        std::getline( infile, line );
    }

    // Read mesh dimensions
    std::getline( infile, line );
    renew_stream( iss, line );
    iss >> section_id;
    iss >> aux_int;
    iss >> aux_int;
    iss >> aux_int;
    iss >> this->nodes_np;
    iss >> elem_total;

    // Allocate space for the mesh variables
    this->x = generate_empty_vector<cusfloat>( this->nodes_np );
    this->y = generate_empty_vector<cusfloat>( this->nodes_np );
    this->z = generate_empty_vector<cusfloat>( this->nodes_np );

    // Read date
    std::getline( infile, line );

    // Loop over nodes to get the position of each of them
    while ( node_count < this->nodes_np )
    {
        // Read node header
        std::getline( infile, line );
        if (  !is_empty_line( line ) )
        {
            renew_stream( iss, line );
            iss >> aux_int;
            iss >> node_id;

            // Check if the current node is the expected one
            if ( ( node_count+1 ) != node_id )
            {
                std::cerr << "IO ERROR" << std::endl;
                std::cerr << "Polyflow format reader."; 
                std::cerr << "Current node number is not the expected one." << std::endl;
                throw std::runtime_error( "" );
            }

            // Read current node positions
            std::getline( infile, line );
            renew_stream( iss, line );

            iss >> this->x[node_count];
            iss >> this->y[node_count];
            iss >> this->z[node_count];

            // Increment node count after reading the node information
            node_count++;
            
        }

    }

    // Get current file position to return back after the
    // elements table inspection
    int elems_table_pos = infile.tellg( );

    // Loop over elements to load the connectivity matrix
    elem_count          = 0;
    elem_valid_count    = 0;
    while ( elem_count < elem_total )
    {
        // Read element header
        std::getline( infile, line );
        
        if ( !is_empty_line( line ) )
        {
            // Read element number and check for maximum nodes per element
            renew_stream( iss, line );
            iss >> aux_int;
            iss >> elem_id;
            iss >> npe;

            if ( ( elem_count+1 ) != elem_id )
            {
                std::cerr << "IO ERROR" << std::endl;
                std::cerr << "Polyflow format reader."; 
                std::cerr << "Current element number is not the expected one." << std::endl;
                throw std::runtime_error( "" );
            }

            if ( npe > mnpe )
            {
                mnpe = npe;
            }

            // Read element configuration line
            std::getline( infile, line );

            // Read element nodes
            std::getline( infile, line );

            // Update element counter
            if ( this->_is_valid_type( npe ) )
            {
                elem_valid_count++;
            }
            elem_count++;

        }

    }

    // Allocate space for the elements data
    this->elems_np  = elem_valid_count;
    this->mnpe      = mnpe;
    this->enrl      = mnpe + 1;
    this->elems     = generate_empty_vector<int>( this->enrl*this->elems_np );

    // Rewind file to the section after nodes definition
    infile.seekg( elems_table_pos, std::ios::beg );

    // Read elements data
    elem_count          = 0;
    elem_valid_count    = 0;
    while ( elem_count < elem_total )
    {
        std::getline( infile, line );
        if ( !is_empty_line( line ) )
        {
            // Read element number and check for maximum nodes per element
            renew_stream( iss, line );
            iss >> aux_int;
            iss >> elem_id;
            iss >> npe;

            if ( ( elem_count+1 ) != elem_id )
            {
                std::cerr << "IO ERROR" << std::endl;
                std::cerr << "Polyflow format reader."; 
                std::cerr << "Current element number is not the expected one." << std::endl;
                throw std::runtime_error( "" );
            }

            // Read element configuration line
            std::getline( infile, line );

            // Read element nodes
            std::getline( infile, line );

            if ( this->_is_valid_type( npe ) )
            {
                renew_stream( iss, line );
                this->elems[ this->enrl*elem_count ] = npe;
                for ( int i=1; i<npe+1; i++ )
                {
                    iss >> this->elems[ this->enrl*elem_count + i ];
                    this->elems[ this->enrl*elem_count + i ]--;
                }

                // Update valid element counter
                elem_valid_count++;
            }

            // Update element counter
            elem_count++;

        }

    }

    // Close file unit
    infile.close( );
}


Mesh::Mesh( 
                            std::string file_path 
            )
{
    // Load mesh
    this->_load_poly_mesh( file_path );

    // Calculate bounding box of the mesh
    this->_calculate_bounding_box( );

    // Create panels for each element
    this->_create_panels( );
}


Mesh::~Mesh( 
                            void 
            )
{
    // Delete source nodes
    if ( this->_is_source_nodes )
    {
        for ( int i=0; i<this->source_nodes_np; i++ )
        {
            delete this->source_nodes[i];
        }
        delete this->source_nodes;
    }
    
    // Delete panels
    for ( int i=0; i<this->elems_np; i++ )
    {
        delete this->panels[i];
    }
    delete [ ] this->panels;

    // Delete elements
    mkl_free( this->elems );

    // Delete nodal positions
    mkl_free( this->x );
    mkl_free( this->y );
    mkl_free( this->z );

}