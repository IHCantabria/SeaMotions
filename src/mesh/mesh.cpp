
// Include general usage libraries
#include <cassert>

// Include local modules
#include "../math/shape_functions.hpp"
#include "mesh.hpp"
#include "../tools.hpp"


void        Mesh::_calculate_bounding_box(
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

    // Set as calculated magnitude
    this->_is_bouding_box = true;

}


void        Mesh::_calculate_fs_centre(
                                                void
                                        )
{
    // Check if the mesh bounding box is available
    if ( !this->_is_bouding_box )
    {
        this->_calculate_bounding_box( );
    }

    // Calculate free surface mesh centre
    this->_fs_centre_x  = ( this->x_max + this->x_min ) / 2.0;
    this->_fs_centre_y  = ( this->y_max + this->y_min ) / 2.0;


    // Set as a calculated magnitudes
    this->_is_fs_centre = true;
}


void        Mesh::calculate_fs_radius(
                                               void
                                       )
{
    // Calculate mean point
    if ( !this->_is_fs_centre )
    {
        this->_calculate_fs_centre( );
    }

    // Calculate maximum distance from the mean position to get the radius
    cusfloat ri         = 0.0;
    this->_fs_radius    = 0.0;

    for ( int i=0; i<this->nodes_np; i++ )
    {
        // Calculate distance from the FS centre to the 
        // ith node
        ri  = std::sqrt(
                            pow2s( this->x[i] - this->_fs_centre_x )
                            +
                            pow2s( this->y[i] - this->_fs_centre_y )
                        );

        // Check if it is the maximum distance
        if ( ri > this->_fs_radius )
        {
            this->_fs_radius = ri;
        }

    }

    // Set as calculated magnitude
    this->_is_fs_radius = true;

}


void        Mesh::_create_panels(
                                               cusfloat*   cog
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

        // Set panel movility
        this->panels[i]->is_move_f = this->_is_move_f;

        // Calculate panel properties
        this->panels[i]->calculate_properties( cog );

        // Set panel type
        if ( std::abs( this->panels[i]->center[2] ) < FS_SEL_THR )
        {
            this->panels[i]->type = LID_PANEL_CODE;
        }
        else
        {
            this->panels[i]->type = DIFFRAC_PANEL_CODE;
        }

        // Calculate integration properties
        this->panels[i]->calculate_integration_properties<NUM_GP>( );

        if ( this->panels[i]->type == LID_PANEL_CODE )
        {
            this->panels[i]->calcualte_free_surface_singularity( );
        }

    }
}


void        Mesh::define_source_nodes(
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
        if ( poly_order == 0 )
        {
            this->source_nodes[count] = new SourceNode(
                                                            this->panels[i],
                                                            poly_order,
                                                            0,
                                                            0,
                                                            position,
                                                            normals_vec
                                                        );
            count++;
        }
        else if ( this->panels[i]->num_nodes == 3 )
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
                                                                    &(normals_vec[6*local_count])
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
                                                                    &(normals_vec[6*local_count])
                                                                );
                    count++;
                }
            }
        }
    }

}


void        Mesh::detect_pc_points(
                                                cusfloat    wl_det_prec
                                    )
{
    // Check if the free surface radius is calculated if not calculate it
    if ( !this->_is_fs_radius )
    {
        this->calculate_fs_radius( );
    }

    // Loop over elements of the mesh to detect the most distanced ones
    // of the centre of the free surface circle
    int*        elems_wl        = generate_empty_vector<int>( this->elems_np * this->enrl );
    int         global_index    = 0;
    int         node_j          = 0;
    cusfloat    r_diff          = 0.0;
    cusfloat    x_diff          = 0.0;
    cusfloat    y_diff          = 0.0;

    for ( int i=0; i<this->elems_np; i++ )
    {
        global_index = i * this->enrl;
        for ( int j=0; j<this->elems[global_index]; j++ )
        {
            // Get node number
            node_j = this->elems[global_index+1+j];

            // Calculate radius from the centre of the mesh to 
            // the jth node
            x_diff = this->x[node_j] - this->_fs_centre_x;
            y_diff = this->y[node_j] - this->_fs_centre_y;
            r_diff = std::sqrt( pow2s( x_diff ) + pow2s( y_diff ) );

            // Check if the node is over the outter line of the free surface
            // mesh
            if ( std::abs( r_diff - this->_fs_radius ) < wl_det_prec )
            {
                elems_wl[global_index]                          += 1;
                elems_wl[global_index+elems_wl[global_index]]   = node_j;
            }
        }
    }

    // Check for those elements that have two nodes on the water line
    int         idx0    = 0;
    int         idx1    = 0;
    int         count   = 0;
    cusfloat    n_mod   = 0.0;
    PanelGeom*  panel_i;
    cusfloat    xd_n    = 0.0;
    cusfloat    yd_n    = 0.0;

    for ( int i=0; i<this->elems_np; i++ )
    {
        global_index = i * this->enrl;
        if ( elems_wl[global_index] > 2 )
        {
            std::stringstream ss;
            ss << "ERROR - INPUT" << std::endl;
            ss << "Element: " << i << " more than two nodes on the WL for the free surface mesh: " << elems_wl[global_index] << std::endl;
            ss << "Nodes: " << elems_wl[global_index+1];
            for ( int j=2; j<elems_wl[global_index]+1; j++ )
            {
                ss << " - " << elems_wl[global_index+j];
            }
            ss << std::endl;
            throw std::runtime_error( ss.str( ).c_str( ) );
        }

        if ( elems_wl[global_index] == 2 )
        {
            panel_i                 = this->panels[i];
            panel_i->is_wl_boundary = true;
            for ( int j=1; j<elems_wl[global_index]+1; j++ )
            {
                node_j                  = elems_wl[global_index+j];
                panel_i->wl_nodes[j-1]  = node_j;
                panel_i->x_wl[j-1]      = this->x[node_j];
                panel_i->y_wl[j-1]      = this->y[node_j];
                panel_i->center_wl[0]   += this->x[node_j];
                panel_i->center_wl[1]   += this->y[node_j];
            }
            panel_i->center_wl[0]       /= 2.0;
            panel_i->center_wl[1]       /= 2.0;

            xd_n                        = panel_i->center_wl[0] - this->_fs_centre_x;
            yd_n                        = panel_i->center_wl[1] - this->_fs_centre_y;
            n_mod                       = std::sqrt( pow2s( xd_n ) + pow2s( yd_n ) );
            panel_i->normal_vec_wl[0]   = xd_n / n_mod;
            panel_i->normal_vec_wl[1]   = yd_n / n_mod;
            panel_i->normal_vec_wl[2]   = 0.0;

            idx0                        = global_index+1;
            idx1                        = global_index+2;
            panel_i->len_wl             = std::sqrt(
                                                        pow2s( this->x[elems_wl[idx1]] - this->x[elems_wl[idx0]] )
                                                        +
                                                        pow2s( this->y[elems_wl[idx1]] - this->y[elems_wl[idx0]] )
                                                    );

            count++;
        }
    }

    // Create a list to pack all the panels that have boundary over the
    // water line
    this->panels_wl_np  = count;
    this->panels_wl     = new PanelGeom*[count];

    count = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        global_index = i * this->enrl;
        if ( elems_wl[global_index] == 2 )
        {
            this->panels_wl[count] = this->panels[i];
            count++;
        }
    }

    // Delete heap memory allocated in the current method
    mkl_free( elems_wl );

}


void        Mesh::detect_wl_points(
                                               cusfloat    wl_det_prec
                                   )
{
    // Loop over elements to check how many points they
    // have on the WL
    int*    elems_wl        = generate_empty_vector<int>( this->elems_np * this->enrl );
    int     global_index    = 0;
    int     node_j          = 0;

    for ( int i=0; i<this->elems_np; i++ )
    {
        global_index = i * this->enrl;
        for ( int j=0; j<this->elems[global_index]; j++ )
        {
            node_j = this->elems[global_index+1+j];
            if ( 
                    ( this->z[node_j] > -wl_det_prec )
                    &&
                    ( this->z[node_j] < wl_det_prec )
                )
            {
                elems_wl[global_index]                          += 1;
                elems_wl[global_index+elems_wl[global_index]]   = node_j;
            }
        }
    }

    // Check for those elements that have two nodes on the water line
    int         idx0    = 0;
    int         idx1    = 0;
    int         count   = 0;
    PanelGeom*  panel_i;
    for ( int i=0; i<this->elems_np; i++ )
    {
        global_index = i * this->enrl;
        if ( ( elems_wl[global_index] > 2 ) && ( this->panels[i]->type == DIFFRAC_PANEL_CODE ) )
        {
            std::stringstream ss;
            ss << "ERROR - INPUT" << std::endl;
            ss << "Element: " << i << " more than two nodes on the WL: " << elems_wl[global_index] << std::endl;
            ss << "Nodes: " << elems_wl[global_index+1];
            for ( int j=2; j<elems_wl[global_index]+1; j++ )
            {
                ss << " - " << elems_wl[global_index+j];
            }
            ss << std::endl;
            throw std::runtime_error( ss.str( ).c_str( ) );
        }

        if ( elems_wl[global_index] == 2 )
        {
            panel_i                 = this->panels[i];
            panel_i->is_wl_boundary = true;
            for ( int j=1; j<elems_wl[global_index]+1; j++ )
            {
                node_j                  = elems_wl[global_index+j];
                panel_i->wl_nodes[j-1]  = node_j;
                panel_i->x_wl[j-1]      = this->x[node_j];
                panel_i->y_wl[j-1]      = this->y[node_j];
                panel_i->center_wl[0]   += this->x[node_j];
                panel_i->center_wl[1]   += this->y[node_j];
            }
            panel_i->center_wl[0]  /= 2.0;
            panel_i->center_wl[1]  /= 2.0;

            idx0                    = global_index+1;
            idx1                    = global_index+2;
            panel_i->len_wl         = std::sqrt(
                                                    pow2s( this->x[elems_wl[idx1]] - this->x[elems_wl[idx0]] )
                                                    +
                                                    pow2s( this->y[elems_wl[idx1]] - this->y[elems_wl[idx0]] )
                                                );

            count++;
        }
    }

    // Create a list to pack all the panels that have boundary over the
    // water line
    this->panels_wl_np  = count;
    this->panels_wl     = new PanelGeom*[count];

    count = 0;
    for ( int i=0; i<this->elems_np; i++ )
    {
        global_index = i * this->enrl;
        if ( elems_wl[global_index] == 2 )
        {
            this->panels_wl[count] = this->panels[i];
            count++;
        }
    }

    // Delete heap memory allocated in the current method
    mkl_free( elems_wl );
}


void        Mesh::get_elem_nodes( 
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


cusfloat    Mesh::get_fs_radius(
                                        void
                                )
{
    assert( this->_is_fs_radius && "There is no FS radius loaded in Mesh class." );

    return this->_fs_radius;
}


bool        Mesh::_is_valid_type( 
                                               int         elem_type 
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


void        Mesh::_joint_meshes(
                                               std::vector<Mesh*>  meshes
                               )
{
    // Get the cumulative number of elements and nodes
    std::vector<int> elems_np_cum;
    std::vector<int> nodes_np_cum;

    elems_np_cum.push_back( 0 );
    nodes_np_cum.push_back( 0 );

    for ( int i=0; i<static_cast<int>( meshes.size( ) ); i++ )
    {
        elems_np_cum.push_back( 
                                    meshes[i]->elems_np 
                                    +
                                    elems_np_cum[i]
                                );
        nodes_np_cum.push_back( 
                                    meshes[i]->nodes_np
                                    +
                                    nodes_np_cum[i]
                                );
    }

    this->elems_np = elems_np_cum[ meshes.size( ) ];
    this->nodes_np = nodes_np_cum[ meshes.size( ) ];

    // Get maximum nodes per element and it's extended counterpart
    // to allocate space for all the elements in the mesh
    this->enrl = 0;
    this->mnpe = 0;
    for ( int i=0; i<static_cast<int>( meshes.size( ) ); i++ )
    {
        if (  meshes[i]->enrl > this->enrl)
        {
            this->enrl = meshes[i]->enrl;
        }

        if (  meshes[i]->mnpe > this->mnpe)
        {
            this->mnpe = meshes[i]->mnpe;
        }
    }

    // Allocate space for the full elements and nodes lists
    this->elems         = generate_empty_vector<int>( this->elems_np * this->enrl );
    this->panels_type   = generate_empty_vector<int>( this->elems_np );
    this->x             = generate_empty_vector<cusfloat>(  this->nodes_np );
    this->y             = generate_empty_vector<cusfloat>(  this->nodes_np );
    this->z             = generate_empty_vector<cusfloat>(  this->nodes_np );

    // Fill elements and nodes positions lists
    int count_elems     = 0;
    int global_index    = 0;
    int local_index     = 0;
    for ( int i=0; i<static_cast<int>( meshes.size( ) ); i++ )
    {
        // Add elements
        for ( int j=0; j<meshes[i]->elems_np; j++ )
        {
            local_index                         = meshes[i]->enrl * j;
            this->elems[count_elems*this->enrl] = meshes[i]->elems[local_index];
            this->panels_type[count_elems]      = meshes[i]->panels_type[j];

            for ( int k=0; k<meshes[i]->elems[local_index]; k++ )
            {
                this->elems[count_elems*this->enrl+k+1] = meshes[i]->elems[local_index+k+1] + nodes_np_cum[i] ;
            }

            count_elems++;
        }

        // Add node positions
        for ( int j=0; j<meshes[i]->nodes_np; j++ )
        {
            global_index            = nodes_np_cum[i] + j;
            this->x[global_index]   = meshes[i]->x[j];
            this->y[global_index]   = meshes[i]->y[j];
            this->z[global_index]   = meshes[i]->z[j];
        }
    }
}


void        Mesh::_load_poly_mesh( 
                                               std::string file_path,
                                               std::string body_name
                               )
{
    // Define auxiliar variable to help in the file parsing
    int                 a0                  = 0;
    int                 a1                  = 0;
    int                 aux_int;
    std::string         aux_str;
    std::string         _body_name;
    int                 elem_count          = 0;
    int                 elem_valid_count    = 0;
    int                 elem_total          = 0;
    int                 elem_id             = 0;
    int                 header_code         = 0;
    std::istringstream  iss;
    int                 items_np            = 0;
    std::string         line;
    int                 mnpe                = 0;
    int                 node_count          = 0;
    int                 node_id             = 0;
    int                 npe                 = 0;
    int                 section_id          = 0;
    int                 sel_count           = 0;

    // Lower case body name to compare with the 
    // names read from files
    str_to_lower( &body_name );

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

    // Storage maximum elements per node values and 
    // it's extended counterpart to storage the number of 
    // elements per node in the same vector
    this->mnpe      = mnpe;
    this->enrl      = mnpe + 1;

    // Allocate space for the local variables to tempora storage the elements
    // and to storage the selection indexes
    int* _elems     = generate_empty_vector<int>( this->enrl * elem_valid_count );
    int* _sel_elems = generate_empty_vector<int>( elem_valid_count );

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
                _elems[ this->enrl*elem_count ] = npe;
                for ( int i=1; i<npe+1; i++ )
                {
                    iss >> _elems[ this->enrl*elem_count + i ];
                    _elems[ this->enrl*elem_count + i ]--;
                }

                // Update valid element counter
                elem_valid_count++;
            }

            // Update element counter
            elem_count++;

        }

    }

    // Read number of elements describing the target body
    while ( true )
    {
        // Get new line
        std::getline( infile, line );

        // Check if the end of the file has been reached
        if ( infile.eof( ) > 0 )
        {
            break;
        }

        // Check if there is a body header line
        if ( !is_empty_line( line ) )
        {
            // Read element number and check for maximum nodes per element
            renew_stream( iss, line );
            iss >> header_code;
            iss >> items_np;
            iss >> aux_int;

            if ( header_code == 21 )
            {
                // Loop until find a new line with at least one element
                while ( true )
                {
                    // Get new line
                    std::getline( infile, line );

                    // Check if there is a body header line
                    if ( !is_empty_line( line ) )
                    {
                        // Read element number and check for maximum nodes per element
                        renew_stream( iss, line );
                        iss >> _body_name;

                        break;
                    }
                }

                // Check if the body name is equal to the target one
                str_to_lower( &_body_name );
                if ( body_name.compare( _body_name ) == 0 )
                {
                    sel_count = 0;
                    while ( true )
                    {
                        // Get new line
                        std::getline( infile, line );

                        // Check if there is a body header line
                        if ( !is_empty_line( line ) )
                        {
                            // Read element number and check for maximum nodes per element
                            renew_stream( iss, line );
                            iss >> header_code;

                            if ( header_code > 20 )
                            {
                                break;
                            }
                            else
                            {
                                // Renew stream for the current line
                                renew_stream( iss, line );

                                // Loop over line elements
                                while ( iss >> a0 >> a1 )
                                {
                                    _sel_elems[sel_count] = a1-1;
                                    sel_count++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if ( sel_count == 0 )
    {
        std::stringstream ss;
        ss << std::endl;
        ss << "ERROR - INPUT" << std::endl;
        ss << "Not possible to find region with name: '" << body_name << "'";
        ss << " in file: " << file_path << std::endl;
        std::cerr << ss.str( ) << std::endl;
        throw std::runtime_error( ss.str( ).c_str( ) );
    }

    // Allocate space for the connectivity matrix
    this->elems_np  = sel_count;
    this->elems     = generate_empty_vector<int>( this->enrl * this->elems_np );

    // Get selected elements for the target body
    for ( int i=0; i<sel_count; i++ )
    {
        for ( int j=0; j<this->enrl; j++ )
        {
            this->elems[i*this->enrl+j] = _elems[_sel_elems[i]*this->enrl+j];
        }
    }

    // Close file unit
    infile.close( );


    // Deallocate heap memory space associated to the current method
    mkl_free( _elems );
    mkl_free( _sel_elems );
}


void        Mesh::_load_simply_mesh( 
                                               std::string file_path,
                                               std::string body_name
                                    )
{
    // Define auxiliar variable to help in the file parsing
    int                 aux_int;
    std::string         aux_str;
    std::string         _body_name;
    std::istringstream  iss;
    std::string         line;
    int                 npe                 = 0;

    // Lower case body name to compare with the 
    // names read from files
    str_to_lower( &body_name );

    // Storage maximum elements per node values and 
    // it's extended counterpart to storage the number of 
    // elements per node in the same vector
    this->mnpe      = 4;
    this->enrl      = this->mnpe + 1;

    // Open file unit
    std::ifstream infile( file_path );
    CHECK_FILE_UNIT_STATUS( infile, file_path );

    // Read mesh dimensions
    std::getline( infile, line );
    renew_stream( iss, line );
    iss >> aux_str;
    iss >> this->nodes_np;

    std::getline( infile, line );
    renew_stream( iss, line );
    iss >> aux_str;
    iss >> this->elems_np;

    // Allocate space for the mesh variables
    this->x         = generate_empty_vector<cusfloat>( this->nodes_np );
    this->y         = generate_empty_vector<cusfloat>( this->nodes_np );
    this->z         = generate_empty_vector<cusfloat>( this->nodes_np );
    this->elems     = generate_empty_vector<int>( this->enrl * this->elems_np );

    // Loop over number of nodes to read them from file
    std::getline( infile, line );
    for ( int i=0; i<this->nodes_np; i++ )
    {
        // Read line from disk
        std::getline( infile, line );
        renew_stream( iss, line );

        // Process line
        iss >> aux_int;
        iss >> aux_int;
        iss >> aux_int;
        iss >> this->x[i];
        iss >> this->y[i];
        iss >> this->z[i];
    }

    // Loop over number of elements to read them from file
    std::getline( infile, line );
    for ( int i=0; i<this->elems_np; i++ )
    {
        // Read line from disk
        std::getline( infile, line );
        renew_stream( iss, line );

        // Process line
        iss >> aux_int;
        iss >> npe;
        this->elems[this->enrl*i+0] = npe;

        for ( int j=0; j<npe; j++ )
        {
            iss >> this->elems[this->enrl*i+1+j];
            this->elems[this->enrl*i+1+j] -= 1;
        }
    }

    // Close file unit
    infile.close( );

}


Mesh::Mesh( 
                                        std::string file_path,
                                        std::string body_name,
                                        cusfloat*   cog,
                                        bool        is_fix,
                                        int         panel_type
            )
{
    // Storage the required input attributes
    this->_is_move_f = static_cast<cusfloat>( !is_fix );

    // Load mesh
    std::string file_ext = get_fipath_extension( file_path );

    if ( file_ext.compare( ".poly" ) == 0 )
    {
        this->_load_poly_mesh( file_path, body_name );
    }
    else if ( file_ext.compare( ".symplymesh.dat" ) )
    {
        this->_load_simply_mesh( file_path, body_name );
    }
    else
    {
        std::cerr << "ERROR - Mesh file extension: " << file_ext << " is not valid." << std::endl; 
        throw std::runtime_error( "Mesh file extension is not valid." );
    }

    // Generate vector with the panels type
    this->set_all_panels_type( panel_type );

    // Calculate bounding box of the mesh
    this->_calculate_bounding_box( );

    // Create panels for each element
    this->_create_panels( cog );
    
}


Mesh::Mesh(
                                        std::vector<Mesh*>  meshes,
                                        cusfloat*           cog,
                                        bool                is_fix
            )
{
    // Storage the required input attributes
    this->_is_move_f = static_cast<cusfloat>( !is_fix );

    // Joint input meshes in a single one
    this->_joint_meshes( meshes );

    // Calculate bounding box of the mesh
    this->_calculate_bounding_box( );

    // Create panels for each element
    this->_create_panels( cog );
}


Mesh::~Mesh( 
                                        void 
            )
{
    std::cout << "Deleting mesh object..." << std::endl;
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
    delete [ ] this->panels_wl;

    // Delete elements
    mkl_free( this->elems );

    // Delete panel type
    mkl_free( this->panels_type );

    // Delete nodal positions
    mkl_free( this->x );
    mkl_free( this->y );
    mkl_free( this->z );

}


void        Mesh::set_all_panels_type(
                                               int panel_type
                                       )
{
    this->panels_type = generate_empty_vector<int>( this->elems_np );
    for ( int i=0; i<this->elems_np; i++ )
    {
        this->panels_type[i] = panel_type;
    }
}