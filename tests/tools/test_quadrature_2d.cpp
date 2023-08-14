
// Include general usage libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/gauss.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/topology.hpp"
#include "../../src/tools.hpp"


// Define reference area to check the results
cusfloat REF_AREA = 0.5477;


struct Mesh
{
private:
    int valid_elem_type[2]  = { 3, 4 };
    int valid_elem_type_np  = 2;

bool is_valid_type( int elem_type )
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

public:
    // Define class attributes
    int*        elems       = nullptr;
    int         elems_np    = 0;
    int         enrl        = 0;
    int         mnpe        = 0;
    int         nodes_np    = 0;
    cusfloat*   x           = nullptr;
    cusfloat*   y           = nullptr;
    cusfloat*   z           = nullptr;

    // Define class constructor and destructor
    Mesh( std::string file_path )
    {
        this->load_poly_mesh( file_path );
    }

    ~Mesh( void )
    {
        // Delete nodal positions
        mkl_free( this->x );
        mkl_free( this->y );
        mkl_free( this->z );
    }

    // Define class methods
    void get_elem_nodes( 
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

    void load_poly_mesh( std::string file_path )
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
                if ( this->is_valid_type( npe ) )
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

                if ( this->is_valid_type( npe ) )
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

};


cusfloat calculate_area( int np, cusfloat* xn, cusfloat* yn )
{
    // Define gauss points for the integration
    const int gp_np = 3;
    cusfloat gp_roots[gp_np], gp_weights[gp_np];
    get_gauss_legendre( gp_np, gp_roots, gp_weights );

    // Loop over gauss points to perform the area integration
    cusfloat int_value = 0.0;
    for ( int i=0; i<gp_np; i++ )
    {
        for ( int j=0; j<gp_np; j++ )
        {
            int_value += gp_weights[i]*gp_weights[j]*jacobi_det_2d( np, xn, yn, gp_roots[i], gp_roots[j] );
        }
    }

    return int_value;
}


void launch_integration( std::string msh_fipath )
{
    // Load mesh
    Mesh msh( msh_fipath );

    // Define local variables
    int         npe     = 0;
    cusfloat    xn[4]   = { 0.0, 0.0, 0.0, 0.0 };
    cusfloat    yn[4]   = { 0.0, 0.0, 0.0, 0.0 };
    cusfloat    zn[4]   = { 0.0, 0.0, 0.0, 0.0 };

    // Loop over mesh elements definition to calculate
    // the cumulative area of all of them
    cusfloat area = 0.0;
    for ( int i=0; i<msh.elems_np; i++ )
    {
        // Get element nodes
        msh.get_elem_nodes( i, npe, xn, yn, zn );

        // Calcualte area of the current element
        area += calculate_area( npe, xn, yn );
    }

    // Compare total area with the reference value
    if ( !assert_scalar_equality( area, REF_AREA, 1e-4 ) )
    {
        std::cerr << "test_quadrature_2d failed!" << std::endl;
        throw std::runtime_error( " " );
    }

}


int main( int argc, char* argv[] )
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 2))
    {
        return 1;
    }

    std::string tri_fipath( argv[1] );
    std::string quad_fipath( argv[2] );

    // Launch test for triangular elements
    launch_integration( tri_fipath );

    // Launch test for quadrilateral elements
    launch_integration( quad_fipath );

    return 0;
}