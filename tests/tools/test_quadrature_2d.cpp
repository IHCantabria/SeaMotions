
// Include general usage libraries
#include <iostream>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/math/gauss.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/topology.hpp"
#include "../../src/mesh/mesh.hpp"
#include "../../src/tools.hpp"


// Define reference area to check the results
cusfloat REF_AREA = 0.5477;


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