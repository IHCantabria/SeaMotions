
// Include local modules
#include "../../src/containers/panel_geom.hpp"
#include "../../src/containers/source_node.hpp"
#include "../../src/interfaces/gwf_interface.hpp"
#include "../../src/math/integration.hpp"


int main( void )
{
    // Define problem properties
    cusfloat    ang_freq    = 2 * PI / 7.0;
    cusfloat    cog[3]      = { 0.0, 0.0, -5.0 };
    cusfloat    grav_acc    = 9.81;
    int         poly_order  = 0;
    cusfloat    water_depth = 1000.0;

    // Define panel
    PanelGeom   panel;

    panel.num_nodes = 4;

    panel.x[0] = -5.0;
    panel.x[1] = -5.0;
    panel.x[2] = -5.0;
    panel.x[3] = -5.0;

    panel.y[0] = -5.0;
    panel.y[1] = -5.0;
    panel.y[2] = 5.0;
    panel.y[3] = 5.0;

    panel.z[0] = -10.0;
    panel.z[1] = 0.0;
    panel.z[2] = 0.0;
    panel.z[3] = -10.0;

    panel.calculate_properties( cog );

    // Calculate source nodes over panel
    panel.calculate_source_nodes(    
                                    poly_order,
                                    cog
                                );

    // Get source nodes position over the panel
    cusfloat* position    = generate_empty_vector<cusfloat>( 3 );
    cusfloat* normals_vec = generate_empty_vector<cusfloat>( 3 );
    panel.get_source_nodes_data( 
                                    position,
                                    normals_vec
                                );

    // Create source node
    SourceNode source_node(
                                &panel,
                                0,
                                0,
                                0,
                                position,
                                normals_vec
                            );

    /// Create Function to integrate potential value
    GWFInterface*   green_interf    = new   GWFInterface(
                                                            &source_node,
                                                            1.0,
                                                            panel.center,
                                                            ang_freq,
                                                            water_depth,
                                                            grav_acc
                                                        );

    auto target_fcn =   [green_interf]
                        ( 
                            cusfloat    xi,
                            cusfloat    eta,
                            cusfloat    X,
                            cusfloat    Y,
                            cusfloat    Z
                        ) -> cuscomplex
                        {
                            return (*green_interf)( xi, eta, X, Y, Z );
                        };

    // Integrate potential
    std::cout << "Adaptive Quadrature" << std::endl;
    cuscomplex  pot_i   = adaptive_quadrature_panel(
                                                        &panel,
                                                        target_fcn,
                                                        0.001,
                                                        10
                                                    );
    std::cout << " --> Done!" << std::endl;

    // Delete heap memory
    mkl_free( position );
    mkl_free( normals_vec );

    return 0;
}