
// Include local modules
#include "panel_geom.hpp"
#include "../math/gauss_t.hpp"
#include "../math/topology.hpp"


template<int NGP>
void PanelGeom::calculate_integration_properties( void )
{
    cusfloat    gp_global[3]     = { 0.0, 0.0, 0.0 };
    for ( int i=0; i<NGP*NGP; i++ )
    {
        // Get local panel points in global coordinates
        this->local_to_global( 
                                    GaussPointsT<2, NGP>::roots_x[i],
                                    GaussPointsT<2, NGP>::roots_y[i],
                                    gp_global
                                );
        
        this->gauss_points_global_x[i]  = gp_global[0];
        this->gauss_points_global_y[i]  = gp_global[1];
        this->gauss_points_global_z[i]  = gp_global[2];

        // Calculate jacobi determinant for each Gauss point
        this->jac_det_gauss_points[i]   = jacobi_det_2d( 
                                                            this->num_nodes,
                                                            this->xl,
                                                            this->yl,
                                                            GaussPointsT<2, NGP>::roots_x[i],
                                                            GaussPointsT<2, NGP>::roots_y[i]
                                                        );

    }
}