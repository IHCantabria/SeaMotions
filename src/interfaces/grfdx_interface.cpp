
// Include local modules
#include "grfdx_interface.hpp"
#include "../math/shape_functions.hpp"


GRFDxInterface::GRFDxInterface( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    water_depth
                                )
{
    // Storage panels memory address, water depth and 
    // gravitational acceleration to use along the class
    this->_source_i     = source_i;
    this->_source_j     = source_j;
    this->_water_depth  = water_depth;

}


cuscomplex  GRFDxInterface::operator()( 
                                        cusfloat xi,
                                        cusfloat eta,
                                        cusfloat x,
                                        cusfloat y,
                                        cusfloat z
                                    )
{
    // Calculate horizontal radius
    cusfloat    R       = std::sqrt(
                                        pow2s( this->_source_j->position[0] - x )
                                        +
                                        pow2s( this->_source_j->position[1] - y )
                                    );

    if ( R < ZEROTH_EPS )
    {
        R = ZEROTH_EPS;
    }

    // Calculate Green function derivatives
    cuscomplex  dG_dX( 0.0, 0.0 );

    if ( R < 1e-4 )
    {
        cusfloat    eps     = 1e-6;
        cusfloat    x0      = x - eps / 2.0;
        cusfloat    x1      = x + eps / 2.0;

        cuscomplex  G0      = G_integral_steady(
                                                    this->_source_j->position[0],
                                                    this->_source_j->position[1],
                                                    this->_source_j->position[2],
                                                    x0,
                                                    y,
                                                    z,
                                                    this->_water_depth
                                                );

        cuscomplex  G1      = G_integral_steady(
                                                    this->_source_j->position[0],
                                                    this->_source_j->position[1],
                                                    this->_source_j->position[2],
                                                    x1,
                                                    y,
                                                    z,
                                                    this->_water_depth
                                                );

                    dG_dX   = ( G1 - G0 ) / eps;
    }
    else
    {
        cuscomplex  dG_dR   = G_integral_steady_dr(
                                                        R,
                                                        this->_source_j->position[2],
                                                        z,
                                                        this->_water_depth
                                                    );

        // Calculate X and Y cartesian coordinates derivatives
        cusfloat    dX      = this->_source_j->position[0] - x;
        cuscomplex  dG_dX   = dG_dR * dX / R;
    }
    
    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source_i->p_order,
                                            this->_source_i->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * dG_dX;
}