
// Include local modules
#include "grfdy_interface.hpp"
#include "../math/shape_functions.hpp"


GRFDyInterface::GRFDyInterface( 
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


cuscomplex  GRFDyInterface::operator()( 
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
    cuscomplex  dG_dY( 0.0, 0.0 );
    
    if ( R < 1e-4 )
    {
        cusfloat    eps     = 1e-6;
        cusfloat    y0      = y - eps / 2.0;
        cusfloat    y1      = y + eps / 2.0;

        cuscomplex  G0      = G_integral_steady(
                                                    this->_source_j->position[0],
                                                    this->_source_j->position[1],
                                                    this->_source_j->position[2],
                                                    x,
                                                    y0,
                                                    z,
                                                    this->_water_depth
                                                );
        
        cuscomplex  G1      = G_integral_steady(
                                                    this->_source_j->position[0],
                                                    this->_source_j->position[1],
                                                    this->_source_j->position[2],
                                                    x,
                                                    y1,
                                                    z,
                                                    this->_water_depth
                                                );

                    dG_dY   = ( G1 - G0 ) / eps;
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
        cusfloat    dY      = this->_source_j->position[1] - y;
        cuscomplex  dG_dY   = dG_dR * dY / R;
    }
    
    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source_i->p_order,
                                            this->_source_i->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * dG_dY;
}