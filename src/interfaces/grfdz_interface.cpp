
// Include local modules
#include "grfdz_interface.hpp"
#include "../math/shape_functions.hpp"


GRFDzInterface::GRFDzInterface( 
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


cuscomplex  GRFDzInterface::operator()( 
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
    cuscomplex  dG_dZ   = G_integral_steady_dz(
                                                    R,
                                                    this->_source_j->position[2],
                                                    z,
                                                    this->_water_depth
                                                );
    
    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source_i->p_order,
                                            this->_source_i->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * dG_dZ;
}