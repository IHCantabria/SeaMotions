
// Include local modules
#include "gwfdz_interface.hpp"
#include "../math/shape_functions.hpp"
#include "../math/math_tools.hpp"
#include "../math/math_tools.hpp"


GWFDzInterface::GWFDzInterface( 
                                    SourceNode* source_i,
                                    SourceNode* source_j,
                                    cusfloat    ang_freq,
                                    cusfloat    water_depth,
                                    cusfloat    grav_acc
                            )
{
    // Storage panels memory address, water depth and 
    // gravitational acceleration to use along the class
    this->_grav_acc     = grav_acc;
    this->_source_i     = source_i;
    this->_source_j     = source_j;
    this->_water_depth  = water_depth;

    // Initialize object features
    this->_initialize( ang_freq );
}


cuscomplex  GWFDzInterface::operator()(
                                            cusfloat    xi,
                                            cusfloat    eta,
                                            cusfloat    x,
                                            cusfloat    y,
                                            cusfloat    z
                                        )
{
    // Calculate horizontal radius
    cusfloat    R       = std::sqrt(
                                        pow2s( this->_field_point_j[0] - x )
                                        +
                                        pow2s( this->_field_point_j[1] - y )
                                    );
    
    if ( R < ZEROTH_EPS )
    {
        R = ZEROTH_EPS;
    }

    // Calculate Green function derivatives
    cuscomplex  dG_dZ   = G_integral_wave_dz(
                                                R,
                                                this->_field_point_j[2],
                                                z,
                                                this->_water_depth,
                                                *(this->_wave_data),
                                                *(this->_integrals_db)
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