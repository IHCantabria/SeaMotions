
// Include local modules
#include "gwfdx_interface.hpp"
#include "../math/shape_functions.hpp"


GWFDxInterface::GWFDxInterface( 
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


cuscomplex  GWFDxInterface::operator()(
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
    
    cuscomplex  dG_dX( 0.0, 0.0 );
    if ( R < GREEN_ZEROTH_DR )
    {
        // Define
        cusfloat    x0      =   x - GREEN_DR_EPS / 2.0;
        cusfloat    x1      =   x + GREEN_DR_EPS / 2.0;

        cuscomplex  G0      =   G_integral_wave(
                                                    this->_field_point_j[0],
                                                    this->_field_point_j[1],
                                                    this->_field_point_j[2],
                                                    x0,
                                                    y,
                                                    0.0,
                                                    this->_water_depth,
                                                    *(this->_wave_data),
                                                    *(this->_integrals_db)
                                                );

        cuscomplex  G1      =   G_integral_wave(
                                                    this->_field_point_j[0],
                                                    this->_field_point_j[1],
                                                    this->_field_point_j[2],
                                                    x1,
                                                    y,
                                                    0.0,
                                                    this->_water_depth,
                                                    *(this->_wave_data),
                                                    *(this->_integrals_db)
                                                );
        
        dG_dX   =   ( G1 - G0 ) / GREEN_DR_EPS;
    }
    else
    {
        // Calculate Green function derivatives
        cuscomplex  dG_dR   = G_integral_wave_dr(
                                                    R,
                                                    this->_field_point_j[2],
                                                    z,
                                                    this->_water_depth,
                                                    *(this->_wave_data),
                                                    *(this->_integrals_db)
                                                );
        
        // Calculate X and Y cartesian coordinates derivatives
        cusfloat    dX      = this->_field_point_j[0] - x;
                    dG_dX   = dG_dR * dX / R;
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