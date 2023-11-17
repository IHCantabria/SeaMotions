
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


GWFDxInterface::~GWFDxInterface( void )
{
    this->_clear_heap( );
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
                                        pow2s( this->_source_j->position[0] - x )
                                        +
                                        pow2s( this->_source_j->position[1] - y )
                                    );
    
    if ( R < ZEROTH_EPS )
    {
        R = ZEROTH_EPS;
    }

    // Calculate Green function derivatives
    cuscomplex  dG_dR   = G_integral_wave_dr(
                                                R,
                                                this->_source_j->position[2],
                                                z,
                                                this->_water_depth,
                                                *(this->_wave_data),
                                                *(this->_integrals_db)
                                            );

    // Calculate X and Y cartesian coordinates derivatives
    cusfloat    dX      = this->_source_j->position[0] - x;
    cuscomplex  dG_dX   = dG_dR * dX / R;

    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source_i->p_order,
                                            this->_source_i->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * dG_dX * this->_source_j->normal_vec[0];
}