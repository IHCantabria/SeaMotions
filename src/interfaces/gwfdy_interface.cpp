
// Include local modules
#include "gwfdy_interface.hpp"
#include "../math/shape_functions.hpp"


GWFDyInterface::GWFDyInterface( 
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


GWFDyInterface::~GWFDyInterface( void )
{
    this->_clear_heap( );
}


cuscomplex  GWFDyInterface::operator()(
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
    cusfloat    dY      = this->_source_j->position[1] - y;
    cuscomplex  dG_dY   = dG_dR * dY / R;
    
    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source_i->p_order,
                                            this->_source_i->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * dG_dY * this->_source_j->normal_vec[1];
}