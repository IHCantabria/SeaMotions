
// Include local modules
#include "grfdn_interface.hpp"
#include "../math/shape_functions.hpp"


GRFDnInterface::GRFDnInterface( 
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


GRFDnInterface::~GRFDnInterface( void )
{
    std::cerr << "Calling GRFDnInterface destructor..." << std::endl;
}


cuscomplex  GRFDnInterface::operator()( 
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
    cuscomplex  dG_dR   = G_integral_steady_dr(
                                                    R,
                                                    this->_source_j->position[2],
                                                    z,
                                                    this->_water_depth
                                                );
    cuscomplex  dG_dZ   = G_integral_steady_dz(
                                                    R,
                                                    this->_source_j->position[2],
                                                    z,
                                                    this->_water_depth
                                                );

    // Calculate X and Y cartesian coordinates derivatives
    cusfloat    dX      = this->_source_j->position[0] - x;
    cuscomplex  dG_dX   = dG_dR * dX / R;
    cusfloat    dY      = this->_source_j->position[1] - y;
    cuscomplex  dG_dY   = dG_dR * dY / R;

    // Calculate normal derivate
    cuscomplex  dG_dn   =   (
                                dG_dX * this->_source_j->normal_vec[0]
                                +
                                dG_dY * this->_source_j->normal_vec[1]
                                +
                                dG_dZ * this->_source_j->normal_vec[2]
                            );
    
    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source_i->p_order,
                                            this->_source_i->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * dG_dn;
}


void    GRFDnInterface::set_field_point_j(
                                            cusfloat*   fp
                                        )
{
    this->_field_point_j[0] = fp[0];
    this->_field_point_j[1] = fp[1];
    this->_field_point_j[2] = fp[2];
}


void    GRFDnInterface::set_source_i(
                                            SourceNode* source_node
                                    )
{
    this->_source_i = source_node;
}


void    GRFDnInterface::set_source_j(
                                            SourceNode* source_node
                                    )
{
    this->_source_j = source_node;
}