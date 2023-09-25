
// Include local modules
#include "../green/pulsating_fin_depth.hpp"
#include "grf_interface.hpp"
#include "../math/math_tools.hpp"
#include "../math/shape_functions.hpp"


GRFInterface::GRFInterface(
                            SourceNode*         source_in,
                            cuscomplex          source_value_in,
                            cusfloat*           field_point_in,
                            cusfloat            water_depth_in
                        )
{
    // Storage necessary input class arguements into class attributes
    this->_field_point  = field_point_in;
    this->_source       = source_in;
    this->_source_value = source_value_in;
    this->_water_depth  = water_depth_in;

}


GRFInterface::~GRFInterface(
                                void
                            )
{
}


cuscomplex  GRFInterface::operator()(
                                        cusfloat    xi,
                                        cusfloat    eta,
                                        cusfloat    x,
                                        cusfloat    y,
                                        cusfloat    z
                                    )
{
    // Calculate Green function value
    cuscomplex  green_val   =   G_integral_steady(
                                                    this->_field_point[0],
                                                    this->_field_point[1],
                                                    this->_field_point[2],
                                                    x,
                                                    y,
                                                    z,
                                                    this->_water_depth
                                                );
    
    // Get green function value
    cusfloat    shp_val     =   shape_functions(
                                                    this->_source->p_order,
                                                    this->_source->q_order,
                                                    xi,
                                                    eta
                                                );

    return this->_source_value * shp_val * green_val;
}


void GRFInterface::set_field_point(
                                        cusfloat* field_point
                                    )
{
    this->_field_point[0] = field_point[0];
    this->_field_point[1] = field_point[1];
    this->_field_point[2] = field_point[2];
}


void GRFInterface::set_source(
                                    SourceNode* source,
                                    cuscomplex  source_val
                                )
{
    this->_source       = source;
    this->_source_value = source_val;
}