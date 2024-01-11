
// Include local modules
#include "../green/pulsating_fin_depth.hpp"
#include "gwf_interface.hpp"
#include "../math/math_tools.hpp"
#include "../math/shape_functions.hpp"


GWFInterface::GWFInterface(
                            SourceNode*         source_in,
                            cuscomplex          source_value_in,
                            cusfloat*           field_point_in,
                            cusfloat            ang_freq,
                            cusfloat            water_depth_in,
                            cusfloat            grav_acc_in
                        )
{
    // Storage necessary input class arguements into class attributes
    this->_grav_acc     = grav_acc_in;
    this->_source       = source_in;
    this->_source_value = source_value_in;
    this->_water_depth  = water_depth_in;

    copy_vector( 3, field_point_in, this->_field_point );
    
    // Define wave dispersion data
    this->_wave_data    = new WaveDispersionFO(
                                                    ang_freq,
                                                    30,
                                                    water_depth_in,
                                                    grav_acc_in
                                                );
    this->_wave_data->calculate_john_terms( );

    // Define Integrals database
    this->_integrals_db = new IntegralsDb( );

    // Fold for current frequency and water depth
    cusfloat H = pow2s( ang_freq ) * water_depth_in / grav_acc_in;
    this->_integrals_db->fold_h( H );
}


GWFInterface::~GWFInterface(
                                void
                            )
{
    delete this->_wave_data;
    delete this->_integrals_db;
}


cuscomplex  GWFInterface::operator()(
                                        cusfloat    xi,
                                        cusfloat    eta,
                                        cusfloat    x,
                                        cusfloat    y,
                                        cusfloat    z
                                    )
{
    // Calculate Green function value
    cuscomplex  green_val   =   G_integral_wave(
                                                    this->_field_point[0],
                                                    this->_field_point[1],
                                                    this->_field_point[2],
                                                    x,
                                                    y,
                                                    z,
                                                    this->_water_depth,
                                                    *(this->_wave_data),
                                                    *(this->_integrals_db)
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


void GWFInterface::set_ang_freq(
                                    cusfloat ang_freq
                                )
{
    // Delete previous instance of WaveDispersionFO class
    delete this->_wave_data;

    // Create new WaveDispersionFO class
    this->_wave_data = new WaveDispersionFO(
                                                    ang_freq,
                                                    30,
                                                    this->_water_depth,
                                                    this->_grav_acc
                                                );
    this->_wave_data->calculate_john_terms( );

    // Fold for current frequency and water depth
    cusfloat H = pow2s( ang_freq ) * this->_water_depth / this->_grav_acc;
    this->_integrals_db->fold_h( H );
}


void GWFInterface::set_field_point(
                                        cusfloat* field_point
                                    )
{
    this->_field_point[0] = field_point[0];
    this->_field_point[1] = field_point[1];
    this->_field_point[2] = field_point[2];
}


void GWFInterface::set_source_i(
                                    SourceNode* source,
                                    cuscomplex  source_val
                                )
{
    this->_source       = source;
    this->_source_value = source_val;
}