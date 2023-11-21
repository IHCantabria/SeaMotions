
// Include local modules
#include "gwfdn_interface.hpp"
#include "../math/shape_functions.hpp"


void    GWFDnInterface::_clear_heap(
                                            void
                                    )
{
    std::cerr << "Calling GWFDnInterface destructor..." << std::endl;
    delete this->_integrals_db;
    delete this->_wave_data;
}


void    GWFDnInterface::_initialize(
                                            cusfloat    ang_freq
                                    )
{
    // Calculate wave numbers
    this->_wave_data    = new WaveDispersionData( 
                                                    ang_freq,
                                                    30,
                                                    this->_water_depth,
                                                    this->_grav_acc
                                                );
    this->_wave_data->calculate_john_terms( );

    // Load integrals database
    this->_integrals_db = new IntegralsDb( );
    
    // Fold for current frequency and water depth
    cusfloat H = pow2s( ang_freq ) * this->_water_depth / this->_grav_acc;
    this->_integrals_db->fold_h( H );
}


GWFDnInterface::GWFDnInterface( 
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


GWFDnInterface::~GWFDnInterface( void )
{
    this->_clear_heap( );
}


cuscomplex  GWFDnInterface::operator()( 
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
    cuscomplex  dG_dR   = G_integral_wave_dr(
                                                R,
                                                this->_source_j->position[2],
                                                z,
                                                this->_water_depth,
                                                *(this->_wave_data),
                                                *(this->_integrals_db)
                                            );
    cuscomplex  dG_dZ   = G_integral_wave_dz(
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


void    GWFDnInterface::set_ang_freq(
                                            cusfloat    ang_freq
                                    )
{
    // Delete previous wave data instance
    delete this->_wave_data;

    // Create new wave data for the new angular frequency
    this->_wave_data = new WaveDispersionData( 
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


void    GWFDnInterface::set_field_point_j(
                                            cusfloat*   fp
                                        )
{
    this->_field_point_j[0] = fp[0];
    this->_field_point_j[1] = fp[1];
    this->_field_point_j[2] = fp[2];
}


void    GWFDnInterface::set_source_i(
                                            SourceNode* source_node
                                    )
{
    this->_source_i = source_node;
}


void    GWFDnInterface::set_source_j(
                                            SourceNode* source_node
                                    )
{
    this->_source_j = source_node;
}