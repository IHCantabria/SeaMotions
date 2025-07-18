
// Include local modules
#include "gwfdn_interface_t.hpp"
#include "../math/shape_functions.hpp"
#include "../math/math_tools.hpp"
#include "../math/math_interface.hpp"


template<std::size_t N>
void    GWFDnInterfaceT<N>::_clear_heap(
                                            void
                                    )
{
    std::cerr << "Calling GWFDnInterfaceT destructor..." << std::endl;
    delete this->_integrals_db;
    delete this->_wave_data;
}


template<std::size_t N>
void    GWFDnInterfaceT<N>::_initialize(
                                            cusfloat    ang_freq
                                    )
{
    // Calculate wave numbers
    this->_wave_data    = new WaveDispersionFO( 
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


template<std::size_t N>
GWFDnInterfaceT<N>::GWFDnInterfaceT( 
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


template<std::size_t N>
GWFDnInterfaceT<N>::~GWFDnInterfaceT( void )
{
    this->_clear_heap( );
}


template<std::size_t N>
void        GWFDnInterfaceT<N>::operator()( 
                                            cusfloat*   xi,
                                            cusfloat*   eta,
                                            cusfloat*   x,
                                            cusfloat*   y,
                                            cusfloat*   z,
                                            cuscomplex* G,
                                            cuscomplex* dG_dn,
                                            cuscomplex* dG_dx,
                                            cuscomplex* dG_dy,
                                            cuscomplex* dG_dz
                                        )
{
    // Allocate space for local variables
    cuscomplex  G[N];
    cusfloat    R[N];
    cuscomplex  dG_dR[N];
    cusfloat    dX[N];
    cusfloat    dY[N];

    // Calculate horizontal radius
    cusfloat source_x = this->_source_j->position[0];
    cusfloat source_y = this->_source_j->position[1];
    for ( std::size_t i=0; i<N; i++ )
    {
        dX[i] = source_x - x[i];
        dY[i] = source_y - y[i];
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        R[i] = ( pow2s( dX[i] ) + pow2s( dY[i] ) );
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        R[i] = std::max( R[i], 1e-12 );
    }

    lv_sqrt( N, R, R );

    // Calculate Green function values
    wave_term_fin_depth_integral( 
                                    R, 
                                    this->_source_j->position[2],
                                    z,
                                    this->_water_depth,
                                    this->_bessel_factory, 
                                    *(this->_wave_data),
                                    G,
                                    dG_dR,
                                    dG_dz
                                );
    
    // Calculate X and Y cartesian coordinates derivatives
    for ( std::size_t i=0; i<N; i++ )
    {
        dG_dx[i] = dG_dR[i] * dX[i] / R[i];
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        dG_dy[i] = dG_dR[i] * dY[i] / R[i];
    }

    // Calculate normal derivate
    cusfloat nx = this->_source_j->normal_vec[0];
    cusfloat ny = this->_source_j->normal_vec[1];
    cusfloat nz = this->_source_j->normal_vec[2];

    for ( std::size_t i=0; i<N; i++ )
    {
        dG_dn[i] = (
                        dG_dx[i] * nx
                        +
                        dG_dy[i] * ny
                        +
                        dG_dz[i] * nz
                    );
    }
    // cuscomplex  dG_dn   =   (
    //                             dG_dX * this->_source_j->normal_vec[0]
    //                             +
    //                             dG_dY * this->_source_j->normal_vec[1]
    //                             +
    //                             dG_dZ * this->_source_j->normal_vec[2]
    //                         );
    
    // Get local shape function value
    // cusfloat    shp_val = shape_functions( 
    //                                         this->_source_i->p_order,
    //                                         this->_source_i->q_order,
    //                                         xi, 
    //                                         eta 
    //                                     );
}


template<std::size_t N>
void    GWFDnInterfaceT<N>::set_ang_freq(
                                            cusfloat    ang_freq
                                    )
{
    // Delete previous wave data instance
    delete this->_wave_data;

    // Create new wave data for the new angular frequency
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


template<std::size_t N>
void    GWFDnInterfaceT<N>::set_field_point(
                                            cusfloat*   fp
                                        )
{
    this->_field_point_j[0] = fp[0];
    this->_field_point_j[1] = fp[1];
    this->_field_point_j[2] = fp[2];
}


template<std::size_t N>
void    GWFDnInterfaceT<N>::set_source_i(
                                            SourceNode* source_node,
                                            cuscomplex  source_value
                                    )
{
    this->_source_i     = source_node;
    this->_source_value = source_value;
}


template<std::size_t N>
void    GWFDnInterfaceT<N>::set_source_j(
                                            SourceNode* source_node
                                    )
{
    this->_source_j = source_node;
}