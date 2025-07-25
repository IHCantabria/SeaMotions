
// Include local modules
#include "gwfcns_interface_t.hpp"
#include "../math/shape_functions.hpp"
#include "../math/math_tools.hpp"
#include "../math/math_interface.hpp"


template<std::size_t N>
void    GWFcnsInterfaceT<N>::_initialize(
                                            cusfloat    ang_freq
                                    )
{
    // Calculate wave numbers
    this->_wave_data.update_full( 
                                    ang_freq,
                                    this->_water_depth,
                                    this->_grav_acc
                                );

}


template<std::size_t N>
GWFcnsInterfaceT<N>::GWFcnsInterfaceT( 
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
template<auto Kernel>
void        GWFcnsInterfaceT<N>::operator()( 
                                            cusfloat*   xi,
                                            cusfloat*   eta,
                                            cusfloat*   x,
                                            cusfloat*   y,
                                            cusfloat*   z
                                        )
{
    // Calculate horizontal radius
    cusfloat source_x = this->_source_j->position[0];
    cusfloat source_y = this->_source_j->position[1];
    for ( std::size_t i=0; i<N; i++ )
    {
        this->_dX[i] = source_x - x[i];
        this->_dY[i] = source_y - y[i];
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_R[i] = ( pow2s( this->_dX[i] ) + pow2s( this->_dY[i] ) );
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_R[i] = std::max( this->_R[i], 1e-12 );
    }

    lv_sqrt( N, this->_R, this->_R );

    // Calculate Green function values
    Kernel(
                this->_R, 
                this->_source_j->position[2],
                z,
                this->_water_depth,
                this->_bessel_factory, 
                this->_wave_data,
                this->G,
                this->_dG_dR,
                this->dG_dz
            );
    
    // Calculate X and Y cartesian coordinates derivatives
    for ( std::size_t i=0; i<N; i++ )
    {
        this->dG_dx[i] = this->_dG_dR[i] * this->_dX[i] / this->_R[i];
    }

    for ( std::size_t i=0; i<N; i++ )
    {
        this->_dG_dy[i] = this->_dG_dR[i] * this->_dY[i] / this->_R[i];
    }

    // Calculate normal derivate
    cusfloat nx = this->_source_j->normal_vec[0];
    cusfloat ny = this->_source_j->normal_vec[1];
    cusfloat nz = this->_source_j->normal_vec[2];

    for ( std::size_t i=0; i<N; i++ )
    {
        this->dG_dn[i] = (
                            this->dG_dx[i] * nx
                            +
                            this->dG_dy[i] * ny
                            +
                            this->dG_dz[i] * nz
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
void    GWFcnsInterfaceT<N>::set_ang_freq(
                                            cusfloat    ang_freq
                                    )
{
    // Create new wave data for the new angular frequency
    this->_wave_data.update_full( 
                                    ang_freq,
                                    this->_water_depth,
                                    this->_grav_acc 
                                );
}


template<std::size_t N>
void    GWFcnsInterfaceT<N>::set_field_point(
                                            cusfloat*   fp
                                        )
{
    this->_field_point_j[0] = fp[0];
    this->_field_point_j[1] = fp[1];
    this->_field_point_j[2] = fp[2];
}


template<std::size_t N>
void    GWFcnsInterfaceT<N>::set_source_i(
                                            SourceNode* source_node,
                                            cuscomplex  source_value
                                    )
{
    this->_source_i     = source_node;
    this->_source_value = source_value;
}


template<std::size_t N>
void    GWFcnsInterfaceT<N>::set_source_j(
                                            SourceNode* source_node
                                        )
{
    this->_source_j = source_node;
}