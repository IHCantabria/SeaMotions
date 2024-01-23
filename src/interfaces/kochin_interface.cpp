
// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "kochin_interface.hpp"

#include "../math/shape_functions.hpp"
#include "../math/special_math.hpp"
#include "../waves/waves_common.hpp"


KochinInterface::KochinInterface(
                                    SourceNode* source,
                                    cusfloat    wave_number,
                                    cusfloat    water_depth,
                                    int         l_order,
                                    bool        is_cos
                                )
{
    // Storage input attributes
    this->_is_cos       = is_cos;
    this->_l_order      = l_order;
    this->_source       = source;
    this->_wave_num     = wave_number;
    this->_water_depth  = water_depth;
}


void        KochinInterface::set_is_cos(
                                            bool is_cos
                                        )
{
    this->_is_cos = is_cos;
}


cuscomplex  KochinInterface::operator()( 
                                            cusfloat xi,
                                            cusfloat eta,
                                            cusfloat x,
                                            cusfloat y,
                                            cusfloat z
                                        )
{
    // Calculate polar coordinates
    cusfloat    r_horz  = std::sqrt( 
                                        pow2s( x )
                                        +
                                        pow2s( y )
                                    );
    
    cusfloat    alpha   = std::atan2( y, x );

    // Calculate kochin integral function
    cuscomplex  kf  = (
                            wave_vertical_profile_fo( this->_wave_num, this->_water_depth, z )
                            *
                            this->_ep_l
                            *
                            this->_cf_l
                            *
                            std::cyl_bessel_j( 
                                                static_cast<cusfloat>( this->_l_order ),
                                                this->_wave_num * r_horz
                                            )
                        );

    if ( this->_is_cos )
    {
        kf *= cos( this->_l_order * alpha );
    }
    else
    {
        kf *= sin( this->_l_order * alpha );
    }


    kf *= - 1.0 / 4.0 / PI;
    

    // Get local shape function value
    cusfloat    shp_val = shape_functions( 
                                            this->_source->p_order,
                                            this->_source->q_order,
                                            xi, 
                                            eta 
                                        );

    return shp_val * kf;
}


void        KochinInterface::set_l_order( 
                                            int l_order 
                                        )
{
    this->_l_order  = l_order;
    this->_cf_l     = std::pow( cuscomplex( 0.0, -1.0 ), this->_l_order );
    this->_ep_l     = ep_n( l_order );
}


void        KochinInterface::set_source(
                                            SourceNode* source
                                        )
{
    this->_source   = source;
}