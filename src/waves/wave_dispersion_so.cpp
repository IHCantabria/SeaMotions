
// Include local modules
#include "wave_dispersion_so.hpp"

#include "waves_common.hpp"


void    WaveDispersionSO::set_new_data(
                                            cusfloat w0_in,
                                            cusfloat w1_in,
                                            cusfloat head_0_in,
                                            cusfloat head_1_in
                                        )
{
    this->_update_status(
                            w0_in,
                            w1_in,
                            head_0_in,
                            head_1_in
                        );
}


void    WaveDispersionSO::_update_status(
                                            cusfloat    w0_in,
                                            cusfloat    w1_in,
                                            cusfloat    head_0_in,
                                            cusfloat    head_1_in
                                        )
{
    // Storage incoming attributes
    this->head_0        = head_0_in;
    this->head_1        = head_1_in;
    this->w0            = w0_in;
    this->w1            = w1_in;

    // Calculate wave numbers
    this->k0            = w2k( this->w0, this->water_depth, this->grav_acc );
    this->k1            = w2k( this->w1, this->water_depth, this->grav_acc );

    this->k0_deep_water = pow2s( this->w0 ) / this->grav_acc;
    this->k1_deep_water = pow2s( this->w1 ) / this->grav_acc;

    this->k0_vec[0]     = this->k0 * std::cos( this->head_0 );
    this->k0_vec[1]     = this->k0 * std::sin( this->head_0 );

    this->k1_vec[0]     = this->k1 * std::cos( this->head_1 );
    this->k1_vec[1]     = this->k1 * std::sin( this->head_1 );

    // Calculate sum and difference terms
    this->w_diff        = this->w0 - this->w1;
    this->w_sum         = this->w0 + this->w1;

    this->k_diff_mod    = std::sqrt(
                                        pow2s( this->k0 )
                                        +
                                        pow2s( this->k1 )
                                        -
                                        2 * this->k0 * this->k1 * std::cos( this->head_0 - this->head_1 )
                                    );
    
    this->k_sum_mod     = std::sqrt(
                                        pow2s( this->k0 )
                                        +
                                        pow2s( this->k1 )
                                        +
                                        2 * this->k0 * this->k1 * std::cos( this->head_0 - this->head_1 )
                                    );

    sv_sub( 
                2,
                this->k0_vec, 
                this->k1_vec,
                this->k_diff_vec
            );
    sv_add( 
                2, 
                this->k0_vec,
                this->k1_vec,
                this->k_sum_vec
            );

    // Calculate equivalent angular frequency for the wave lenghts difference
    this->w_k_diff  = k2w( this->k_diff_mod, this->water_depth, this->grav_acc );
    this->w_k_sum   = k2w( this->k_sum_mod, this->water_depth, this->grav_acc );

}


WaveDispersionSO::WaveDispersionSO(
                                        cusfloat    a0_in,
                                        cusfloat    a1_in,
                                        cusfloat    w0_in,
                                        cusfloat    w1_in,
                                        cusfloat    head_0_in,
                                        cusfloat    head_1_in,
                                        cusfloat    water_depth_in,
                                        cusfloat    grav_acc_in
                                    )
{
    // Storage general usage input arguments
    this->a0            = a0_in;
    this->a1            = a1_in;
    this->grav_acc      = grav_acc_in;
    this->water_depth   = water_depth_in;

    // Update object status
    this->_update_status( 
                            w0_in,
                            w1_in,
                            head_0_in,
                            head_1_in
                        );
}