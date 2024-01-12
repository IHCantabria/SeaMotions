
#ifndef __wave_dispersion_so_hpp
#define __wave_dispersion_so_hpp

// Include local modules
#include "../config.hpp"


struct WaveDispersionSO
{
private:
    // Define private class methods
    void    _update_status(
                                cusfloat    w0_in,
                                cusfloat    w1_in,
                                cusfloat    head_0_in,
                                cusfloat    head_1_in
                            );

public:
    // Define public class attributes
    cusfloat    a0              = 0.0;
    cusfloat    a1              = 0.0;
    cusfloat    grav_acc        = 0.0;
    cusfloat    head_0          = 0.0;
    cusfloat    head_1          = 0.0;
    cusfloat    k0              = 0.0;
    cusfloat    k0_deep_water   = 0.0;
    cusfloat    k0_vec[2]       = { 0.0, 0.0 };
    cusfloat    k1              = 0.0;
    cusfloat    k1_deep_water   = 0.0;
    cusfloat    k1_vec[2]       = { 0.0, 0.0 };
    cusfloat    k_diff_mod      = 0.0;
    cusfloat    k_diff_vec[2]   = { 0.0, 0.0 };
    cusfloat    k_sum_mod       = 0.0;
    cusfloat    k_sum_vec[2]    = { 0.0, 0.0 };
    cusfloat    w0              = 0.0;
    cusfloat    w1              = 0.0;
    cusfloat    water_depth     = 0.0;
    cusfloat    w_diff          = 0.0;
    cusfloat    w_k_diff        = 0.0;
    cusfloat    w_k_sum         = 0.0;
    cusfloat    w_sum           = 0.0;

    // Define class constructors and destructor
    WaveDispersionSO( ) = default;

    WaveDispersionSO(
                        cusfloat    a0_in,
                        cusfloat    a1_in,
                        cusfloat    w0_in,
                        cusfloat    w1_in,
                        cusfloat    head_0_in,
                        cusfloat    head_1_in,
                        cusfloat    water_depth,
                        cusfloat    grav_acc
                    );

    // Define class methods
    cusfloat    get_w_ds(
                                bool    is_diff
                        );

    void        set_new_data(
                                cusfloat w0_in,
                                cusfloat w1_in,
                                cusfloat head_0_in,
                                cusfloat head_1_in
                        );

};

#endif