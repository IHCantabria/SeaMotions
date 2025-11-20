
#ifndef __simulation_data_hpp
#define __simulation_data_hpp

// Include local modules
#include "../config.hpp"
#include "../inout/input.hpp"
#include "mpi_config.hpp"


struct SimulationData
{
private:
    // Declare private class attributes
    Input*          _input                      = nullptr;
    bool            _is_mdrift                  = false;
    bool            _is_qtf_base_freq           = false;
    bool            _is_qtf_data                = false;
    bool            _is_qtf_direct_freq         = false;
    bool            _is_qtf_indirect_freq       = false;
    MpiConfig*      _mpi_config                 = nullptr;

public:
    // Declare class attributes
    cusfloat*       added_mass                          = nullptr;
    cusfloat*       added_mass_p0                       = nullptr;
    cusfloat*       damping_rad                         = nullptr;
    cusfloat*       damping_rad_p0                      = nullptr;
    int             dofs_np                             = 0;
    cuscomplex*     froude_krylov                       = nullptr;
    cuscomplex*     froude_krylov_p0                    = nullptr;
    int             heads_np                            = 0;
    int             hydmech_np                          = 0;
    cusfloat*       hydrostiff_p0                       = nullptr;
    cuscomplex*     intensities                         = nullptr;
    cuscomplex*     intensities_p0                      = nullptr;
    
    cuscomplex*     mdrift                              = nullptr;
    cuscomplex*     mdrift_acc                          = nullptr;
    cuscomplex*     mdrift_bern                         = nullptr;
    cuscomplex*     mdrift_body_pot_fk                  = nullptr;
    cuscomplex*     mdrift_body_pot_raddif              = nullptr;
    cuscomplex*     mdrift_body_pot_total               = nullptr;
    cuscomplex*     mdrift_body_vel_x_fk                = nullptr;
    cuscomplex*     mdrift_body_vel_y_fk                = nullptr;
    cuscomplex*     mdrift_body_vel_z_fk                = nullptr;
    cuscomplex*     mdrift_body_vel_x_raddif            = nullptr;
    cuscomplex*     mdrift_body_vel_y_raddif            = nullptr;
    cuscomplex*     mdrift_body_vel_z_raddif            = nullptr;
    cuscomplex*     mdrift_body_vel_x_total             = nullptr;
    cuscomplex*     mdrift_body_vel_y_total             = nullptr;
    cuscomplex*     mdrift_body_vel_z_total             = nullptr;
    cuscomplex*     mdrift_fs_pot_fk                    = nullptr;
    cuscomplex*     mdrift_fs_pot_raddif                = nullptr;
    cuscomplex*     mdrift_fs_pot_total                 = nullptr;
    cuscomplex*     mdrift_fs_vel_x_fk                  = nullptr;
    cuscomplex*     mdrift_fs_vel_y_fk                  = nullptr;
    cuscomplex*     mdrift_fs_vel_z_fk                  = nullptr;
    cuscomplex*     mdrift_fs_vel_x_raddif              = nullptr;
    cuscomplex*     mdrift_fs_vel_y_raddif              = nullptr;
    cuscomplex*     mdrift_fs_vel_z_raddif              = nullptr;
    cuscomplex*     mdrift_fs_vel_x_total               = nullptr;
    cuscomplex*     mdrift_fs_vel_y_total               = nullptr;
    cuscomplex*     mdrift_fs_vel_z_total               = nullptr;
    cuscomplex*     mdrift_kochin_pert_cos              = nullptr;
    cuscomplex*     mdrift_kochin_pert_sin              = nullptr;
    cuscomplex*     mdrift_kochin_rad_cos               = nullptr;
    cuscomplex*     mdrift_kochin_rad_sin               = nullptr;
    cuscomplex*     mdrift_mom                          = nullptr;
    cuscomplex*     mdrift_pc_pot_fk                    = nullptr;
    cuscomplex*     mdrift_pc_pot_raddif                = nullptr;
    cuscomplex*     mdrift_pc_pot_total                 = nullptr;
    cuscomplex*     mdrift_pc_vel_x_fk                  = nullptr;
    cuscomplex*     mdrift_pc_vel_y_fk                  = nullptr;
    cuscomplex*     mdrift_pc_vel_z_fk                  = nullptr;
    cuscomplex*     mdrift_pc_vel_x_raddif              = nullptr;
    cuscomplex*     mdrift_pc_vel_y_raddif              = nullptr;
    cuscomplex*     mdrift_pc_vel_z_raddif              = nullptr;
    cuscomplex*     mdrift_pc_vel_x_total               = nullptr;
    cuscomplex*     mdrift_pc_vel_y_total               = nullptr;
    cuscomplex*     mdrift_pc_vel_z_total               = nullptr;
    cuscomplex*     mdrift_wl                           = nullptr;
    cuscomplex*     mdrift_wl_pot_fk                    = nullptr;
    cuscomplex*     mdrift_wl_pot_raddif                = nullptr;
    cuscomplex*     mdrift_wl_pot_total                 = nullptr;
    cuscomplex*     mdrift_wl_rel_we                    = nullptr;
    cuscomplex*     mdrift_wl_vel_x_fk                  = nullptr;
    cuscomplex*     mdrift_wl_vel_y_fk                  = nullptr;
    cuscomplex*     mdrift_wl_vel_z_fk                  = nullptr;
    cuscomplex*     mdrift_wl_vel_x_raddif              = nullptr;
    cuscomplex*     mdrift_wl_vel_y_raddif              = nullptr;
    cuscomplex*     mdrift_wl_vel_z_raddif              = nullptr;
    cuscomplex*     mdrift_wl_vel_x_total               = nullptr;
    cuscomplex*     mdrift_wl_vel_y_total               = nullptr;
    cuscomplex*     mdrift_wl_vel_z_total               = nullptr;
    cuscomplex*     mdrift_wl_we                        = nullptr;
    cuscomplex*     panels_potential                    = nullptr;
    cuscomplex*     panels_potential_p0                 = nullptr;
    cuscomplex*     panels_pressure                     = nullptr;
    cuscomplex*     panels_pressure_p0                  = nullptr;
    cuscomplex*     potential_secord_force              = nullptr;
    cuscomplex*     qtf                                 = nullptr;
    int             qtf_body_heads_np                   = 0;
    cuscomplex*     qtf_body_pot_raddif_freq            = nullptr;
    int             qtf_body_raddif_np                  = 0;
    cuscomplex*     qtf_body_vel_x_fk_freq              = nullptr;
    cuscomplex*     qtf_body_vel_y_fk_freq              = nullptr;
    cuscomplex*     qtf_body_vel_z_fk_freq              = nullptr;
    cuscomplex*     qtf_body_vel_x_raddif_freq          = nullptr;
    cuscomplex*     qtf_body_vel_y_raddif_freq          = nullptr;
    cuscomplex*     qtf_body_vel_z_raddif_freq          = nullptr;
    cuscomplex*     qtf_body_vel_x_total_freq           = nullptr;
    cuscomplex*     qtf_body_vel_y_total_freq           = nullptr;
    cuscomplex*     qtf_body_vel_z_total_freq           = nullptr;
    cuscomplex*     qtf_diff_acc                        = nullptr;
    cuscomplex*     qtf_diff_acc_freqs                  = nullptr;
    cuscomplex*     qtf_diff_bern                       = nullptr;
    cuscomplex*     qtf_diff_bern_freqs                 = nullptr;
    cuscomplex*     qtf_diff_body_force_p0              = nullptr;
    cuscomplex*     qtf_diff_body_force_freqs_p0        = nullptr;
    cuscomplex*     qtf_diff_freqs                      = nullptr;
    cuscomplex*     qtf_diff_froude_krylov_fo_p0        = nullptr;
    cuscomplex*     qtf_diff_froude_krylov_fo_freqs_p0  = nullptr;
    cuscomplex*     qtf_diff_fs_far_field_p0            = nullptr;
    cuscomplex*     qtf_diff_fs_far_field_freqs_p0      = nullptr;
    cuscomplex*     qtf_diff_fs_near_field_p0           = nullptr;
    cuscomplex*     qtf_diff_fs_near_field_freqs_p0     = nullptr;
    cuscomplex*     qtf_diff_mom                        = nullptr;
    cuscomplex*     qtf_diff_mom_freqs                  = nullptr;
    cuscomplex*     qtf_diff_secord_force               = nullptr;
    cuscomplex*     qtf_diff_secord_force_freqs         = nullptr;
    cuscomplex*     qtf_diff_wl                         = nullptr;
    cuscomplex*     qtf_diff_wl_freqs                   = nullptr;
    int             qtf_fs_heads_np                     = 0;
    cuscomplex*     qtf_fs_pot_fk_freq                  = nullptr;
    cuscomplex*     qtf_fs_pot_raddif_freq              = nullptr;
    cuscomplex*     qtf_fs_pot_total_freq               = nullptr;
    int             qtf_fs_raddif_np                    = 0;
    cuscomplex*     qtf_fs_vel_x_fk_freq                = nullptr;
    cuscomplex*     qtf_fs_vel_y_fk_freq                = nullptr;
    cuscomplex*     qtf_fs_vel_z_fk_freq                = nullptr;
    cuscomplex*     qtf_fs_vel_x_total_freq             = nullptr;
    cuscomplex*     qtf_fs_vel_y_total_freq             = nullptr;
    cuscomplex*     qtf_fs_vel_z_total_freq             = nullptr;
    int             qtf_kochin_heads_np                 = 0;
    int             qtf_kochin_rad_np                   = 0;
    cuscomplex*     qtf_kochin_pert_cos_freqs           = nullptr;
    cuscomplex*     qtf_kochin_pert_sin_freqs           = nullptr;
    cuscomplex*     qtf_kochin_rad_cos_freqs            = nullptr;
    cuscomplex*     qtf_kochin_rad_sin_freqs            = nullptr;
    int             qtf_pc_heads_np                     = 0;
    cuscomplex*     qtf_pc_pot_total_freq               = nullptr;
    int             qtf_pc_raddif_np                    = 0;
    cuscomplex*     qtf_pc_vel_x_total_freq             = nullptr;
    cuscomplex*     qtf_pc_vel_y_total_freq             = nullptr;
    cuscomplex*     qtf_pc_vel_z_total_freq             = nullptr;
    int             qtf_np                              = 0;
    cuscomplex*     qtf_raos_freq                       = nullptr;
    cuscomplex*     qtf_sum_acc                         = nullptr;
    cuscomplex*     qtf_sum_acc_freqs                   = nullptr;
    cuscomplex*     qtf_sum_bern                        = nullptr;
    cuscomplex*     qtf_sum_bern_freqs                  = nullptr;
    cuscomplex*     qtf_sum_body_force_p0               = nullptr;
    cuscomplex*     qtf_sum_body_force_freqs_p0         = nullptr;
    cuscomplex*     qtf_sum_freqs                       = nullptr;
    cuscomplex*     qtf_sum_froude_krylov_fo_p0         = nullptr;
    cuscomplex*     qtf_sum_froude_krylov_fo_freqs_p0   = nullptr;
    cuscomplex*     qtf_sum_fs_far_field_p0             = nullptr;
    cuscomplex*     qtf_sum_fs_far_field_freqs_p0       = nullptr;
    cuscomplex*     qtf_sum_fs_near_field_p0            = nullptr;
    cuscomplex*     qtf_sum_fs_near_field_freqs_p0      = nullptr;
    cuscomplex*     qtf_sum_mom                         = nullptr;
    cuscomplex*     qtf_sum_mom_freqs                   = nullptr;
    cuscomplex*     qtf_sum_secord_force                = nullptr;
    cuscomplex*     qtf_sum_secord_force_freqs          = nullptr;
    cuscomplex*     qtf_sum_wl                          = nullptr;
    cuscomplex*     qtf_sum_wl_freqs                    = nullptr;
    int             qtf_wl_heads_np                     = 0;
    cuscomplex*     qtf_wl_pot_raddif_freq              = nullptr;
    cuscomplex*     qtf_wl_pot_total_freq               = nullptr;
    int             qtf_wl_raddif_np                    = 0;
    cuscomplex*     qtf_wl_rel_we_total_freq            = nullptr;
    cuscomplex*     qtf_wl_vel_x_total_freq             = nullptr;
    cuscomplex*     qtf_wl_vel_y_total_freq             = nullptr;
    cuscomplex*     qtf_wl_vel_z_total_freq             = nullptr;
    cuscomplex*     raos                                = nullptr;
    cusfloat*       structural_mass_p0                  = nullptr;
    cuscomplex*     sysmat                              = nullptr;
    cuscomplex*     sysmat_steady                       = nullptr;
    cuscomplex*     wave_diffrac                        = nullptr;
    cuscomplex*     wave_diffrac_p0                     = nullptr;
    int             wave_exc_np                         = 0;
    cuscomplex*     wave_exc_p0                         = nullptr;

    // Declare class constructors and destructor
    SimulationData( ) = default;

    SimulationData(
                        Input*      input_in,
                        MpiConfig*  mpi_config_in,
                        int         bodies_np,
                        int         dofs_np,
                        int         heads_np,
                        int         rows_np
                    );

    ~SimulationData( 
                        void
                    );

    // Declare public class methods
    void    add_mean_drift_data(
                                        int body_panels_tnp,
                                        int wl_panels_tnp,
                                        int body_gp_np,
                                        int wl_gp_np
                                );

    void    add_qtf_base_data(
                                        int body_panels_tnp,
                                        int body_gp_np,
                                        int wl_panels_tnp,
                                        int wl_gp_np,
                                        int freqs_np
                            );

    void    add_qtf_data(
                                        int freqs_np
                        );

    void    add_qtf_direct_data(
                                        int body_panels_tnp,
                                        int body_gp_np,
                                        int fs_panels_tnp,
                                        int fs_gp_np,
                                        int wl_panels_tnp,
                                        int wl_gp_np,
                                        int pc_panels_tnp,
                                        int pc_gp_np,
                                        int freqs_np
                                );

    void    add_qtf_indirect_data(
                                        int body_panels_tnp,
                                        int body_gp_np,
                                        int fs_panels_tnp,
                                        int fs_gp_np,
                                        int wl_panels_tnp,
                                        int wl_gp_np,
                                        int freqs_np
                                );

    int     get_heads_np(
                                        int panels_tnp,
                                        int body_gp_np
                        );

    int     get_raddif_np(
                                        int panels_tnp,
                                        int body_gp_np
                        );

    void    storage_qtf_base_freq(
                                        int         freq_num,
                                        cuscomplex* qtf_body_vel_x_total,
                                        cuscomplex* qtf_body_vel_y_total,
                                        cuscomplex* qtf_body_vel_z_total,
                                        cuscomplex* raos,
                                        cuscomplex* qtf_wl_we_total
                                );

    void    storage_qtf_direct_freq(
                                        int         freq_num
                                    );

    void    storage_qtf_indirect_freq(
                                        int         freq_num,
                                        cuscomplex* qtf_body_pot_raddif,
                                        cuscomplex* qtf_body_vel_x_raddif,
                                        cuscomplex* qtf_body_vel_y_raddif,
                                        cuscomplex* qtf_body_vel_z_raddif,
                                        cuscomplex* qtf_fs_pot_fk,
                                        cuscomplex* qtf_fs_pot_raddif,
                                        cuscomplex* qtf_fs_pot_total,
                                        cuscomplex* qtf_fs_vel_x_raddif,
                                        cuscomplex* qtf_fs_vel_y_raddif,
                                        cuscomplex* qtf_fs_vel_z_raddif,
                                        cuscomplex* qtf_fs_vel_x_total,
                                        cuscomplex* qtf_fs_vel_y_total,
                                        cuscomplex* qtf_fs_vel_z_total,
                                        cuscomplex* qtf_wl_pot_raddif,
                                        cuscomplex* qtf_wl_vel_x_total,
                                        cuscomplex* qtf_wl_vel_y_total,
                                        cuscomplex* qtf_wl_vel_z_total,
                                        cuscomplex* qtf_kochin_pert_cos,
                                        cuscomplex* qtf_kochin_pert_sin,
                                        cuscomplex* qtf_kochin_rad_cos,
                                        cuscomplex* qtf_kochin_rad_sin
                                    );
    
};

#endif