
// Include general usage libraries
#include <cassert>

// Include local modules
#include "../math/math_tools.hpp"
#include "simulation_data.hpp"


void    SimulationData::add_mean_drift_data(
                                                int body_panels_tnp,
                                                int wl_panels_tnp,
                                                int body_gp_np,
                                                int wl_gp_np
                                            )
{
    int body_raddif_np  =   this->get_raddif_np( body_panels_tnp,  body_gp_np );
    int body_heads_np   =   this->get_heads_np( body_panels_tnp,  body_gp_np );
    int wl_raddif_np    =   this->get_raddif_np( wl_panels_tnp, wl_gp_np );
    int wl_heads_np     =   this->get_heads_np( wl_panels_tnp,  wl_gp_np );
    if ( this->_mpi_config->is_root( ) )
    {
        this->mdrift                        = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->mdrift_wl                     = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->mdrift_bern                   = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->mdrift_acc                    = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->mdrift_mom                    = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->mdrift_body_vel_x_fk          = generate_empty_vector<cuscomplex>( body_heads_np );
        this->mdrift_body_vel_y_fk          = generate_empty_vector<cuscomplex>( body_heads_np );
        this->mdrift_body_vel_z_fk          = generate_empty_vector<cuscomplex>( body_heads_np );
        this->mdrift_body_vel_x_raddif      = generate_empty_vector<cuscomplex>( body_raddif_np );
        this->mdrift_body_vel_y_raddif      = generate_empty_vector<cuscomplex>( body_raddif_np );
        this->mdrift_body_vel_z_raddif      = generate_empty_vector<cuscomplex>( body_raddif_np );
        this->mdrift_body_vel_x_total       = generate_empty_vector<cuscomplex>( body_heads_np );
        this->mdrift_body_vel_y_total       = generate_empty_vector<cuscomplex>( body_heads_np );
        this->mdrift_body_vel_z_total       = generate_empty_vector<cuscomplex>( body_heads_np );
        this->mdrift_wl_rel_we              = generate_empty_vector<cuscomplex>( wl_heads_np );
        this->mdrift_wl_we_fk               = generate_empty_vector<cuscomplex>( wl_heads_np );
        this->mdrift_wl_we_raddif           = generate_empty_vector<cuscomplex>( wl_raddif_np );
        this->mdrift_wl_we_total            = generate_empty_vector<cuscomplex>( wl_heads_np );
    }
    this->_is_mdrift = true;
}


void    SimulationData::add_qtf_base_data(
                                                int body_panels_tnp,
                                                int body_gp_np,
                                                int wl_panels_tnp,
                                                int wl_gp_np,
                                                int freqs_np
                                        )
{
    this->qtf_body_raddif_np    = this->get_raddif_np( body_panels_tnp,  body_gp_np );
    this->qtf_body_heads_np     = this->get_heads_np( body_panels_tnp,  body_gp_np );
    int body_heads_freq_np      = this->qtf_body_heads_np * freqs_np;

    this->qtf_wl_raddif_np      = this->get_raddif_np( wl_panels_tnp, wl_gp_np );
    this->qtf_wl_heads_np       = this->get_heads_np( wl_panels_tnp,  wl_gp_np );
    int wl_freqs_heads_np       = this->qtf_wl_heads_np * freqs_np;

    int wex_freq_np             = this->wave_exc_np * freqs_np;

    if ( this->_mpi_config->is_root( ) )
    {
        this->qtf_body_vel_x_total_freq     = generate_empty_vector<cuscomplex>( body_heads_freq_np );
        this->qtf_body_vel_y_total_freq     = generate_empty_vector<cuscomplex>( body_heads_freq_np );
        this->qtf_body_vel_z_total_freq     = generate_empty_vector<cuscomplex>( body_heads_freq_np );
        this->qtf_raos_freq                 = generate_empty_vector<cuscomplex>( wex_freq_np );
        this->qtf_wl_we_total_freq          = generate_empty_vector<cuscomplex>( wl_freqs_heads_np );
    }
    this->_is_qtf_base_freq = true;
}


void    SimulationData::add_qtf_data(
                                        int freqs_np
                                    )
{
    int qtf_freq_np = this->qtf_np * pow2s( freqs_np );
    if ( this->_mpi_config->is_root( ) )
    {
        // Define common variables to storage QTF values
        this->qtf                       = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_diff_acc              = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_diff_bern             = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_diff_mom              = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_diff_secord_force     = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_diff_wl               = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_sum_acc               = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_sum_bern              = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_sum_mom               = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_sum_secord_force      = generate_empty_vector<cuscomplex>( this->qtf_np );
        this->qtf_sum_wl                = generate_empty_vector<cuscomplex>( this->qtf_np );

        // Define variables used for Indirect method
        if ( this->_input->out_qtf_so_model == 1 )
        {
            this->qtf_diff_froude_krylov_fo_p0  = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_diff_body_force_p0        = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_diff_fs_near_field_p0     = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_diff_fs_far_field_p0      = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_sum_froude_krylov_fo_p0   = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_sum_body_force_p0         = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_sum_fs_near_field_p0      = generate_empty_vector<cuscomplex>( this->qtf_np );
            this->qtf_sum_fs_far_field_p0       = generate_empty_vector<cuscomplex>( this->qtf_np );
        }

        if ( this->_input->out_qtf_comp )
        {
            // Define common variables to storage QTF values
            this->qtf_diff_acc_freqs            = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_diff_bern_freqs           = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_diff_freqs                = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_diff_mom_freqs            = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_diff_secord_force_freqs   = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_diff_wl_freqs             = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_sum_acc_freqs             = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_sum_bern_freqs            = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_sum_freqs                 = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_sum_mom_freqs             = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_sum_secord_force_freqs    = generate_empty_vector<cuscomplex>( qtf_freq_np );
            this->qtf_sum_wl_freqs              = generate_empty_vector<cuscomplex>( qtf_freq_np );

            // Define variables used for Indirect method
            if ( this->_input->out_qtf_so_model == 1 )
            {
                this->qtf_diff_froude_krylov_fo_freqs_p0    = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_diff_body_force_freqs_p0          = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_diff_fs_near_field_freqs_p0       = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_diff_fs_far_field_freqs_p0        = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_sum_froude_krylov_fo_freqs_p0     = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_sum_body_force_freqs_p0           = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_sum_fs_near_field_freqs_p0        = generate_empty_vector<cuscomplex>( qtf_freq_np );
                this->qtf_sum_fs_far_field_freqs_p0         = generate_empty_vector<cuscomplex>( qtf_freq_np );
            }
        }
    }
    this->_is_qtf_data = true;
}


void    SimulationData::add_qtf_indirect_data(
                                                int body_panels_tnp,
                                                int body_gp_np,
                                                int wl_panels_tnp,
                                                int wl_gp_np,
                                                int freqs_np
                                            )
{
    assert( this->_is_qtf_base_freq && "It is required to load qtf base data prior to load qtf indirect data." );

    int body_raddif_freq_np     = this->qtf_body_raddif_np * freqs_np;
    int body_heads_freq_np      = this->qtf_body_heads_np * freqs_np;

    int wl_freqs_raddif_np      = this->qtf_wl_raddif_np * freqs_np;
    int wl_freqs_heads_np       = this->qtf_wl_heads_np * freqs_np;

    int wex_freq_np             = this->wave_exc_np * freqs_np;

    if ( this->_mpi_config->is_root( ) )
    {
        this->mdrift_wl_vel_x_fk            = generate_empty_vector<cuscomplex>( this->qtf_wl_heads_np );
        this->mdrift_wl_vel_y_fk            = generate_empty_vector<cuscomplex>( this->qtf_wl_heads_np );
        this->mdrift_wl_vel_z_fk            = generate_empty_vector<cuscomplex>( this->qtf_wl_heads_np );
        this->mdrift_wl_vel_x_raddif        = generate_empty_vector<cuscomplex>( this->qtf_wl_raddif_np );
        this->mdrift_wl_vel_y_raddif        = generate_empty_vector<cuscomplex>( this->qtf_wl_raddif_np );
        this->mdrift_wl_vel_z_raddif        = generate_empty_vector<cuscomplex>( this->qtf_wl_raddif_np );
        this->mdrift_wl_vel_x_total         = generate_empty_vector<cuscomplex>( this->qtf_wl_heads_np );
        this->mdrift_wl_vel_y_total         = generate_empty_vector<cuscomplex>( this->qtf_wl_heads_np );
        this->mdrift_wl_vel_z_total         = generate_empty_vector<cuscomplex>( this->qtf_wl_heads_np );
        this->qtf_body_pot_raddif_freq      = generate_empty_vector<cuscomplex>( body_raddif_freq_np );
        this->qtf_wl_pot_raddif_freq        = generate_empty_vector<cuscomplex>( wl_freqs_raddif_np );
        this->qtf_wl_vel_x_total_freq       = generate_empty_vector<cuscomplex>( wl_freqs_heads_np );
        this->qtf_wl_vel_y_total_freq       = generate_empty_vector<cuscomplex>( wl_freqs_heads_np );
        this->qtf_wl_vel_z_total_freq       = generate_empty_vector<cuscomplex>( wl_freqs_heads_np );
    }
    this->_is_qtf_indirect_freq = true;
}


int     SimulationData::get_heads_np(
                                                int panels_tnp,
                                                int body_gp_np
                                    )
{
    return this->heads_np * panels_tnp * body_gp_np;
}


int     SimulationData::get_raddif_np(
                                                int panels_tnp,
                                                int body_gp_np
                                    )
{
    return ( this->heads_np + this->dofs_np ) * panels_tnp * body_gp_np;
}


SimulationData::SimulationData(
                                    Input*      input_in,
                                    MpiConfig*  mpi_config_in,
                                    int         bodies_np_in,
                                    int         dofs_np_in,
                                    int         heads_np_in,
                                    int         rows_local_np,
                                    int         rows_np
                                )
{
    // Storage input arguments into class attributes
    this->dofs_np           = dofs_np_in;
    this->heads_np          = heads_np_in;
    this->hydmech_np        = pow2s( dofs_np_in * bodies_np_in );
    this->_input            = input_in;
    this->_mpi_config       = mpi_config_in;
    this->qtf_np            = pow2s( heads_np_in ) * bodies_np_in * dofs_np_in;
    this->wave_exc_np       = heads_np_in * bodies_np_in * dofs_np_in;

    // Allocate space variables used in all the processes
    this->added_mass        = generate_empty_vector<cusfloat>( this->hydmech_np );
    this->damping_rad       = generate_empty_vector<cusfloat>( this->hydmech_np );
    this->froude_krylov     = generate_empty_vector<cuscomplex>( this->wave_exc_np );
    this->raos              = generate_empty_vector<cuscomplex>( this->wave_exc_np );
    this->intensities       = generate_empty_vector<cuscomplex>( ( dofs_np_in + heads_np_in ) * rows_local_np );
    this->panels_potential  = generate_empty_vector<cuscomplex>( ( dofs_np_in + heads_np_in ) * rows_np );
    this->wave_diffrac      = generate_empty_vector<cuscomplex>( this->wave_exc_np );

    // Allocate space for variables used only on root processor
    if ( this->_mpi_config->is_root( ) )
    {
        this->added_mass_p0         = generate_empty_vector<cusfloat>( hydmech_np );
        this->damping_rad_p0        = generate_empty_vector<cusfloat>( hydmech_np );
        this->froude_krylov_p0      = generate_empty_vector<cuscomplex>( wave_exc_np );
        this->hydrostiff_p0         = generate_empty_vector<cusfloat>( hydmech_np );
        this->structural_mass_p0    = generate_empty_vector<cusfloat>( hydmech_np );
        this->wave_diffrac_p0       = generate_empty_vector<cuscomplex>( wave_exc_np );
        this->wave_exc_p0           = generate_empty_vector<cuscomplex>( wave_exc_np );
    }
}


SimulationData::~SimulationData(
                                    void
                                )
{
    mkl_free( this->added_mass );
    mkl_free( this->damping_rad );
    mkl_free( this->froude_krylov );
    mkl_free( this->raos );
    mkl_free( this->intensities );
    mkl_free( this->panels_potential );
    mkl_free( this->sysmat );
    mkl_free( this->sysmat_steady );
    mkl_free( this->wave_diffrac );

    if ( this->_mpi_config->is_root( ) )
    {
        mkl_free( this->added_mass_p0 );
        mkl_free( this->damping_rad_p0 );
        mkl_free( this->froude_krylov_p0 );
        mkl_free( this->hydrostiff_p0 );
        mkl_free( this->structural_mass_p0 );
        mkl_free( this->wave_diffrac_p0 );
        mkl_free( this->wave_exc_p0 );
    }

    if ( this->_is_mdrift )
    {
        if ( this->_mpi_config->is_root( ) )
        {
            mkl_free( this->mdrift );
            mkl_free( this->mdrift_wl );
            mkl_free( this->mdrift_bern );
            mkl_free( this->mdrift_acc );
            mkl_free( this->mdrift_mom );
            mkl_free( this->mdrift_body_vel_x_fk );
            mkl_free( this->mdrift_body_vel_y_fk );
            mkl_free( this->mdrift_body_vel_z_fk );
            mkl_free( this->mdrift_body_vel_x_raddif );
            mkl_free( this->mdrift_body_vel_y_raddif );
            mkl_free( this->mdrift_body_vel_z_raddif );
            mkl_free( this->mdrift_body_vel_x_total );
            mkl_free( this->mdrift_body_vel_y_total );
            mkl_free( this->mdrift_body_vel_z_total );
            mkl_free( this->mdrift_wl_rel_we );
            mkl_free( this->mdrift_wl_we_fk );
            mkl_free( this->mdrift_wl_we_raddif );
            mkl_free( this->mdrift_wl_we_total );
        }
    }

    if ( this->_is_qtf_base_freq )
    {
        if ( this->_mpi_config->is_root( ) )
        {
            mkl_free( this->qtf_body_vel_x_total_freq );
            mkl_free( this->qtf_body_vel_y_total_freq );
            mkl_free( this->qtf_body_vel_z_total_freq );
            mkl_free( this->qtf_raos_freq             );
            mkl_free( this->qtf_wl_we_total_freq      );
            mkl_free( this->qtf_body_vel_x_total_freq );
            mkl_free( this->qtf_body_vel_y_total_freq );
            mkl_free( this->qtf_body_vel_z_total_freq );
        }
    }

    if ( this->_is_qtf_data )
    {
        if ( this->_mpi_config->is_root( ) )
        {
            mkl_free( this->qtf                   );
            mkl_free( this->qtf_diff_acc          );
            mkl_free( this->qtf_diff_bern         );
            mkl_free( this->qtf_diff_mom          );
            mkl_free( this->qtf_diff_secord_force );
            mkl_free( this->qtf_diff_wl           );
            mkl_free( this->qtf_sum_acc           );
            mkl_free( this->qtf_sum_bern          );
            mkl_free( this->qtf_sum_mom           );
            mkl_free( this->qtf_sum_secord_force  );
            mkl_free( this->qtf_sum_wl            );

            if ( this->_input->out_qtf_so_model == 1 )
            {
                mkl_free( this->qtf_diff_froude_krylov_fo_p0 );
                mkl_free( this->qtf_diff_body_force_p0       );
                mkl_free( this->qtf_diff_fs_near_field_p0    );
                mkl_free( this->qtf_diff_fs_far_field_p0     );
                mkl_free( this->qtf_sum_froude_krylov_fo_p0  );
                mkl_free( this->qtf_sum_body_force_p0        );
                mkl_free( this->qtf_sum_fs_near_field_p0     );
                mkl_free( this->qtf_sum_fs_far_field_p0      );
            }

            if ( this->_input->out_qtf_comp )
            {
                mkl_free( this->qtf_diff_acc_freqs          );
                mkl_free( this->qtf_diff_bern_freqs         );
                mkl_free( this->qtf_diff_freqs              );
                mkl_free( this->qtf_diff_mom_freqs          );
                mkl_free( this->qtf_diff_secord_force_freqs );
                mkl_free( this->qtf_diff_wl_freqs           );
                mkl_free( this->qtf_sum_acc_freqs           );
                mkl_free( this->qtf_sum_bern_freqs          );
                mkl_free( this->qtf_sum_freqs               );
                mkl_free( this->qtf_sum_mom_freqs           );
                mkl_free( this->qtf_sum_secord_force_freqs  );
                mkl_free( this->qtf_sum_wl_freqs            );

                if ( this->_input->out_qtf_so_model == 1 )
                {
                    mkl_free( this->qtf_diff_froude_krylov_fo_freqs_p0 );
                    mkl_free( this->qtf_diff_body_force_freqs_p0       );
                    mkl_free( this->qtf_diff_fs_near_field_freqs_p0    );
                    mkl_free( this->qtf_diff_fs_far_field_freqs_p0     );
                    mkl_free( this->qtf_sum_froude_krylov_fo_freqs_p0  );
                    mkl_free( this->qtf_sum_body_force_freqs_p0        );
                    mkl_free( this->qtf_sum_fs_near_field_freqs_p0     );
                    mkl_free( this->qtf_sum_fs_far_field_freqs_p0      );
                }
            }
        }
    }

    if ( this->_is_qtf_indirect_freq )
    {
        mkl_free( this->mdrift_wl_vel_x_fk          );
        mkl_free( this->mdrift_wl_vel_y_fk          );
        mkl_free( this->mdrift_wl_vel_z_fk          );
        mkl_free( this->mdrift_wl_vel_x_raddif      );
        mkl_free( this->mdrift_wl_vel_y_raddif      );
        mkl_free( this->mdrift_wl_vel_z_raddif      );
        mkl_free( this->mdrift_wl_vel_x_total       );
        mkl_free( this->mdrift_wl_vel_y_total       );
        mkl_free( this->mdrift_wl_vel_z_total       );
        mkl_free( this->qtf_body_pot_raddif_freq    );
        mkl_free( this->qtf_wl_pot_raddif_freq      );
        mkl_free( this->qtf_wl_vel_x_total_freq     );
        mkl_free( this->qtf_wl_vel_y_total_freq     );
        mkl_free( this->qtf_wl_vel_z_total_freq     );
    }
}


void    SimulationData::storage_qtf_base_freq(
                                                int         freq_num,
                                                cuscomplex* qtf_body_vel_x_total,
                                                cuscomplex* qtf_body_vel_y_total,
                                                cuscomplex* qtf_body_vel_z_total,
                                                cuscomplex* raos,
                                                cuscomplex* qtf_wl_we_total
                                            )
{
    if ( this->_mpi_config->is_root( ) )
    {
        // Storage body data
        int idx0    = freq_num * this->qtf_body_heads_np;
        int idx1    = 0;

        for ( int i=0; i<this->qtf_body_heads_np; i++ )
        {
            idx1 = idx0 + i;
            this->qtf_body_vel_x_total_freq[idx1]   = qtf_body_vel_x_total[i];
            this->qtf_body_vel_y_total_freq[idx1]   = qtf_body_vel_y_total[i];
            this->qtf_body_vel_z_total_freq[idx1]   = qtf_body_vel_z_total[i];
        }

        // Storage raos data
        idx0    = freq_num * this->wave_exc_np;
        idx1    = 0;

        for ( int i=0; i<this->wave_exc_np; i++ )
        {
            idx1                        = idx0 + i;
            this->qtf_raos_freq[idx1]   = raos[i];
        }

        // Storage WL data
        idx0    = freq_num * this->qtf_wl_heads_np;
        idx1    = 0;

        for ( int i=0; i<this->qtf_wl_heads_np; i++ )
        {
            idx1                                = idx0 + i;
            this->qtf_wl_we_total_freq[idx1]    = qtf_wl_we_total[i];
        }
    }
}


void    SimulationData::storage_qtf_indirect_freq(
                                                    int         freq_num,
                                                    cuscomplex* qtf_body_pot_raddif,
                                                    cuscomplex* qtf_wl_pot_raddif,
                                                    cuscomplex* qtf_wl_vel_x_total,
                                                    cuscomplex* qtf_wl_vel_y_total,
                                                    cuscomplex* qtf_wl_vel_z_total
                                                )
{
    if ( this->_mpi_config->is_root( ) )
    {
        // Storage body velocity raddiation data
        int idx0    = freq_num * this->qtf_body_raddif_np;
        int idx1    = 0;

        for ( int i=0; i<this->qtf_body_raddif_np; i++ )
        {
            idx1 = idx0 + i;
            this->qtf_body_pot_raddif_freq[idx1] = qtf_body_pot_raddif[i];
        }

        // Storage WL velocity raddiation data
        idx0    = freq_num * this->qtf_wl_raddif_np;
        idx1    = 0;

        for ( int i=0; i<this->qtf_wl_raddif_np; i++ )
        {
            idx1 = idx0 + i;
            this->qtf_wl_pot_raddif_freq[idx1] = qtf_wl_pot_raddif[i];
        }

        // Storage WL velocity total data
        idx0    = freq_num * this->qtf_wl_heads_np;
        idx1    = 0;

        for ( int i=0; i<this->qtf_wl_heads_np; i++ )
        {
            idx1 = idx0 + i;
            this->qtf_wl_vel_x_total_freq[idx1]  = qtf_wl_vel_x_total[i];
            this->qtf_wl_vel_y_total_freq[idx1]  = qtf_wl_vel_y_total[i];
            this->qtf_wl_vel_z_total_freq[idx1]  = qtf_wl_vel_z_total[i];
        }
    }
}