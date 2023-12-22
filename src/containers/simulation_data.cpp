
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
    int body_raddif_np  =   ( this->heads_np + this->dofs_np ) * body_panels_tnp * body_gp_np;
    int body_heads_np   =   this->heads_np * body_panels_tnp * body_gp_np;
    int wl_raddif_np    =   ( this->heads_np + this->dofs_np ) * wl_panels_tnp * wl_gp_np;
    int wl_heads_np     =   this->heads_np * wl_panels_tnp * wl_gp_np;
    if ( this->_mpi_config->is_root( ) )
    {
        this->mdrift                    = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->qtf_wl                    = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->qtf_bern                  = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->qtf_acc                   = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->qtf_mom                   = generate_empty_vector<cuscomplex>( this->wave_exc_np );
        this->qtf_body_vel_x_fk         = generate_empty_vector<cuscomplex>( body_heads_np );
        this->qtf_body_vel_y_fk         = generate_empty_vector<cuscomplex>( body_heads_np );
        this->qtf_body_vel_z_fk         = generate_empty_vector<cuscomplex>( body_heads_np );
        this->qtf_body_vel_x_raddif     = generate_empty_vector<cuscomplex>( body_raddif_np );
        this->qtf_body_vel_y_raddif     = generate_empty_vector<cuscomplex>( body_raddif_np );
        this->qtf_body_vel_z_raddif     = generate_empty_vector<cuscomplex>( body_raddif_np );
        this->qtf_body_vel_x_total      = generate_empty_vector<cuscomplex>( body_heads_np );
        this->qtf_body_vel_y_total      = generate_empty_vector<cuscomplex>( body_heads_np );
        this->qtf_body_vel_z_total      = generate_empty_vector<cuscomplex>( body_heads_np );
        this->qtf_wl_rel_we             = generate_empty_vector<cuscomplex>( wl_heads_np );
        this->qtf_wl_we_fk              = generate_empty_vector<cuscomplex>( wl_heads_np );
        this->qtf_wl_we_raddif          = generate_empty_vector<cuscomplex>( wl_raddif_np );
        this->qtf_wl_we_total           = generate_empty_vector<cuscomplex>( wl_heads_np );
        this->potential_secord_force    = generate_empty_vector<cuscomplex>( this->wave_exc_np );
    }
    this->_is_mdrift = true;
}


SimulationData::SimulationData(
                                    int         bodies_np_in,
                                    int         dofs_np_in,
                                    int         heads_np_in,
                                    int         rows_local_np,
                                    int         rows_np,
                                    MpiConfig*  mpi_config_in
                                )
{
    // Storage input arguments into class attributes
    this->dofs_np           = dofs_np_in;
    this->heads_np          = heads_np_in;
    this->hydmech_np        = pow2s( dofs_np_in * bodies_np_in );
    this->_mpi_config       = mpi_config_in;  
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
            mkl_free( this->qtf_wl );
            mkl_free( this->qtf_bern );
            mkl_free( this->qtf_acc );
            mkl_free( this->qtf_mom );
            mkl_free( this->qtf_body_vel_x_fk );
            mkl_free( this->qtf_body_vel_y_fk );
            mkl_free( this->qtf_body_vel_z_fk );
            mkl_free( this->qtf_body_vel_x_raddif );
            mkl_free( this->qtf_body_vel_y_raddif );
            mkl_free( this->qtf_body_vel_z_raddif );
            mkl_free( this->qtf_body_vel_x_total );
            mkl_free( this->qtf_body_vel_y_total );
            mkl_free( this->qtf_body_vel_z_total );
            mkl_free( this->qtf_wl_rel_we );
            mkl_free( this->qtf_wl_we_fk );
            mkl_free( this->qtf_wl_we_raddif );
            mkl_free( this->qtf_wl_we_total );
            mkl_free( this->potential_secord_force );
        }
    }
}