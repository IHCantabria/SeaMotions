
#ifndef __simulation_data_hpp
#define __simulation_data_hpp

// Include local modules
#include "../config.hpp"
#include "mpi_config.hpp"


struct SimulationData
{
private:
    // Declare private class attributes
    bool            _is_mdrift              = false;
    MpiConfig*      _mpi_config             = nullptr;

public:
    // Declare class attributes
    cusfloat*       added_mass              = nullptr;
    cusfloat*       added_mass_p0           = nullptr;
    cusfloat*       damping_rad             = nullptr;
    cusfloat*       damping_rad_p0          = nullptr;
    int             dofs_np                 = 0;
    cuscomplex*     froude_krylov           = nullptr;
    cuscomplex*     froude_krylov_p0        = nullptr;
    int             heads_np                = 0;
    int             hydmech_np              = 0;
    cusfloat*       hydrostiff_p0           = nullptr;
    cuscomplex*     intensities             = nullptr;
    cuscomplex*     mdrift                  = nullptr;
    cuscomplex*     qtf_wl                  = nullptr;
    cuscomplex*     qtf_bern                = nullptr;
    cuscomplex*     qtf_acc                 = nullptr;
    cuscomplex*     qtf_mom                 = nullptr;
    cuscomplex*     qtf_body_vel_x_fk       = nullptr;
    cuscomplex*     qtf_body_vel_y_fk       = nullptr;
    cuscomplex*     qtf_body_vel_z_fk       = nullptr;
    cuscomplex*     qtf_body_vel_x_raddif   = nullptr;
    cuscomplex*     qtf_body_vel_y_raddif   = nullptr;
    cuscomplex*     qtf_body_vel_z_raddif   = nullptr;
    cuscomplex*     qtf_body_vel_x_total    = nullptr;
    cuscomplex*     qtf_body_vel_y_total    = nullptr;
    cuscomplex*     qtf_body_vel_z_total    = nullptr;
    cuscomplex*     qtf_wl_rel_we           = nullptr;
    cuscomplex*     qtf_wl_we               = nullptr;
    cuscomplex*     qtf_wl_we_fk            = nullptr;
    cuscomplex*     qtf_wl_we_raddif        = nullptr;
    cuscomplex*     qtf_wl_we_total         = nullptr;
    cuscomplex*     panels_potential        = nullptr;
    cuscomplex*     potential_secord_force  = nullptr;
    cuscomplex*     raos                    = nullptr;
    cusfloat*       structural_mass_p0      = nullptr;
    cuscomplex*     sysmat                  = nullptr;
    cuscomplex*     sysmat_steady           = nullptr;
    cuscomplex*     wave_diffrac            = nullptr;
    cuscomplex*     wave_diffrac_p0         = nullptr;
    int             wave_exc_np             = 0;
    cuscomplex*     wave_exc_p0             = nullptr;

    // Declare class constructors and destructor
    SimulationData( ) = default;

    SimulationData(
                        int         bodies_np,
                        int         dofs_np,
                        int         heads_np,
                        int         rows_local_np,
                        int         rows_np,
                        MpiConfig*  mpi_config_in
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
    
};

#endif