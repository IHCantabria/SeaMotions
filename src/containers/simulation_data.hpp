
#ifndef __simulation_data_hpp
#define __simulation_data_hpp

// Include local modules
#include "../config.hpp"
#include "mpi_config.hpp"


struct SimulationData
{
private:
    // Declare private class attributes
    bool            _is_mdrift          = false;
    MpiConfig*      _mpi_config         = nullptr;

public:
    // Declare class attributes
    cusfloat*       added_mass          = nullptr;
    cusfloat*       added_mass_p0       = nullptr;
    cusfloat*       damping_rad         = nullptr;
    cusfloat*       damping_rad_p0      = nullptr;
    cuscomplex*     froude_krylov       = nullptr;
    cuscomplex*     froude_krylov_p0    = nullptr;
    int             hydmech_np          = 0;
    cusfloat*       hydrostiff_p0       = nullptr;
    cuscomplex*     intensities         = nullptr;
    cuscomplex*     mdrift_rel_we       = nullptr;
    cuscomplex*     mdrift_we           = nullptr;
    cuscomplex*     mdrift_we_pot_total = nullptr;
    cuscomplex*     panels_potential    = nullptr;
    cuscomplex*     raos                = nullptr;
    cusfloat*       structural_mass_p0  = nullptr;
    cuscomplex*     sysmat              = nullptr;
    cuscomplex*     sysmat_steady       = nullptr;
    cuscomplex*     wave_diffrac        = nullptr;
    cuscomplex*     wave_diffrac_p0     = nullptr;
    int             wave_exc_np         = 0;
    cuscomplex*     wave_exc_p0         = nullptr;

    // Declare class constructors and destructor
    SimulationData( ) = default;

    SimulationData(
                        int         bodies_np,
                        int         dofs_np,
                        int         heads_np,
                        int         rows_local_np,
                        int         cols_local_np,
                        int         rows_np,
                        MpiConfig*  mpi_config_in
                    );

    ~SimulationData( 
                        void
                    );

    // Declare public class methods
    void    add_mean_drift_data(
                                    int mdrift_np
                                );
    
};

#endif