
#ifndef __fds_hydromechanics_hpp
#define __fds_hydromechanics_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../interfaces/hmf_interface.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_hydromechanic_coeffs_lin( 
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panels_pot,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad
                                            );


void    calculate_hydromechanic_coeffs_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                HMFInterface*   hmf_interf,
                                                cusfloat        ang_freq,
                                                cusfloat*       added_mass,
                                                cusfloat*       damping_rad
                                            );

#endif