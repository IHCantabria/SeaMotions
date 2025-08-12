
#ifndef __fds_diffraction_hpp
#define __fds_diffraction_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
// #include "../../interfaces/hmf_interface.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_diffraction_forces_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     panel_pot,
                                                cusfloat        w,
                                                cuscomplex*     wave_diffrac,
                                                cuscomplex*     panel_pressure
                                        );


// void    calculate_diffraction_forces_nlin(
//                                                 Input*          input,
//                                                 MpiConfig*      mpi_config,
//                                                 MeshGroup*      mesh_gp,
//                                                 HMFInterface*   hmf_interf,
//                                                 cusfloat        w,
//                                                 cuscomplex*     wave_diffrac
//                                         );

#endif