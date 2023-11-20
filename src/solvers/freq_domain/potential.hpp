
#ifndef __fds_potential_hpp
#define __fds_potential_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_influence_potmat_steady(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        pot_gp
                                            );


void    calculate_influence_potmat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        pot_gp
                                    );


void    calculate_potpanel_raddif_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     intensities,
                                                cusfloat        ang_freq
                                        );


void    calculate_potpanel_total_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     intensities,
                                                cuscomplex*     raos,
                                                MLGCmpx*        pot_gp,
                                                cuscomplex*     potpanel_total
                                    );

#endif