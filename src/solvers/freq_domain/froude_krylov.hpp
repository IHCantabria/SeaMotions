
#ifndef __fds_froude_krylov_hpp
#define __fds_froude_krylov_hpp

// Include local modules
#include "../../inout/input.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_froude_krylov_fo(
                                    Input*          input,
                                    MpiConfig*      mpi_config,
                                    MeshGroup*      mesh_gp,
                                    cusfloat        ang_freq,
                                    cuscomplex*     froude_krylov
                                );

void    calculate_froude_krylov_so(
                                    Input*          input,
                                    MeshGroup*      mesh_gp,
                                    cusfloat        ang_freq_i,
                                    cusfloat        ang_freq_j,
                                    int             qtf_type,
                                    cuscomplex*     froude_krylov
                                );

#endif