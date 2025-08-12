
#ifndef __fds_gf_intensities_hpp
#define __fds_gf_intensities_hpp

// Include local modules
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../interfaces/grfdn_interface.hpp"
#include "../../interfaces/gwfcns_interface_t.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_gf_intensity_sysmat(
                                                    Input*                      input,
                                                    MpiConfig*                  mpi_config,
                                                    SclCmpx*                    scl,
                                                    MeshGroup*                  mesh_gp,
                                                    GWFcnsInterfaceT<NUM_GP2>&  gwf_interf,
                                                    cusfloat                    w,
                                                    MLGCmpx*                    sources_mlg,
                                                    MLGCmpx*                    potential_mlg,
                                                    cuscomplex*                 sources_int
                                        );


void    calculate_gf_intensity_steady_sysmat_lin(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    cuscomplex*     sysmat
                                                );


void    calculate_gf_intensity_steady_sysmat_nlin(
                                                    Input*          input,
                                                    SclCmpx*        scl,
                                                    MeshGroup*      mesh_gp,
                                                    GRFDnInterface* grf_interf,
                                                    cuscomplex*     sysmat
                                                );

#endif