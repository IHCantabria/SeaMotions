
#ifndef __fds_gf_intensities_hpp
#define __fds_gf_intensities_hpp

// Include local modules
#include "../../inout/input.hpp"
#include "../../interfaces/grfdn_interface.hpp"
#include "../../interfaces/gwfdn_interface.hpp"
#include "../../math/scalapack_solver.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_gf_intensity_sysmat(
                                                Input*          input,
                                                SclCmpx*        scl,
                                                MeshGroup*      mesh_gp,
                                                GWFDnInterface* gwf_interf,
                                                cusfloat        w,
                                                cuscomplex*     sysmat_steady,
                                                cuscomplex*     sysmat,
                                                cuscomplex*     sources_int
                                        );


void    calculate_gf_intensity_steady_sysmat(
                                                Input*          input,
                                                SclCmpx*        scl,
                                                MeshGroup*      mesh_gp,
                                                GRFDnInterface* grf_interf,
                                                cuscomplex*     sysmat
                                            );

#endif