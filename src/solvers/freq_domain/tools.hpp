
#ifndef __fds_tools_hpp
#define __fds_tools_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_fields_raddif_lin(
                                                Input*          input,
                                                cuscomplex*     intensities,
                                                MLGCmpx*        pot_gp
                                    );


void    define_gauss_points_diffrac_panels(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        mat_gp
                                            );


void    define_gauss_points_wl(
                                                Input*      input,
                                                MeshGroup*  mesh_gp,
                                                MLGCmpx*    mat_gp
                                );


#endif