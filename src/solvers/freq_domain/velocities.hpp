
#ifndef __fds_velocities_hpp
#define __fds_velocities_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../interfaces/gwfdx_interface.hpp"
#include "../../interfaces/gwfdy_interface.hpp"
#include "../../interfaces/gwfdz_interface.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_raddif_velocity_mat(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                GWFDxInterface* gwfdx_interf,
                                                GWFDyInterface* gwfdy_interf,
                                                GWFDzInterface* gwfdz_interf,
                                                MLGCmpx*        vel_x_gp,
                                                MLGCmpx*        vel_y_gp,
                                                MLGCmpx*        vel_z_gp
                                        );


void    calculate_raddif_velocity_mat_steady(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        vel_x_gp,
                                                MLGCmpx*        vel_y_gp,
                                                MLGCmpx*        vel_z_gp
                                            );


void    calculate_raddif_velocity_mat_wave(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                MLGCmpx*        vel_x_gp,
                                                MLGCmpx*        vel_y_gp,
                                                MLGCmpx*        vel_z_gp
                                            );


void    calculate_velocities_total(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     intensities,
                                                cuscomplex*     raos,
                                                MLGCmpx*        vel_x_gp,
                                                MLGCmpx*        vel_y_gp,
                                                MLGCmpx*        vel_z_gp,
                                                cuscomplex*     vel_x_total,
                                                cuscomplex*     vel_y_total,
                                                cuscomplex*     vel_z_total
                                    );

#endif