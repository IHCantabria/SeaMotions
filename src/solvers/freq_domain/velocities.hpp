
#ifndef __fds_velocities_hpp
#define __fds_velocities_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_raddif_velocity_mat_steady(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    MLGCmpx*        vel_x_gp,
                                                    MLGCmpx*        vel_y_gp,
                                                    MLGCmpx*        vel_z_gp
                                            );


void    calculate_raddif_velocity_mat_steady_nlin(
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