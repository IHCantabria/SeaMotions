
#ifndef __qtf_hpp
#define __qtf_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_second_order_force(
                                        Input*          input,
                                        MpiConfig*      mpi_config,
                                        MeshGroup*      mesh_gp,
                                        cuscomplex*     mdrift_rel_we_i,
                                        cuscomplex*     mdrift_rel_we_j,
                                        cuscomplex*     raos_i,
                                        cuscomplex*     raos_j,
                                        cuscomplex*     vel_x_i,
                                        cuscomplex*     vel_y_i,
                                        cuscomplex*     vel_z_i,
                                        cuscomplex*     vel_x_j,
                                        cuscomplex*     vel_y_j,
                                        cuscomplex*     vel_z_j,
                                        cusfloat        ang_freq_i,
                                        cusfloat        ang_freq_j,
                                        cuscomplex*     qtf_values
                                    );


#endif