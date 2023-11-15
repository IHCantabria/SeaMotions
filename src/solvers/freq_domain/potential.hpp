
#ifndef __fds_potential_hpp
#define __fds_potential_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/mpi_config.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_influence_potmat_steady(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat*       field_points,
                                                int             field_points_np,
                                                cuscomplex*     inf_pot_mat
                                            );


void    calculate_influence_potmat(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     inf_pot_steady,
                                                cusfloat*       field_points,
                                                int             field_points_np,
                                                cuscomplex*     inf_pot_total
                                    );


void    calculate_potpanel_raddif_lin(
                                                Input*          input,
                                                cuscomplex*     inf_pot_mat,
                                                int             rows_np,
                                                int             cols_np,
                                                int             start_col,
                                                cuscomplex*     sources,
                                                cuscomplex*     panel_pot
                                    );


void    calculate_potpanel_raddif_nlin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cuscomplex*     sources,
                                                cusfloat        ang_freq
                                        );


void    calculate_potpanel_total_lin(
                                                Input*          input,
                                                MpiConfig*      mpi_config,
                                                MeshGroup*      mesh_gp,
                                                cusfloat        ang_freq,
                                                cuscomplex*     intensities,
                                                cuscomplex*     pot_steady_sysmat,
                                                cuscomplex*     pot_sysmat,
                                                cusfloat*       pot_fp,
                                                int*            pot_fp_cnp,
                                                int             pot_fp_nb,
                                                cuscomplex*     raos,
                                                cuscomplex*     potpanel_total
                                    );

#endif