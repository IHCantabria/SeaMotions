
#ifndef __fds_tools_hpp
#define __fds_tools_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void        calculate_field_point_rot(
                                                cuscomplex*     raos_trans,
                                                cuscomplex*     raos_rot,
                                                cusfloat*       field_point,
                                                cusfloat*       cog,
                                                cuscomplex*     point_disp
                                    );


void        calculate_field_point_vel_rot(
                                                cuscomplex*     raos_trans,
                                                cuscomplex*     raos_rot,
                                                cusfloat*       field_point,
                                                cusfloat*       cog,
                                                cusfloat        ang_freq,
                                                cuscomplex*     point_disp
                                        );


std::string compose_dof_path( 
                                                std::string     base_path,
                                                int             dofs_num,
                                                int             ang_freq_num
                            );


void        define_gauss_points_diffrac_panels(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        mat_gp
                                            );


void        define_gauss_points_wl(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                MLGCmpx*        mat_gp
                                );


void        storage_radiation_potential( 
                                                std::string     dof_base_path,
                                                int             field_points_np,
                                                cuscomplex*     rad_potential
                                        );


#endif