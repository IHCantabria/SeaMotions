
#ifndef __qtf_indirect_method_hpp
#define __qtf_indirect_method_hpp

// Include local modules
#include "../../containers/simulation_data.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"
#include "../../waves/wave_dispersion_so.hpp"


void        calculate_qtf_indirect_body_term(
                                                Input*          input,
                                                MeshGroup*      mesh_gp,
                                                int             freq_pos_i,
                                                int             freq_pos_j,
                                                int             qtf_type,
                                                SimulationData* sim_data,
                                                MLGCmpx*        body_gp,
                                                MLGCmpx*        wl_gp,
                                                cuscomplex*     qtf_body_force
                                            );


void        calculate_qtf_indirect_fs_near_term(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    int             freq_pos_i,
                                                    int             freq_pos_j,
                                                    int             qtf_type,
                                                    SimulationData* sim_data,
                                                    MLGCmpx*        body_gp,
                                                    cuscomplex*     qtf_fs_force
                                                );


void        calculate_qtf_indirect_fs_far_term(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    int             freq_pos_i,
                                                    int             freq_pos_j,
                                                    int             qtf_type,
                                                    SimulationData* sim_data,
                                                    MLGCmpx*        body_gp,
                                                    cuscomplex*     qtf_fs_force
                                                );


cuscomplex  calculate_r0_integral(
                                                    cusfloat            R,
                                                    WaveDispersionSO*   wdso,
                                                    int                 l_order,
                                                    int                 qtf_type
                                );


cuscomplex  calculate_r1_integral(
                                                    cusfloat            R,
                                                    WaveDispersionSO*   wdso,
                                                    int                 l_order,
                                                    int                 qtf_type
                                );


cuscomplex  calculate_theta_integral(
                                                    Input*      input,
                                                    cusfloat    beta,
                                                    int         l_order,
                                                    int         qtf_type,
                                                    int         theta_type,
                                                    cuscomplex* kochin_cos_pert_j,
                                                    cuscomplex* kochin_sin_pert_j,
                                                    cuscomplex* kochin_cos_rad_i,
                                                    cuscomplex* kochin_sin_rad_i,
                                                    cuscomplex* kochin_cos_rad_j,
                                                    cuscomplex* kochin_sin_rad_j,
                                                    cuscomplex* body_force
                                    );


void        calculate_secord_force_indirect(
                                                    Input*          input,
                                                    MeshGroup*      mesh_gp,
                                                    int             freq_pos_i,
                                                    int             freq_pos_j,
                                                    int             qtf_type,
                                                    MLGCmpx*        body_gp,
                                                    MLGCmpx*        wl_gp,
                                                    SimulationData* sim_data
                                            );

#endif