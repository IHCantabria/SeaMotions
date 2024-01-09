
#ifndef __qtf_hpp
#define __qtf_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/simulation_data.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"

// Set module definitions
#define QTF_DIFF_CODE   0
#define QTF_SUM_CODE    1


// Declare module functions
cuscomplex  calculate_qtf_diff_term(
                                        cuscomplex c0,
                                        cuscomplex c1
                                    );

void        calculate_qtf_terms_force(
                                        Input*          input,
                                        MeshGroup*      mesh_gp,
                                        int             qtf_type,
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
                                        cuscomplex*     phi_2_force,
                                        cusfloat        ang_freq_i,
                                        cusfloat        ang_freq_j,
                                        cuscomplex*     qtf_values,
                                        cuscomplex*     qtf_wl,
                                        cuscomplex*     qtf_bern,
                                        cuscomplex*     qtf_acc,
                                        cuscomplex*     qtf_mom,
                                        MLGCmpx*        pot_gp,
                                        MLGCmpx*        vel_gp,
                                        bool            is_multi_head
                                    );


#endif