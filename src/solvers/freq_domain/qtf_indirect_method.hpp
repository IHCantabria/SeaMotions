
#ifndef __qtf_indirect_method_hpp
#define __qtf_indirect_method_hpp

// Include local modules
#include "../../inout/input.hpp"
#include "../../mesh/mesh_group.hpp"


void    calculate_secord_force_indirect(
                                            Input*      input,
                                            MeshGroup*  mesh_gp,
                                            cusfloat    ang_freq_i,
                                            cusfloat    ang_freq_j,
                                            bool        is_diff,
                                            cuscomplex* froude_krylov,
                                            cuscomplex* body_force,
                                            cuscomplex* fs_near_field,
                                            cuscomplex* fs_far_field,
                                            cuscomplex* secord_force_total
                                        );

#endif