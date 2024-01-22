
#ifndef __kochin_hpp
#define __kochin_hpp

// Include local modules
#include "../containers/simulation_data.hpp"
#include "../inout/input.hpp"
#include "../interfaces/kochin_interface.hpp"
#include "../mesh/mesh_group.hpp"


void    calculate_kochin_coefficients(
                                        Input*              input,
                                        MeshGroup*          mesh_gp,
                                        KochinInterface*    kochin,
                                        cuscomplex*         sources,
                                        cuscomplex*         cos_coeff,
                                        cuscomplex*         sin_coeff
                                    );


void    calculate_kochin_pert_coeffs(
                                        Input*          input,
                                        MeshGroup*      mesh_gp,
                                        int             freq_pos,
                                        SimulationData* sim_data
                                    );


void    calculate_kochin_rad_coeffs(
                                        Input*          input,
                                        MeshGroup*      mesh_gp,
                                        int             freq_pos,
                                        SimulationData* sim_data
                                    );

#endif