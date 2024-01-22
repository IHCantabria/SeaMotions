
#ifndef __kochin_hpp
#define __kochin_hpp

// Include local modules
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
#endif