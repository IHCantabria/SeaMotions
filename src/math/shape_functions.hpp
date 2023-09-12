
#ifndef __shape_functions_hpp
#define __shape_functions_hpp

// Include local modules
#include "../config.hpp"


int         dofs_rectangular_region(
                                        int         poly_order
                                    );

int         dofs_triangular_region(
                                        int         poly_order
                                    );

cusfloat    shape_functions(
                                        int         p_order,
                                        int         q_order,
                                        cusfloat    ,
                                        cusfloat    
                           );

#endif