
#ifndef __fds_tools_hpp
#define __fds_tools_hpp

// Include local modules
#include "../../config.hpp"
#include "../../containers/matlin_group.hpp"
#include "../../inout/input.hpp"


void    calculate_fields_raddif_lin(
                                        Input*          input,
                                        cuscomplex*     intensities,
                                        MLGCmpx*        pot_gp
                                    );


#endif