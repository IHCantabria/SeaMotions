
#ifndef __fds_raos_hpp
#define __fds_raos_hpp

// Include local modules
#include "../../config.hpp"
#include "../../inout/input.hpp"


void    calculate_raos(
                            Input*          input,
                            cusfloat*       structural_mass,
                            cusfloat*       added_mass,
                            cusfloat*       damping_rad,
                            cusfloat*       hydstiffness,
                            cuscomplex*     wave_diffrac,
                            cuscomplex*     froude_krylov,
                            cusfloat        ang_freq,
                            cuscomplex*     rao
                        );

#endif