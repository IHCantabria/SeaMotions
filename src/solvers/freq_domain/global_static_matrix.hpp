
#pragma once

// Include local modules
#include "../../hydrostatics.hpp"
#include "../../inout/input.hpp"


void    calculate_global_hydstiffness(
                                                Input*          input,
                                                Hydrostatics**  hydrostatics,
                                                cusfloat*       hydstiffness
                                    );


void    calculate_global_structural_mass(
                                                Input*          input,
                                                cusfloat*       structural_mass_p0
                                        );