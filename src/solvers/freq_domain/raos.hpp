
/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

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