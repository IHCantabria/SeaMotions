
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

// Include local modules
#include "raos.hpp"


void   calculate_rao_disp(
                                cuscomplex*     rao_cog,
                                const cusfloat* radius,
                                cuscomplex*     rao_disp
                            )
{
    rao_disp[0] = rao_cog[0] + rao_cog[4] * radius[2] - rao_cog[5] * radius[1];
    rao_disp[1] = rao_cog[1] + rao_cog[5] * radius[0] - rao_cog[3] * radius[2];
    rao_disp[2] = rao_cog[2] + rao_cog[3] * radius[1] - rao_cog[4] * radius[0];
    rao_disp[3] = rao_cog[3];
    rao_disp[4] = rao_cog[4];
    rao_disp[5] = rao_cog[5];
}