
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
#include "../config.hpp"



template<typename T>    void    euler_local_to_global(
                                                                T           alpha,
                                                                T           beta,
                                                                T           gamma,
                                                                T*          rot_mat
                                                    );


template<typename T>    void    euler_local_to_global_disp(
                                                                T*          dofs_trans,
                                                                T*          dofs_rot,
                                                                cusfloat*   radius,
                                                                T*          displacement
                                                            );


template<typename T>    void    _rot_x(
                                                                T*          mat,
                                                                T           alpha
                                        );


template<typename T>    void    _rot_y(
                                                                T*          mat,
                                                                T           beta
                                        );


template<typename T>    void    _rot_z(
                                                                T*          mat,
                                                                T           gamma
                                        );


#include "euler_transforms.txx"
