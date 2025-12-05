
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


template<int NV, int NGP>   inline  void    jacobi_det_2d_gp(
                                                                    cusfloat*   xn,
                                                                    cusfloat*   yn,
                                                                    cusfloat*   jac_det
                                                            );


template<int NV, int NGP>   inline  void    jacobi_mat_2d_gp(
                                                                    cusfloat*   xn,
                                                                    cusfloat*   yn,
                                                                    cusfloat*   jac_mat
                                                            );

                                                    
template<int NV, int NGP>   inline  void    shape_fcn_2d_gp( 
                                                                    cusfloat*   N
                                                            );


template<int NV, int NGP>   inline  void    shape_fcn_der_2d_gp(
                                                                    cusfloat*   dNdxi,
                                                                    cusfloat*   dNdeta
                                                                );


template<int NV, int NGP>   inline  void    shape_fcn_deta_2d_gp( 
                                                                    cusfloat*   N
                                                                );


template<int NV, int NGP>   inline  void    shape_fcn_dxi_2d_gp( 
                                                                    cusfloat*   N
                                                                );


#include "topology_gauss.txx"