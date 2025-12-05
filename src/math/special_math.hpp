
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


cusfloat    besseli0( 
                                    cusfloat x 
                    );

cusfloat    besseli1( 
                                    cusfloat x
                    );

cusfloat    besselj0( 
                                    cusfloat x
                    );

cusfloat    besselj1( 
                                    cusfloat x
                    );

cusfloat    besseljn_cos_int( 
                                    cusfloat alpha, 
                                    cusfloat beta, 
                                    cusfloat nu
                           );

cusfloat    besseljn_cos_kernel(
                                    cusfloat alpha,
                                    cusfloat beta,
                                    cusfloat nu,
                                    cusfloat x
                                );

cuscomplex  besseljn_expi_int( 
                                    cusfloat alpha,
                                    cusfloat beta,
                                    cusfloat nu
                            );

cusfloat    besseljn_sin_int( 
                                    cusfloat alpha, 
                                    cusfloat beta, 
                                    cusfloat nu
                            );

cusfloat    besseljn_sin_kernel(
                                    cusfloat alpha,
                                    cusfloat beta,
                                    cusfloat nu,
                                    cusfloat x
                                );

cusfloat    besselk0( 
                                    cusfloat x 
                    );

cusfloat    besselk1( 
                                    cusfloat x
                    );

cusfloat    bessely0( 
                                    cusfloat x
                    );

cusfloat    bessely1( 
                                    cusfloat x
                    );

cusfloat    ep_n( 
                                    int n 
                );

cusfloat    expint_i( 
                                    cusfloat x
                    );

cusfloat    legendre_poly_raw( 
                                    int n,
                                    cusfloat x 
                            );

cusfloat    legendre_poly_der_raw( 
                                    int n,
                                    cusfloat x
                                );

cusfloat    polynomial_f0( 
                                    cusfloat x
                        );

cusfloat    polynomial_f1( 
                                    cusfloat x 
                        );

cusfloat    polynomial_th0( 
                                    cusfloat x
                            );

cusfloat    polynomial_th1( 
                                    cusfloat x
                        );

cusfloat    psi_fun( 
                                    int n
                    );

cusfloat    rational_fraction_f0( 
                                    cusfloat x
                                );

cusfloat    rational_fraction_f1( 
                                    cusfloat x
                                );

cusfloat    rational_fraction_th0( 
                                    cusfloat x
                                );

cusfloat    rational_fraction_th1( 
                                    cusfloat x
                                );

cusfloat    struve0( 
                                    cusfloat x
                    );

cusfloat    struve1( 
                                    cusfloat x
                    );
