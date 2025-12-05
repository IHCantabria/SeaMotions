
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

// Import local modules
#include "sparse_containers.hpp"


class SparseSolver
{
public:
    // Declare class attributes
    double      ddum;           // Double dummy
    int         idum;           // Integer dummy
    int         iparm[64];      // PARDISO control parameters
    int         maxfct;         // Maximum number of numerical factorizations.
    int         mnum;           // Which factorization to use
    int         msglvl;         // Print statistical information
    int         mtype;          // Matrix type. (Symmetric, Unsymmetric, etc.)
    int         nrhs = 1;       // Number of right hand sides. By default 1
    int         phase;          // Solver execution step
    int*        pt[64];         // Internal pointer used by PARDISO solver
    CSRMatrix*  sysmat=nullptr; // System matrix used to solve the system of equations
    bool        verbose;        // Signal to control the print out to screen

    // Declare class constructors and destructors
    SparseSolver(CSRMatrix* sysmat_in, bool verbose_in);

    // Declare class methods
    void initialize(void);
    void solve(double* rhs, double* sol);
};