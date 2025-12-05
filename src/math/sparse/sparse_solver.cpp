
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

// Include general usage libraries
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include "mkl.h"

// Include local modules
#include "sparse_solver.hpp"


SparseSolver::SparseSolver(CSRMatrix* sysmat_in, bool verbose_in)
{
    // Storage 
    this->sysmat = sysmat_in;
    this->verbose = verbose_in;

    // Intialize solver
    this->initialize();
}


void SparseSolver::initialize(void)
{
    // Define method name in a string in order to use it
    // along the screen print-outs of the method
    std::string method_name = "INITIALIZATION";
     
    // Setup PARDISO control parameters
    for (int i = 0; i<64; i++ )
    {
        this->iparm[i] = 0;
    }
    this->iparm[0] = 1;         // No solver default
    this->iparm[1] = 2;         // Fill-in reordering from METIS
    this->iparm[3] = 0;         // No iterative-direct algorithm
    this->iparm[4] = 0;         // No user fill-in reducing permutation
    this->iparm[5] = 0;         // Write solution into x
    this->iparm[6] = 0;         // Not in use
    this->iparm[7] = 2;         // Max numbers of iterative refinement step
    this->iparm[8] = 0;         // Not in use
    this->iparm[9] = 13;        // Perturb the pivot elements with 1E-13
    this->iparm[10] = 1;        // Use nonsymmetric permutation and scaling MPS
    this->iparm[11] = 0;        // Conjugate transposed/transpose solve
    this->iparm[12] = 1;        // Maximum weighted matching algorithm is switched-on (default for non-symmetric)
    this->iparm[13] = 0;        // Output: Number of perturbed pivots
    this->iparm[14] = 0;        // Not in use
    this->iparm[15] = 0;        // Not in use
    this->iparm[16] = 0;        // Not in use
    this->iparm[17] = -1;       // Output: Number of nonzeros in the factor LU
    this->iparm[18] = -1;       // Output: Mflops for LU factorization
    this->iparm[19] = 0;        // Output: Numbers of CG Iterations
    this->iparm[34] = 1;        // Set Zero indexing
    this->maxfct = 1;           // Maximum number of numerical factorizations
    this->mnum = 1;             // Which factorization to use
    this->msglvl = 0;           // Print statistical information
    this->mtype = 11;           // Matrix type. Real Unsymetric: 11

    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (int i = 0; i < 64; i++ )
    {
        this->pt[i] = 0;
    }

    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    std::string section_name = "Reordering and Symbolic Factorization";
    if ( this->verbose )
    {
        std::cout << method_name << ": " << section_name << std::endl;
    }

    int error = 0;
    this->phase = 11;
    PARDISO(
        this->pt,                   // Internal solver memory pointer
        &(this->maxfct),            // Maximum number of numerical factorizations.
        &(this->mnum),              // Which factorization to use
        &(this->mtype),             // Matrix type. (Symmetric, Unsymmetric, etc.)
        &(this->phase),             // Solver execution step
        &(this->sysmat->rows_np),   // Number of rows, the matrix is suposed to be squared
        this->sysmat->values,       // Vector containing the non-zero values
        this->sysmat->row_index_cum,// Vector containing the cumulative elements along the rows
        this->sysmat->col_index,    // Vector containing the column index 
        &(this->idum),              // Integer dummy
        &(this->nrhs),              // Number of right hand sides
        this->iparm,                // PARDISO control parameters
        &(this->msglvl),            // Print statistical information
        &(this->ddum),              // Double dummy
        &(this->ddum),              // Double dummy
        &error                      // Integer that contains the status code
        );
    
    if ( error != 0 )
    {
        std::cout << method_name << ": " << section_name << std::endl;
        std::cout << " - ERROR during symbolic factorization: " << error << std::endl;
        exit(101);
    }

    if ( this->verbose )
    {
        std::cout << method_name << ": " << section_name << " - Reordering completed ... " << std::endl;
        std::cout << method_name << ": " << section_name << " - Number of nonzeros in factors = " << iparm[17] << std::endl;
        std::cout << method_name << ": " << section_name << " - Number of factorization MFLOPS = " << iparm[18] << std::endl;
        std::cout << method_name << ": " << section_name << " -> Done" << std::endl;
    }
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    section_name = "Numerical factorization";
    if ( this->verbose )
    {
        std::cout << method_name << ": " << section_name << std::endl;
    }

    phase = 22;
    PARDISO(
            this->pt,                   // Internal solver memory pointer
            &(this->maxfct),            // Maximum number of numerical factorizations.
            &(this->mnum),              // Which factorization to use
            &(this->mtype),             // Matrix type. (Symmetric, Unsymmetric, etc.)
            &(this->phase),             // Solver execution step
            &(this->sysmat->rows_np),   // Number of rows, the matrix is suposed to be squared
            this->sysmat->values,       // Vector containing the non-zero values
            this->sysmat->row_index_cum,// Vector containing the cumulative elements along the rows
            this->sysmat->col_index,    // Vector containing the column index 
            &(this->idum),              // Integer dummy
            &(this->nrhs),              // Number of right hand sides
            this->iparm,                // PARDISO control parameters
            &(this->msglvl),            // Print statistical information
            &(this->ddum),              // Double dummy
            &(this->ddum),              // Double dummy
            &error                      // Integer that contains the status code
            );
    
    if ( error != 0 )
    {
        std::cout << method_name << ": " << section_name << std::endl;
        std::cout << " - ERROR during numerical factorization: " << error << std::endl;
        exit (102);
    }

    if ( this->verbose )
    {
        std::cout << method_name << ": " << section_name << " - Factorization completed ... " << std::endl;
        std::cout << method_name << ": " << section_name << " - Done " << std::endl;
    }

}


void SparseSolver::solve(double* rhs, double* sol)
{
    // Define method name in a string in order to use it
    // along the screen print-outs of the method
    std::string method_name = "SOLVE";

    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    std::string section_name = "Back substitution and iterative refinement";
    if ( this->verbose )
    {
        std::cout << method_name << ": " << section_name << std::endl;
    }

    // Configure solver to solve the system of equations
    this->phase = 33;

    // Configure solver to resolve the matrix without performing
    // any operation on it
    this->iparm[11] = 0;


    // Solve system of equation using the previously decomposed
    // system
    int error = 0;
    PARDISO(
            this->pt,                   // Internal solver memory pointer
            &(this->maxfct),            // Maximum number of numerical factorizations.
            &(this->mnum),              // Which factorization to use
            &(this->mtype),             // Matrix type. (Symmetric, Unsymmetric, etc.)
            &(this->phase),             // Solver execution step
            &(this->sysmat->rows_np),   // Number of rows, the matrix is suposed to be squared
            this->sysmat->values,       // Vector containing the non-zero values
            this->sysmat->row_index_cum,// Vector containing the cumulative elements along the rows    
            this->sysmat->col_index,    // Vector containing the column index 
            &(this->idum),              // Integer dummy
            &(this->nrhs),              // Number of right hand sides
            this->iparm,                // PARDISO control parameters
            &(this->msglvl),            // Print statistical information
            rhs,                        // Right hand side vector to solve the system of equations
            sol,                        // Solution of the system of equations
            &error                      // Integer that contains the status code
            );
    
    if ( error != 0 )
    {
        std::cout << method_name << ": " << section_name << std::endl;
        std::cout << " - ERROR during solution: " << error << std::endl;
        exit (103);
    }
        
}