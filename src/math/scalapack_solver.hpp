
#ifndef scalapack_solver_hpp__
#define scalapack_solver_hpp__

// Include general usage libraries
#include <iostream>
#include "mpi.h"

// Include general usage scientific libraries
#include "mkl.h"
#include "mkl_scalapack.h"
#include "mkl_blacs.h"
#include "mkl_pblas.h"

// Include local libraries
#include "../../src/math/math_interface.hpp"
#include "../../src/mpi_interface.hpp"


template <class T>
class ScalapackSolver
{
public:
    // Declare public variables
    MKL_INT* descA;
    MKL_INT* descB;
    MKL_INT iam = -1;
    MKL_INT info;
    MKL_INT izero = 0;
    char layout='R'; // Block cyclic, Row major processor mapping
    MKL_INT lddA;
    MKL_INT proc_rank;
    MKL_INT num_block_size;
    MKL_INT num_cols_rhs = 1;
    MKL_INT num_cols_local;
    MKL_INT num_procs;
    MKL_INT num_procs_blas;
    MKL_INT num_procs_col;
    MKL_INT num_procs_row;
    MKL_INT num_rows;
    MKL_INT num_rows_local;
    MPI_Comm rhs_comm;
    MPI_Group rhs_group;
    MKL_INT start_col;
    MKL_INT start_row;
    MKL_INT end_col;
    MKL_INT end_row;
    MKL_INT ictxt, myrow, mycol;
    MKL_INT zero = 0;

    // Declare Constructors
    ScalapackSolver(MKL_INT num_rows, MKL_INT num_procs, MKL_INT proc_rank);

    // Declare Class Methdos
    void GenerateRhsComm(void);
    T* GetGlobalRhs(T* subrhs, T* sol_vec);
    void Initialize(void);
    void Solve(T* subsysmat, T* subrhs);
};


template <class T>
void ScalapackSolver<T>::GenerateRhsComm(void)
{
    MKL_INT* cols_position = NULL;
    cols_position = new MKL_INT[num_procs];
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&mycol, 1, MPI_INT, cols_position, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MKL_INT new_ranks[num_procs_row];

    MKL_INT count = 0;
    for (MKL_INT i=0; i<num_procs; i++)
    {
        if (cols_position[i] == 0)
        {
            new_ranks[count] = i;
            count++;
        }
    }

    // Construct a group containing all of the prime ranks in world_group
    MPI_Group_incl(world_group, num_procs_row, new_ranks, &rhs_group);
    MPI_Comm_create_group(MPI_COMM_WORLD, rhs_group, 0, &rhs_comm);
    
    // Delete matrixes
    delete [] cols_position;
}


template <class T>
T* ScalapackSolver<T>::GetGlobalRhs(T* subrhs, T* sol_vec)
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (rhs_comm != MPI_COMM_NULL)
    {
        cmpi_gather<T>( subrhs, num_rows_local, sol_vec, num_rows_local, 0, rhs_comm );
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return sol_vec;
}


template <class T>
void ScalapackSolver<T>::Initialize(void)
{
    // Initialize BLACS
    blacs_pinfo(&iam, &num_procs_blas) ; // BLACS rank and world size
    blacs_get(&zero, &zero, &ictxt ); // -> Create context
    blacs_gridinit(&ictxt, &layout, &num_procs_row, &num_procs_col); // Context -> Initialize the grid
    blacs_gridinfo(&ictxt, &num_procs_row, &num_procs_col, &myrow, &mycol); // Context -> Context grid info (# procs row/col, current procs row/col)

    // Compute size of the local matrixes
    num_rows_local = numroc(&num_rows, &num_block_size, &myrow, &izero, &num_procs_row); // My proc -> row of local A
    num_cols_local = numroc(&num_rows, &num_block_size, &mycol, &izero, &num_procs_col); // My proc -> col of local A

    // Create matrix descriptor 
    descA = new MKL_INT [9];
    descB = new MKL_INT [9];
    lddA = num_rows_local > 1 ? num_rows_local : 1;
    descinit_(descA, &num_rows, &num_rows, &num_block_size, &num_block_size, &izero, &izero, &ictxt, &lddA, &info);
    descinit_(descB, &num_rows, &num_cols_rhs, &num_block_size, &num_block_size, &izero, &izero, &ictxt, &lddA, &info);

    // descA[0] = 1; // descriptor type
    // descA[1] = ictxt; // blacs context
    // descA[2] = num_rows; // global number of rows
    // descA[3] = num_rows; // global number of columns
    // descA[4] = num_block_size; // row block size
    // descA[5] = num_block_size; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    // descA[6] = myrow; // initial process row(DEFINED 0)
    // descA[7] = mycol; // initial process column (DEFINED 0)
    // descA[8] = num_rows_local; // leading dimension of local array

    // descB[0] = 1; // descriptor type
    // descB[1] = ictxt; // blacs context
    // descB[2] = num_rows; // global number of rows
    // descB[3] = 1; // global number of columns
    // descB[4] = num_block_size; // row block size
    // descB[5] = num_block_size; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    // descB[6] = myrow; // initial process row(DEFINED 0)
    // descB[7] = mycol; // initial process column (DEFINED 0)
    // descB[8] = num_rows_local; // leading dimension of local array

    // Distribute matrix
    MKL_INT* num_rows_div = NULL;
    MKL_INT* num_cols_div = NULL;
    MKL_INT* rows_position = NULL;
    MKL_INT* cols_position = NULL;
    if (proc_rank == 0)
    {
        num_rows_div = new MKL_INT[num_procs];
        num_cols_div = new MKL_INT[num_procs];
        rows_position = new MKL_INT[num_procs];
        cols_position = new MKL_INT[num_procs];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&num_rows_local, 1, MPI_INT, num_rows_div, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&num_cols_local, 1, MPI_INT, num_cols_div, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&myrow, 1, MPI_INT, rows_position, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&mycol, 1, MPI_INT, cols_position, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate variables
    MKL_INT* start_cols = NULL;
    MKL_INT* start_rows = NULL;
    if (proc_rank == 0)
    {
        start_cols = new MKL_INT [num_procs];
        start_rows = new MKL_INT [num_procs];

        for (MKL_INT i=0; i<num_procs; i++)
        {
            // Calculate rows
            start_rows[i] = 1;
            for (MKL_INT j=0; j<rows_position[i]; j++)
            {
                start_rows[i] += num_rows_div[j];
            }

            // Calculate cols
            start_cols[i] = 1;
            for (MKL_INT j=0; j<cols_position[i]; j++)
            {
                start_cols[i] += num_cols_div[j];
            }
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(start_rows, 1, MPI_INT, &start_row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(start_cols, 1, MPI_INT, &start_col, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    this->end_row = this->start_row + this->num_rows_local;
    this->end_row = ( this->end_row > this->num_rows ) ? this->num_rows : this->end_row;
    this->end_col = this->start_col + this->num_cols_local;
    this->end_col = ( this->end_col > this->num_rows ) ? this->num_rows : this->end_col;

    // Delete vectors
    if (proc_rank == 0)
    {
        delete [] num_rows_div;
        delete [] num_cols_div;
        delete [] rows_position;
        delete [] cols_position;
        delete [] start_cols;
        delete [] start_rows;
    }

    this->GenerateRhsComm();
}


template <class T>
ScalapackSolver<T>::ScalapackSolver(MKL_INT num_rows_inc, MKL_INT num_procs_inc, MKL_INT proc_rank_inc)
{
    // Assing variables
    proc_rank = proc_rank_inc;
    num_rows = num_rows_inc;
    num_procs = num_procs_inc;

    // Calculate processors per row and column
    if (num_procs < 1)
    {
        throw "Number of processors cannot be less than 1\n";
    }
    else if (num_procs == 1)
    {
        num_procs_row = 1;
        num_procs_col = 1;
    }
    else if ( num_procs > num_rows )
    {
        throw "ERROR SCALAPACK:\n More processors than rows\n";
    }
    else if (num_procs < 4)
    {
        num_procs_row = 1;
        num_procs_col = num_procs;
    }
    else
    {
        num_procs_row = 1;
        num_procs_col = num_procs;
    }

    num_block_size = num_rows/num_procs_col + num_rows%num_procs_col;
    num_block_size = num_block_size > 0 ? num_block_size: 1;

    // Initalize solver
    this->Initialize();

}


template <class T>
void ScalapackSolver<T>::Solve(T* subsysmat, T* subrhs)
{
    MKL_INT startrow = 1;
    MKL_INT startcol = 1;
    MKL_INT* ipiv = new MKL_INT[num_rows_local + num_block_size];
    pgesv<T>(&num_rows, &num_cols_rhs, subsysmat, &startrow, &startcol, descA, ipiv, subrhs, &startrow, &startrow, descB, &info);

    if ( info != 0 )
    {
        std::cerr << "ERROR - Scalapack Solver" << std::endl;
        std::cout << "Scalapack solver finished abnormally with code: " << info << std::endl;
        throw std::runtime_error( "" );
    }
}

#endif // scalapack_real_solver_hpp__