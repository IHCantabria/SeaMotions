
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


template <class T>
class ScalapackSolver
{
public:
    // Declare public variables
    int* descA;
    int* descB;
    int iam = -1;
    int info;
    int izero = 0;
    char layout='R'; // Block cyclic, Row major processor mapping
    int lddA;
    int proc_rank;
    int num_block_size;
    int num_cols_rhs = 1;
    int num_cols_local;
    int num_procs;
    int num_procs_blas;
    int num_procs_col;
    int num_procs_row;
    int num_rows;
    int num_rows_local;
    MPI_Comm rhs_comm;
    MPI_Group rhs_group;
    int start_col;
    int start_row;
    int ictxt, myrow, mycol;
    int zero = 0;

    // Declare Constructors
    ScalapackSolver(int num_rows, int num_procs, int proc_rank);

    // Declare Class Methdos
    void GenerateRhsComm(void);
    T* GetGlobalRhs(T* subrhs, T* sol_vec);
    void Initialize(void);
    void Solve(T* subsysmat, T* subrhs);
};


template <class T>
void ScalapackSolver<T>::GenerateRhsComm(void)
{
    int* cols_position = NULL;
    cols_position = new int[num_procs];
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather(&mycol, 1, MPI_INT, cols_position, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    int new_ranks[num_procs_row];

    int count = 0;
    for (int i=0; i<num_procs; i++)
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
        MPI_Gather(subrhs, num_rows_local, MPI_DOUBLE, sol_vec, num_rows_local, MPI_DOUBLE, 0, rhs_comm);
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
    descA = new int [9];
    descB = new int [9];
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
    int* num_rows_div = NULL;
    int* num_cols_div = NULL;
    int* rows_position = NULL;
    int* cols_position = NULL;
    if (proc_rank == 0)
    {
        num_rows_div = new int[num_procs];
        num_cols_div = new int[num_procs];
        rows_position = new int[num_procs];
        cols_position = new int[num_procs];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&num_rows_local, 1, MPI_INT, num_rows_div, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&num_cols_local, 1, MPI_INT, num_cols_div, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&myrow, 1, MPI_INT, rows_position, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&mycol, 1, MPI_INT, cols_position, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Allocate variables
    int* start_cols = NULL;
    int* start_rows = NULL;
    if (proc_rank == 0)
    {
        start_cols = new int [num_procs];
        start_rows = new int [num_procs];

        for (int i=0; i<num_procs; i++)
        {
            // Calculate rows
            start_rows[i] = 1;
            for (int j=0; j<rows_position[i]; j++)
            {
                start_rows[i] += num_rows_div[j];
            }

            // Calculate cols
            start_cols[i] = 1;
            for (int j=0; j<cols_position[i]; j++)
            {
                start_cols[i] += num_cols_div[j];
            }
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scatter(start_rows, 1, MPI_INT, &start_row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(start_cols, 1, MPI_INT, &start_col, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Delete vectors
    if (proc_rank == 0)
    {
        for (int i=0; i<num_procs; i++)
        {
            std::cout << "Processor Rank: " << i << " - Start_Row: " << start_rows[i] << " - Start_Col: " << start_cols[i] << "\n";
        }

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
ScalapackSolver<T>::ScalapackSolver(int num_rows_inc, int num_procs_inc, int proc_rank_inc)
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
    else if (num_procs < 4)
    {
        num_procs_row = 1;
        num_procs_col = num_procs;
    }
    else
    {
        num_procs_row = num_procs/2;
        num_procs_col = num_procs/2;
    }

    // Calculate block size
    num_block_size = num_rows/num_procs_row + num_rows%num_procs_row;
    
    // Initalize solver
    this->Initialize();

}


template <class T>
void ScalapackSolver<T>::Solve(T* subsysmat, T* subrhs)
{
    int startrow = 1;
    int startcol = 1;
    int* ipiv = new int[num_rows_local + num_block_size];
    pgesv<T>(&num_rows, &num_cols_rhs, subsysmat, &startrow, &startcol, descA, ipiv, subrhs, &startrow, &startrow, descB, &info);

    std::cout << "ScaLAPACK solution status: " << info << "\n";
}

#endif // scalapack_real_solver_hpp__