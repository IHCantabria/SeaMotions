
#pragma once

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
    MKL_INT*    descA;
    MKL_INT*    descB;
    MPI_Comm    global_comm;
    MPI_Group   global_comm_gp;
    MKL_INT     iam = -1;
    MKL_INT     info;
    MKL_INT*    ipiv = nullptr;
    MKL_INT     izero = 0;
    char        layout='R'; // Block cyclic, Row major processor mapping
    MKL_INT     lddA;
    MKL_INT     proc_rank;
    MKL_INT     proc_root;
    MKL_INT     num_block_size;
    MKL_INT     num_block_size_rhs;
    MKL_INT     num_cols_rhs = 1;
    MKL_INT     num_cols_local;
    MKL_INT     num_procs;
    MKL_INT     num_procs_blas;
    MKL_INT     num_procs_col;
    MKL_INT     num_procs_row;
    MKL_INT     num_rows;
    MKL_INT     num_rows_local;
    MPI_Comm    rhs_comm;
    MPI_Group   rhs_group;
    MKL_INT     start_col;
    MKL_INT     start_col_0;
    MKL_INT     start_row;
    MKL_INT     start_row_0;
    MKL_INT     end_col;
    MKL_INT     end_col_0;
    MKL_INT     end_row;
    MKL_INT     end_row_0;
    MKL_INT     ictxt, myrow, mycol;
    MKL_INT     zero = 0;

    // Declare Constructors
    ScalapackSolver( ) = default;

    ScalapackSolver(
                        MKL_INT num_rows, 
                        MKL_INT num_cols_rhs,  
                        MKL_INT num_procs, 
                        MKL_INT proc_rank,
                        MKL_INT proc_root,
                        MPI_Comm global_comm_inc
                    );

    ~ScalapackSolver(  );

    // Declare Class Methdos
    void    Cond(               T*          subsysmat,
                                cusfloat&   cond 
                );

    void    GenerateRhsComm(    
                                void 
                            );

    T*      GetGlobalRhs(       
                                T* subrhs, 
                                T* sol_vec
                        );
    
    void    Initialize(
                                void
                    );

    void    Solve( 
                                T* subsysmat, 
                                T* subrhs 
                    );
};


template <class T>
void ScalapackSolver<T>::Cond( T* subsysmat_orig, cusfloat& cond )
{
    // Make a local copy of the matrix to avoid factorizing it
    int mat_size = num_rows_local *  num_cols_local;
    T* subsysmat = (T*)mkl_calloc( mat_size, sizeof(T), FLOATING_PRECISION );

    for ( int i=0; i<mat_size; i++ )
    {
        subsysmat[i] = subsysmat_orig[i];
    }

    // Define local variables
    cusfloat    anorm       = 0;
    char        norm        = '1';
    MKL_INT     startrow    = 1;
    MKL_INT     startcol    = 1;
    MKL_INT     lwork       = 4 * num_rows;
    T*          work        = (T*)mkl_calloc( lwork, sizeof(T), FLOATING_PRECISION );
    cusfloat*   iwork       = (cusfloat*)mkl_calloc( lwork, sizeof(cusfloat), FLOATING_PRECISION );

    // Clear pivot vector
    for ( int i=0; i< ( num_rows_local + num_block_size ); i++ )
    {
        this->ipiv[i] = 0;
    }

    // Calculate matrix norm ( 1-norm )
    MPI_Barrier( this->global_comm );
    anorm = plange<T>(&norm, &num_rows_local, &num_cols_local, subsysmat, &startrow, &startcol, descA, iwork);
    MPI_Barrier( this->global_comm );

    if ( info != 0 )
    {
        std::cerr << "ERROR - Scalapack Norm 1" << std::endl;
        std::cerr << "Scalapack norm 1 finished abnormally with code: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    // Calculate LU factorization
    MPI_Barrier( this->global_comm );
    pgetrf<T>(&num_rows_local, &num_rows_local, subsysmat, &startrow, &startcol, descA, ipiv, &info);
    MPI_Barrier( this->global_comm );

    if ( info != 0 )
    {
        std::cerr << "ERROR - Scalapack LU decomposition" << std::endl;
        std::cerr << "Scalapack LU decomposition finished abnormally with code: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    // Calculate Condition Number
    MPI_Barrier( this->global_comm );
    pgecon<T>(&norm, &num_rows_local, subsysmat, &startrow, &startcol, descA, &anorm, &cond, work, &lwork, iwork, &lwork, &info);
    MPI_Barrier( this->global_comm );

    if ( info != 0 )
    {
        std::cerr << "ERROR - Scalapack LU decomposition" << std::endl;
        std::cerr << "Scalapack LU decomposition finished abnormally with code: " << info << std::endl;
        throw std::runtime_error( "" );
    }

    cond = 1 / cond;

    // Clean heap local heap memory
    mkl_free( iwork );
    mkl_free( subsysmat);
    mkl_free( work );
}


template <class T>
void ScalapackSolver<T>::GenerateRhsComm(void)
{
    MKL_INT* cols_position = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );
    
    MPI_Barrier(this->global_comm);
    MPI_Allgather(&mycol, 1, MPI_INT, cols_position, 1, MPI_INT, this->global_comm);
    MPI_Barrier(this->global_comm);

    // Get the group of processes in global_comm
    MPI_Group world_group;
    MPI_Comm_group(this->global_comm, &world_group);

    // Set new ranks to compose the rhs_comm. In this case
    // will be the root process only due to the processors
    // layout configured for Scalapack Solver in this class
    MKL_INT new_ranks[num_procs_row];
    new_ranks[0] = this->proc_root;

    // Construct a group containing all of the prime ranks in world_group
    MPI_Group_incl(world_group, num_procs_row, new_ranks, &rhs_group);
    MPI_Comm_create_group(this->global_comm, rhs_group, 0, &rhs_comm);
    
    // Delete matrixes
    mkl_free( cols_position );
}


template <class T>
T* ScalapackSolver<T>::GetGlobalRhs(T* subrhs, T* sol_vec)
{
    MPI_Barrier(this->global_comm);
    if (rhs_comm != MPI_COMM_NULL)
    {
        cmpi_gather<T>( subrhs, num_rows_local, sol_vec, num_rows_local, 0, rhs_comm );
    }
    MPI_Barrier(this->global_comm);

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
    descA = (MKL_INT*)mkl_calloc( 9, sizeof(MKL_INT), FLOATING_PRECISION );
    descB = (MKL_INT*)mkl_calloc( 9, sizeof(MKL_INT), FLOATING_PRECISION );
    lddA = num_rows_local > 1 ? num_rows_local : 1;
    descinit_(descA, &num_rows, &num_rows, &num_block_size, &num_block_size, &izero, &izero, &ictxt, &lddA, &info);
    descinit_(descB, &num_rows, &num_cols_rhs, &num_block_size, &num_block_size_rhs, &izero, &izero, &ictxt, &lddA, &info);

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
    MKL_INT* num_rows_div   = NULL;
    MKL_INT* num_cols_div   = NULL;
    MKL_INT* rows_position  = NULL;
    MKL_INT* cols_position  = NULL;
    if (proc_rank == 0)
    {
        num_rows_div    = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );
        num_cols_div    = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );
        rows_position   = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );
        cols_position   = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );
    }

    MPI_Barrier(this->global_comm);
    MPI_Gather(&num_rows_local, 1, MPI_INT, num_rows_div, 1, MPI_INT, 0, this->global_comm);
    MPI_Gather(&num_cols_local, 1, MPI_INT, num_cols_div, 1, MPI_INT, 0, this->global_comm);
    MPI_Gather(&myrow, 1, MPI_INT, rows_position, 1, MPI_INT, 0, this->global_comm);
    MPI_Gather(&mycol, 1, MPI_INT, cols_position, 1, MPI_INT, 0, this->global_comm);
    MPI_Barrier(this->global_comm);

    // Allocate variables
    MKL_INT* start_cols = NULL;
    MKL_INT* start_rows = NULL;
    if (proc_rank == 0)
    {
        start_cols = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );
        start_rows = (MKL_INT*)mkl_calloc( num_procs, sizeof(MKL_INT), FLOATING_PRECISION );

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

    MPI_Barrier(this->global_comm);
    MPI_Scatter(start_rows, 1, MPI_INT, &start_row, 1, MPI_INT, 0, this->global_comm);
    MPI_Scatter(start_cols, 1, MPI_INT, &start_col, 1, MPI_INT, 0, this->global_comm);
    MPI_Barrier(this->global_comm);

    // Get end positions of the matrix
    this->end_row = this->start_row + this->num_rows_local;
    this->end_row = ( this->end_row > this->num_rows ) ? this->num_rows+1 : this->end_row;
    this->end_col = this->start_col + this->num_cols_local;
    this->end_col = ( this->end_col > this->num_rows ) ? this->num_rows+1 : this->end_col;
    this->end_row += 1;
    this->end_col += 1;

    // Get 0 reference position start and end intervals
    this->end_col_0     = this->end_col - 2;
    this->end_row_0     = this->end_row - 2;
    this->start_col_0   = this->start_col - 1;
    this->start_row_0   = this->start_row - 1;

    // Delete vectors
    if (proc_rank == 0)
    {
        mkl_free( num_rows_div  );
        mkl_free( num_cols_div  );
        mkl_free( rows_position );
        mkl_free( cols_position );
        mkl_free( start_cols    );
        mkl_free( start_rows    );
    }

    // Generate ipiv vector
    this->ipiv = (MKL_INT*)mkl_calloc( num_rows_local + num_block_size, sizeof(MKL_INT), FLOATING_PRECISION );

    this->GenerateRhsComm();
}


template <class T>
ScalapackSolver<T>::ScalapackSolver(
                                        MKL_INT     num_rows_inc, 
                                        MKL_INT     num_cols_rhs_inc,
                                        MKL_INT     num_procs_inc, 
                                        MKL_INT     proc_rank_inc,
                                        MKL_INT     proc_root_inc,
                                        MPI_Comm    global_comm_inc
                                    )
{
    // Assing variables
    this->global_comm   = global_comm_inc;
    this->proc_rank     = proc_rank_inc;
    this->proc_root     = proc_root_inc;
    this->num_rows      = num_rows_inc;
    this->num_procs     = num_procs_inc;
    this->num_cols_rhs  = num_cols_rhs_inc;

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

    num_block_size      = num_rows/num_procs_col + num_rows%num_procs_col;
    num_block_size      = num_block_size > 0 ? num_block_size: 1;
    num_block_size_rhs  = num_block_size > this->num_cols_rhs ? num_block_size : this->num_cols_rhs;

    // Initalize solver
    this->Initialize();

}


template<typename T>
ScalapackSolver<T>::~ScalapackSolver( )
{
    mkl_free( this->descA );
    mkl_free( this->descB );
    mkl_free( this->ipiv );
}


template <class T>
void ScalapackSolver<T>::Solve(T* subsysmat, T* subrhs)
{
    // Define local variables
    MKL_INT startrow    = 1;
    MKL_INT startcol    = 1;

    // Clear pivot vector
    for ( int i=0; i< ( num_rows_local + num_block_size ); i++ )
    {
        this->ipiv[i] = 0;
    }

    MPI_Barrier( this->global_comm );
    pgesv<T>(&num_rows, &num_cols_rhs, subsysmat, &startrow, &startcol, descA, this->ipiv, subrhs, &startrow, &startrow, descB, &info);
    MPI_Barrier( this->global_comm );

    if ( info != 0 )
    {
        std::cerr << "ERROR - Scalapack Solver" << std::endl;
        std::cout << "Scalapack solver finished abnormally with code: " << info << std::endl;
        throw std::runtime_error( "" );
    }
}

// Define custom types to have more handy ways to 
// work with ScalapackSolver implementations
typedef ScalapackSolver<cuscomplex> SclCmpx;
