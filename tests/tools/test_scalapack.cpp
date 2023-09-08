
// Include general usage libraries
#include <fstream>
#include "mpi.h"

// Include general usage scientific libraries
#include "../../src/containers/mpi_config.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/scalapack_solver.hpp"


template<typename T>
struct RefData
{
private:
    // Declare class attributes
    bool _is_heap = false;

    // Declare class methods
    void _load_data( std::string fipath );

public:
    // Declare class attributes
    int         rows_np = 0;
    T*          sol     = nullptr;
    T*          sysmat  = nullptr;
    T*          sysrhs  = nullptr;

    // Define Class constructors and destructor
    RefData( std::string fipath )
    {
        this->_load_data( fipath );
    }

    ~RefData( void )
    {
        if ( this->_is_heap )
        {
            mkl_free( this->sol );
            mkl_free( this->sysmat );
            mkl_free( this->sysrhs );
        }
    }

};

template<>
void RefData<cusfloat>::_load_data( std::string fipath )
{
    // Define aux variables
    std::string aux_str;

    // Open file unit
    std::ifstream infile( fipath );

    // Read number of rows
    infile >> aux_str >> this->rows_np;

    // Allocate space for the matrix containers
    this->sysmat    = generate_empty_vector<cusfloat>( this->rows_np * this->rows_np );
    this->sysrhs    = generate_empty_vector<cusfloat>( this->rows_np );
    this->sol       = generate_empty_vector<cusfloat>( this->rows_np );

    // Read system matrix
    infile >> aux_str;

    for ( int i=0; i<(this->rows_np*this->rows_np); i++ )
    {
        infile >> this->sysmat[i];
    }

    // Read RHS
    infile >> aux_str;

    for ( int i=0; i<this->rows_np; i++ )
    {
        infile >> this->sysrhs[i];
    }

    // Read solution vector
    infile >> aux_str;

    for ( int i=0; i<this->rows_np; i++ )
    {
        infile >> this->sol[i];
    }

    // Close file unit
    infile.close( );

}


template<>
void RefData<cuscomplex>::_load_data( std::string fipath )
{
    // Define aux variables
    std::string aux_str;
    cusfloat    a=0.0, b=0.0;

    // Open file unit
    std::ifstream infile( fipath );

    // Read number of rows
    infile >> aux_str >> this->rows_np;

    // Allocate space for the matrix containers
    this->sysmat    = generate_empty_vector<cuscomplex>( this->rows_np * this->rows_np );
    this->sysrhs    = generate_empty_vector<cuscomplex>( this->rows_np );
    this->sol       = generate_empty_vector<cuscomplex>( this->rows_np );

    // Read system matrix
    infile >> aux_str;

    for ( int i=0; i<(this->rows_np*this->rows_np); i++ )
    {
        // Read real and imaginary parts from file
        infile >> a >> b;

        // Create complex number
        this->sysmat[i] = cuscomplex( a, b );
    }

    // Read RHS
    infile >> aux_str;

    for ( int i=0; i<this->rows_np; i++ )
    {
        // Read real and imaginary parts from file
        infile >> a >> b;

        // Create complex number
        this->sysrhs[i] = cuscomplex( a, b );
    }

    // Read solution vector
    infile >> aux_str;

    for ( int i=0; i<this->rows_np; i++ )
    {
        // Read real and imaginary parts
        infile >> a >> b;

        // Create complex number
        this->sol[i] = cuscomplex( a, b );
    }

    // Close file unit
    infile.close( );

}


template<typename T>
void fill_sys_chunks( 
                        T*                  full_mat,
                        T*                  sysmat,
                        T*                  full_rhs,
                        T*                  sysrhs,
                        ScalapackSolver<T>& scl
                    )
{
    // Loop over rows
    int col_global_idx  = 0;
    int count           = 0;
    for ( int j=0; j<scl.num_cols_local; j++ )
    {
        for ( int i=0; i<scl.num_rows_local; i++ )
        {
            col_global_idx  = ( scl.start_row -1 + i )*scl.num_rows + scl.start_col - 1 + j;
            sysmat[count] = full_mat[col_global_idx];
            count++;
        }
    }

    for ( int i=0; i<scl.num_rows_local; i++ )
    {
        sysrhs[i] = full_rhs[scl.start_row-1+i];
    }

}

template<typename T>
void launch_test_simple( std::string res_fipath )
{
    // Initializa MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &procs_total
                );

    // Get current processor rank
    int proc_rank  = 0;
    MPI_Comm_rank( 
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Create instance of ScalapackSolver object
    int N = 2;

    ScalapackSolver<T> scl( N, procs_total, proc_rank );

    // Allocate space for the test matrix
    T* sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
    T* sysrhs = generate_empty_vector<T>( scl.num_rows_local );

    // Fill submatrix
    sysmat[0] = 1.0;
    sysmat[1] = 2.0;
    sysmat[2] = 5.0;
    sysmat[3] = 3.0;

    sysrhs[0] = 2.0;
    sysrhs[1] = 3.0;

    // Solve 
    scl.Solve(sysmat, sysrhs);


    // Show solution
    std::cout << "SOLUTION: " << std::endl;
    print_vector<cusfloat>( N, sysrhs, 1, 6 );

    // Deallocate heap memory in each process
    mkl_free( sysmat );
    mkl_free( sysrhs );

    // Close MPI environment
    MPI_Finalize( );
}


template<typename T>
void launch_test_simple_2( std::string res_fipath )
{
    // Initializa MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &procs_total
                );

    // Get current processor rank
    int proc_rank  = 0;
    MPI_Comm_rank( 
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Create instance of ScalapackSolver object
    int N = 4;

    ScalapackSolver<T> scl( N, procs_total, proc_rank );

    // Allocate space for the test matrix
    T* sysmat = nullptr;
    T* sysrhs = nullptr;
    if ( proc_rank == 0 )
    {
        std::cout << "Generating heap memory process 0" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 0 -> Done" << std::endl;
    }
    else if ( proc_rank == 1 )
    {
        std::cout << "Generating heap memory process 1" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 1 -> Done" << std::endl;
    }

    // Fill submatrix
    if ( proc_rank == 0 )
    {
        sysmat[0] = 1.0;
        sysmat[1] = 4.0;
        sysmat[2] = 8.0;
        sysmat[3] = 5.0;
        sysmat[4] = 3.0;
        sysmat[5] = 1.0;
        sysmat[6] = 5.0;
        sysmat[7] = 6.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
        sysrhs[2] = 3.0;
        sysrhs[3] = 6.0;
    }
    else if ( proc_rank == 1 )
    {
        sysmat[0] = 4.0;
        sysmat[1] = 2.0;
        sysmat[2] = 4.0;
        sysmat[3] = 2.0;
        sysmat[4] = 9.0;
        sysmat[5] = 3.0;
        sysmat[6] = 9.0;
        sysmat[7] = 3.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
        sysrhs[2] = 3.0;
        sysrhs[3] = 6.0;
    }

    // Solve 
    // std::cout << "here 2 - Processor: " << proc_rank << std::endl;
    scl.Solve(sysmat, sysrhs);
    // std::cout << "done - Processor: " << proc_rank << std::endl;


    // Show total solution
    cusfloat* total_rhs = nullptr;
    if ( proc_rank == 0 )
    {
        total_rhs = generate_empty_vector<cusfloat>( N );
    }

    scl.GetGlobalRhs( sysrhs, total_rhs );

    if ( proc_rank == 0 )
    {
        std::cout << "Solution: " << std::endl;
        print_vector( N, total_rhs, 1, 6 );
    }

    // Deallocate heap memory in each process
    mkl_free( sysmat );
    mkl_free( sysrhs );
    if ( proc_rank == 0 )
    {
        mkl_free( total_rhs );
    }


    // Close MPI environment
    MPI_Finalize( );
}


template<typename T>
void launch_test_simple_3( std::string res_fipath )
{
    // Initializa MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &procs_total
                );

    // Get current processor rank
    int proc_rank  = 0;
    MPI_Comm_rank( 
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Create instance of ScalapackSolver object
    int N = 4;

    ScalapackSolver<T> scl( N, procs_total, proc_rank );

    // Allocate space for the test matrix
    T* sysmat = nullptr;
    T* sysrhs = nullptr;
    if ( proc_rank == 0 )
    {
        std::cout << "Generating heap memory process 0" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 0 -> Done" << std::endl;
    }
    else if ( proc_rank == 1 )
    {
        std::cout << "Generating heap memory process 1" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 1 -> Done" << std::endl;
    }
    else if ( proc_rank == 2 )
    {
        std::cout << "Generating heap memory process 2" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 1 -> Done" << std::endl;
    }

    // Fill submatrix
    if ( proc_rank == 0 )
    {
        sysmat[0] = 1.0;
        sysmat[1] = 4.0;
        sysmat[2] = 8.0;
        sysmat[3] = 5.0;
        sysmat[4] = 3.0;
        sysmat[5] = 1.0;
        sysmat[6] = 5.0;
        sysmat[7] = 6.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
        sysrhs[2] = 3.0;
        sysrhs[3] = 6.0;
    }
    else if ( proc_rank == 1 )
    {
        sysmat[0] = 4.0;
        sysmat[1] = 2.0;
        sysmat[2] = 4.0;
        sysmat[3] = 2.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
        sysrhs[2] = 3.0;
        sysrhs[3] = 6.0;
    }
    else if ( proc_rank == 2 )
    {
        sysmat[0] = 9.0;
        sysmat[1] = 3.0;
        sysmat[2] = 9.0;
        sysmat[3] = 3.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
        sysrhs[2] = 3.0;
        sysrhs[3] = 6.0;
    }

    // Solve 
    // std::cout << "here 2 - Processor: " << proc_rank << std::endl;
    scl.Solve(sysmat, sysrhs);
    // std::cout << "done - Processor: " << proc_rank << std::endl;


    // Show solution
    std::cout << "SOLUTION: " << std::endl;
    for ( int i=0; i<scl.num_rows_local; i++ )
    {
        std::cout << "ProcRank: " << proc_rank << " - V[" << i << "]: " << sysrhs[i] << std::endl;
    }
    // print_vector<cusfloat>( 2, sysrhs, 1, 6 );

    cusfloat* total_rhs = nullptr;
    if ( proc_rank == 0 )
    {
        total_rhs = generate_empty_vector<cusfloat>( N );
    }

    scl.GetGlobalRhs( sysrhs, total_rhs );

    if ( proc_rank == 0 )
    {
        std::cout << "Solution: " << std::endl;
        print_vector( N, total_rhs, 1, 6 );
    }

    // Deallocate heap memory in each process
    mkl_free( sysmat );
    mkl_free( sysrhs );
    if ( proc_rank == 0 )
    {
        mkl_free( total_rhs );
    }


    // Close MPI environment
    MPI_Finalize( );
}


template<typename T>
void launch_test_simple_4( std::string res_fipath )
{
    // Initializa MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    int procs_total = 0;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &procs_total
                );

    // Get current processor rank
    int proc_rank  = 0;
    MPI_Comm_rank( 
                    MPI_COMM_WORLD,
                    &proc_rank
                );

    // Create instance of ScalapackSolver object
    int N = 4;

    // std::cout << "here 1 - Processor: " << proc_rank << std::endl;
    ScalapackSolver<T> scl( N, procs_total, proc_rank );
    // std::cout << "done - Processor: " << proc_rank << std::endl;

    // Allocate space for the test matrix
    T* sysmat = nullptr;
    T* sysrhs = nullptr;
    if ( proc_rank == 0 )
    {
        std::cout << "Generating heap memory process 0" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 0 -> Done" << std::endl;
    }
    else if ( proc_rank == 1 )
    {
        std::cout << "Generating heap memory process 1" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 1 -> Done" << std::endl;
    }
    else if ( proc_rank == 2 )
    {
        std::cout << "Generating heap memory process 2" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 2 -> Done" << std::endl;
    }
    else if ( proc_rank == 3 )
    {
        std::cout << "Generating heap memory process 3" << std::endl;
        std::cout << "Proc: " << proc_rank << " - NumRowsLocal: " << scl.num_cols_local << " - NumColsLocal: " << scl.num_cols_local << std::endl;
        sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
        sysrhs = generate_empty_vector<T>( scl.num_rows_local );
        std::cout << "Generating heap memory process 3 -> Done" << std::endl;
    }

    // std::cout << "Num rows local: " << scl.num_rows_local << std::endl;
    // std::cout << "Num cols local: " << scl.num_cols_local << std::endl;

    // // Fill submatrix
    // sysmat[0] = 1.0;
    // sysmat[1] = 4.0;
    // sysmat[2] = 8.0;
    // sysmat[3] = 5.0;

    // sysmat[4] = 3.0;
    // sysmat[5] = 1.0;
    // sysmat[6] = 5.0;
    // sysmat[7] = 6.0;

    // sysmat[8] = 4.0;
    // sysmat[9] = 2.0;
    // sysmat[10] = 4.0;
    // sysmat[11] = 2.0;

    // sysmat[12] = 9.0;
    // sysmat[13] = 3.0;
    // sysmat[14] = 9.0;
    // sysmat[15] = 3.0;

    // sysrhs[0] = 8.0;
    // sysrhs[1] = 5.0;
    // sysrhs[2] = 3.0;
    // sysrhs[3] = 6.0;

    if ( proc_rank == 0 )
    {
        sysmat[0] = 1.0;
        sysmat[1] = 4.0;
        sysmat[2] = 3.0;
        sysmat[3] = 1.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
    }
    else if ( proc_rank == 1 )
    {
        sysmat[0] = 4.0;
        sysmat[1] = 2.0;
        sysmat[2] = 9.0;
        sysmat[3] = 3.0;

        sysrhs[0] = 8.0;
        sysrhs[1] = 5.0;
    }
    else if ( proc_rank == 2 )
    {
        sysmat[0] = 8.0;
        sysmat[1] = 5.0;
        sysmat[2] = 5.0;
        sysmat[3] = 6.0;

        sysrhs[0] = 3.0;
        sysrhs[1] = 6.0;
    }
    else if ( proc_rank == 3 )
    {
        sysmat[0] = 4.0;
        sysmat[1] = 2.0;
        sysmat[2] = 9.0;
        sysmat[3] = 3.0;

        sysrhs[0] = 3.0;
        sysrhs[1] = 6.0;
    }

    // Solve 
    // std::cout << "here 2 - Processor: " << proc_rank << std::endl;
    scl.Solve(sysmat, sysrhs);
    // std::cout << "done - Processor: " << proc_rank << std::endl;


    // Show total solution
    cusfloat* total_rhs = nullptr;
    if ( proc_rank == 0 )
    {
        total_rhs = generate_empty_vector<cusfloat>( N );
    }

    scl.GetGlobalRhs( sysrhs, total_rhs );

    if ( proc_rank == 0 )
    {
        std::cout << "Solution: " << std::endl;
        print_vector( N, total_rhs, 1, 6 );
    }

    // Deallocate heap memory in each process
    mkl_free( sysmat );
    mkl_free( sysrhs );
    if ( proc_rank == 0 )
    {
        mkl_free( total_rhs );
    }

    // Close MPI environment
    MPI_Finalize( );
}


template<typename T>
void launch_test( 
                    std::string res_fipath,
                    MpiConfig&  mpi_config
                )
{
    // Read input data
    RefData<T> ref_data( res_fipath );

    // Create instance of ScalapackSolver object
    ScalapackSolver<T> scl( ref_data.rows_np, 1, mpi_config.procs_total, mpi_config.proc_rank );

    // Allocate space for the test matrix
    T* sysmat = generate_empty_vector<T>( scl.num_rows_local*scl.num_cols_local );
    T* sysrhs = generate_empty_vector<T>( scl.num_rows_local );

    // Fill submatrix
    fill_sys_chunks(
                        ref_data.sysmat,
                        sysmat,
                        ref_data.sysrhs,
                        sysrhs,
                        scl
                    );

    // Solve 
    scl.Solve(sysmat, sysrhs);

    // Get total solution
    T* global_sol = nullptr;
    if ( mpi_config.proc_rank == 0 )
    {
        // Allocate memory
        global_sol = generate_empty_vector<T>( scl.num_rows );
    }
        // Get global solution
    scl.GetGlobalRhs( sysrhs, global_sol );

    if ( mpi_config.proc_rank == 0 )
    {
        // Check the solution
        cusfloat eps = 1e-6;
        if ( !assert_vector_equality( 
                                        scl.num_rows, 
                                        global_sol, 
                                        ref_data.sol, 
                                        eps 
                                    ) 
            )
        {
            std::cout << std::endl;
            std::cout << "test test_scalapack failed!" << std::endl;
            std::cout << "Local solution for processor: " << mpi_config.proc_rank;
            std::cout << " is not equal to the reference data." << std::endl;
            std::cout << std::endl;
            for ( int i=0; i<scl.num_rows; i++ )
            {
                std::cout << "i[" << i << "] - global_sol: " << global_sol[i] << " - sol: " << ref_data.sol[i] << std::endl;
            }
            throw std::runtime_error( "" );
        }

    }

    // Deallocate heap memory in each process
    mkl_free( sysmat );
    mkl_free( sysrhs );
    if ( mpi_config.proc_rank == 0 )
    {
        mkl_free( global_sol );
    }
   
}


int main( int argc, char* argv[] )
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 2))
    {
        return 1;
    }

    std::string cusfloat_fipath(argv[1]);
    std::string cuscomplex_fipath(argv[2]);

    // Initialize MPI environment
    MPI_Init( NULL, NULL );

    // Get total number of processors
    MpiConfig mpi_config;
    MPI_Comm_size(
                    MPI_COMM_WORLD,
                    &(mpi_config.procs_total)
                );

    // Get current processor rank
    MPI_Comm_rank( 
                    MPI_COMM_WORLD,
                    &(mpi_config.proc_rank)
                );

    // Launch test for real numbers
    launch_test<cusfloat>( 
                                cusfloat_fipath,
                                mpi_config
                            );

    // Launch test for complex numbers
    launch_test<cuscomplex>( 
                                cuscomplex_fipath,
                                mpi_config
                            );

    // Close MPI environment
    MPI_Finalize( );

    return 0;
}