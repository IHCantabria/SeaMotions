
// Include general usage libraries
#include <iostream>
#include <string>


int main( int argc, char* argv[] )
{
    // Convert input data
    int ref_mpi_build = std::stoi( argv[1] );

    // Check the macro definition to check if there is
    // a MPI build
    #ifdef _MPI_BUILD
    int mpi_build = 1;
    #else
    int mpi_build = 0;
    #endif

    // Check correspondence
    if ( ref_mpi_build != mpi_build )
    {
        std::cerr << std::endl;
        std::cerr << "Test test_mpi_build failed!" << std::endl;
        std::cerr << "There is no correspondence in between the mpi build ";
        std::cerr << "macro definition and the build signal value." << std::endl;
        std::runtime_error( "" );
    }

    return 0;
}