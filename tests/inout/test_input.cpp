
// Include general usage libraries
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// Include local modules
#include "../../src/inout/input.hpp"
#include "../../src/inout/reader.hpp"
#include "../../src/tools.hpp"


int main( int argc, char* argv[] )
{
    // Read command line arguments
    if ( !check_num_cmd_args(argc, 1) )
    {
        return 1;
    }

    std::string case_fopath(argv[1]);

    // Read case inputs
    Input* input = read_input_files( case_fopath );

    // Print out inputs
    input->print( );

    return 0;
}