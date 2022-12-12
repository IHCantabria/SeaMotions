
// Include general usage libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/special_math.hpp"
#include "../../src/tools.hpp"


struct DataRef
{
    int num_points = 0;
    cusfloat* x = nullptr;
    cusfloat* y = nullptr;

    ~DataRef()
    {
        delete [] this->x;
        delete [] this->y;
    }

    void read_data(std::string file_path)
    {
        // Declare local variables
        std::string s0, s1, s2;

        // Open file unit
        std::ifstream infile(file_path);

        // Read number of points
        infile >> s0 >> s1 >> s2;
        infile >> this->num_points;

        // Allocate heap memory for the data
        this->x = new cusfloat [this->num_points];
        this->y = new cusfloat [this->num_points];

        // Read data from file
        infile >> s0 >> s1;
        for (int i=0; i<this->num_points; i++)
        {
            infile >> this->x[i] >> this->y[i];
        }

        // Close file unit
        infile.close();

    }

};


int main(int argc, char* argv[])
{
    // Check input arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path(argv[1]);

    // Read reference data for the test
    DataRef data_ref;
    data_ref.read_data(file_path);

    // Loop over the reference data values and check if the
    // output of the PSI function is correct according with the
    // machine precision defined during compilation
    cusfloat fi = 0.0;
    for (int i=0; i<data_ref.num_points; i++)
    {
        fi = psi_fun(static_cast<cusfloat>(data_ref.x[i]));
        if (std::abs(fi-data_ref.y[i]) > EPS_PRECISION)
        {
            std::cerr << "test_psi_function failed!" << std::endl;
        }
    }

    return 0;
}