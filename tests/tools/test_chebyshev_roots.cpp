
// Include general usage libraries
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// Include local modules
#include "../../src/chebyshev.hpp"
#include "../../src/config.hpp"
#include "../../src/tools.hpp"


struct RefData
{
    int num_polys = 0;
    cusfloat** poly_data = nullptr;
    int* poly_order = nullptr;

    ~RefData()
    {
        // Delete heap memory
        delete[] this->poly_order;
        for (int i=0; i<this->num_polys; i++)
        {
            delete[] this->poly_data[i];
        }
        delete[] this->poly_data;
    }

    int get_maximum_order(void)
    {
        int max_order = 0;
        for (int i=0; i<this->num_polys; i++)
        {
            if (this->poly_order[i] > max_order)
            {
                max_order = this->poly_order[i];
            }
        }

        return max_order;
    }

    void load_data(std::string file_path)
    {
        // Declare auxiliary variables
        std::string aux_str;

        // Open file unit
        std::ifstream infile(file_path);

        // Read number of polys
        infile >> aux_str;
        infile >> this->num_polys;
        
        // Allocate space for the polynomial orders and the
        // polynomials data
        this->poly_data = new cusfloat*[this->num_polys];
        this->poly_order = new int[this->num_polys];

        // Loop over roots polynomial definition to storage the data
        // in the class variables
        for (int i=0; i<this->num_polys; i++)
        {
            // Read poly order i
            infile >> aux_str >> this->poly_order[i];

            // Allocate space to storage the roots of the Chebyshev
            // polynomial of order i
            this->poly_data[i] = new cusfloat[this->poly_order[i]];

            // Read roots data from file
            for (int j=0; j<this->poly_order[i]; j++)
            {
                infile >> this->poly_data[i][j];
            }
        }

        // Close file unit
        infile.close();

        // // Show polynonmial file data
        // std::cout << "Num.Polys: " << this->num_polys << std::endl;
        // std::cout << "Poly orders: " << std::endl;
        // for (int i=0; i<this->num_polys; i++)
        // {
        //     std::cout << this->poly_order[i] << std::endl;
        // }

        // for (int i=0; i<this->num_polys; i++)
        // {
        //     std::cout << "Chebyt " << this->poly_order[i] << std::endl;
        //     for (int j=0; j<this->poly_order[i]; j++)
        //     {
        //         std::cout << j << " " << std::setprecision(15) << this->poly_data[i][j] << std::endl;
        //     }
        // }
    }
};


bool launch_test(RefData &ref_data)
{
    // Define pass flag
    bool pass = true;

    // Get maximum order of the polynomial roots
    // and allocate heap memeory for it
    int max_order = ref_data.get_maximum_order();
    cusfloat* poly_roots = new cusfloat[max_order];

    // Loop over polynomial orders checking the roots
    // values
    cusfloat diff = 0.0;
    for (int i=0; i<ref_data.num_polys; i++)
    {
        // Calculate roots points for the i degree of the 
        // Chebyshev polynomial
        chebyshev_poly_roots(ref_data.poly_order[i], poly_roots);

        // Check roots according to the working precision
        for (int j=0; j<ref_data.poly_order[i]; j++)
        {
            diff = poly_roots[j]-ref_data.poly_data[i][j];
            if (abs(diff) > EPS_PRECISION)
            {
                std::cerr << "Polynomial of order: " << ref_data.poly_order[i];
                std::cerr << " - Failed at root: " << ref_data.poly_data[i][j];
                std::cerr << " - Err: " << diff << std::endl;
                std::cerr << "test_chebyshe _roots failed!" << std::endl;
                return false;
            }
        }
    }

    // Deallocate heap memory
    delete[] poly_roots;

    return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path_chebyshev(argv[1]);

    // Load data from file
    RefData ref_data;
    ref_data.load_data(file_path_chebyshev);

    // Launch test
    int pass = launch_test(ref_data);
    if (!pass)
    {
        return 1;
    }

    return 0;
}