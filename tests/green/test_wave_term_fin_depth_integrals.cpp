
// Include general usage libraries
#include <fstream>
#include <iostream>
#include <string>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating_fin_depth_cheby.hpp"
#include "../../src/tools.hpp"


cusfloat EPS_INTEGRALS = 1e-6;


struct RefData
{
    cusfloat* a;
    cusfloat* b;
    cusfloat* g_int;
    cusfloat* h;
    bool is_loaded = false;
    int num_points = 0;

    ~RefData(void)
    {
        if (this->is_loaded)
        {
            delete [] this->a;
            delete [] this->b;
            delete [] this->g_int;
            delete [] this->h;
        }
    }

    void load_data(std::string file_path)
    {
        // Create auxiliar strings
        std::string s0, s1, s2, s3;

        // Open file unit
        std::ifstream infile(file_path);

        // Read number of points
        infile >> s0 >> this->num_points;

        // Allocate heap memory
        this->a = new cusfloat [this->num_points];
        this->b = new cusfloat [this->num_points];
        this->g_int = new cusfloat [this->num_points];
        this->h = new cusfloat [this->num_points];

        // Read header
        infile >> s0 >> s1 >> s2 >> s3;

        // Loop over points to read data
        for (int i=0; i<this->num_points; i++)
        {
            infile >> this->a[i] >> this->b[i] >> this->h[i] >> this->g_int[i];
        }

        // Activate flag to show that the data is
        // loaded
        this->is_loaded = true;    

        // Close file unit
        infile.close();
    }

    void print_data(void)
    {
        if(this->is_loaded)
        {
            std::cout << "Number of points: " << this->num_points << std::endl;
            std::cout << "   A   " << "   B   " << "   H    " << "   G_int  " << std::endl;
            for (int i=0; i<this->num_points; i++)
            {
                std::cout << this->a[i] << " " << this->b[i] << " ";
                std::cout << this->h[i] << " " << this->g_int[i] << std::endl;
            }
        }
    }
};


bool launch_test_3d(P3* lobj, std::string file_path, std::string test_name)
{
    // Read reference data
    RefData ref_data;
    ref_data.load_data(file_path);

    // Loop over data checking the subroutines
    bool flag = true;
    cusfloat diff = 0.0;
    cusfloat sol_i = 0.0;
    for (int i=0; i<ref_data.num_points; i++)
    {
        // Check using direct method
        sol_i = lobj->get_value_abh(
                                    ref_data.a[i],
                                    ref_data.b[i],
                                    ref_data.h[i]
                                    );
        diff = sol_i-ref_data.g_int[i];
        if (abs(diff)>EPS_INTEGRALS)
        {
            std::cout << "Test " << test_name << " - direct method: has ";
            std::cout << "an error of: " << diff;
            std::cout << " which is above the requested limits." << std::endl;
            flag = false;
            break;
        }

        // Check using simple folding method
        lobj->fold_h(ref_data.h[i]);
        // lobj->fold_b(ref_data.b[i]);
        // sol_i = lobj->get_value_a(ref_data.a[i]);
        sol_i = lobj->get_value_ab(ref_data.a[i], ref_data.b[i]);
        diff = sol_i-ref_data.g_int[i];
        if (abs(diff)>EPS_INTEGRALS)
        {
            std::cout << "Test " << test_name << " - simple fold method: has ";
            std::cout << "an error of: " << diff;
            std::cout << " which is above the requested limits." << std::endl;
            flag = false;
            break;
        }

        // Check using simple folding method
        lobj->fold_h(ref_data.h[i]);
        lobj->fold_b(ref_data.b[i]);
        sol_i = lobj->get_value_a(ref_data.a[i]);
        diff = sol_i-ref_data.g_int[i];
        if (abs(diff)>EPS_INTEGRALS)
        {
            std::cout << "Test " << test_name << " - double fold method: has ";
            std::cout << "an error of: " << diff;
            std::cout << " which is above the requested limits." << std::endl;
            flag = false;
            break;
        }
    }

    return flag;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path_l1(argv[1]);

    // Launch test L1
    P3* l1 = new P3();
    set_data_l1(l1);
    launch_test_3d(l1, file_path_l1, "L1");
    

    // Delete heap memory allocation
    delete l1;

    return 0;
}