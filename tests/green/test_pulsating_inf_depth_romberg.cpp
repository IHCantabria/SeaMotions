
// Include general usage libraries
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

// Include genera usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating.hpp"
#include "../../src/math_tools.hpp"


struct DataTable
{
    cusfloat* int_value = nullptr;
    int num_points_x = 0;
    int num_points_y = 0;
    cusfloat* x = nullptr;
    cusfloat* y = nullptr;

    ~DataTable()
    {
        delete [] this->int_value;
        delete [] this->x;
        delete [] this->y;
    }

    void read_data(std::string file_path)
    {
        // Open file unit
        std::ifstream infile(file_path);

        // Declare local dummy variables
        std::string str_aux;

        // Read number of X points
        infile >> str_aux;
        infile >> this->num_points_x;
        infile >> str_aux;
        this->x = new cusfloat [this->num_points_x];
        for (int i=0; i<this->num_points_x; i++)
        {
            infile >> this->x[i];
        }

        // Read number of Y points
        infile >> str_aux;
        infile >> this->num_points_y;
        infile >> str_aux;
        this->y = new cusfloat [this->num_points_y];
        for (int i=0; i<this->num_points_y; i++)
        {
            infile >> this->y[i];
        }

        // Read number of integral points
        infile >> str_aux;
        int num_int_points = this->num_points_x*this->num_points_y;
        this->int_value = new cusfloat [num_int_points];
        for (int i=0; i<num_int_points; i++)
        {
            infile >> this->int_value[i];
        }

        // Close file unit
        infile.close();
    }
};


template<typename Functor>
void print_t(Functor f)
{
    for (int i=0; i<10; i++)
    {
        std::cout << "t: " << i << " - f(i): " << f(static_cast<cusfloat>(i)) << std::endl;
    }
}


bool launch_test(std::string file_path, std::function<cusfloat(cusfloat, cusfloat, cusfloat)> f)
{
    // Read data from file
    DataTable data_table;
    data_table.read_data(file_path);

    // Loop over reference values to check the integral representation
    cusfloat abs_err = 0.0;
    int count = 0;
    cusfloat fi = 0.0;
    bool pass = true;
    cusfloat X = 0.0, Y = 0.0;
    for (int j=0; j<data_table.num_points_y; j++)
    {
        for (int i=0; i<data_table.num_points_x; i++)
        {
            // Calculate numerical integral
            X = data_table.x[i];
            Y = data_table.y[j];
            fi = romberg_quadrature(
                [X, Y, f](cusfloat t)->cusfloat {return f(X, Y, t);},
                0,
                Y,
                1e-12
            );

            // Check error
            abs_err = fi - data_table.int_value[count];
            if (std::abs(abs_err) > 1e-12)
            {
                std::cerr << "X: " << X << "- Y: " << Y << " - Int(num): " << std::setprecision(15) << fi;
                std::cerr << " - Int(ref): " << std::setprecision(15) << data_table.int_value[count];
                std::cerr << " - Diff: " << std::setprecision(15) << abs_err << " - Line: " << count+10 << std::endl;
                pass = false;
                return pass;
            }

            // Advance counter
            count++;
        }
    }

    return pass;
}


bool launch_test(std::string file_path, std::function<cusfloat(cusfloat, cusfloat)> f)
{
    // Read data from file
    DataTable data_table;
    data_table.read_data(file_path);

    // Loop over reference values to check the integral representation
    cusfloat abs_err = 0.0;
    int count = 0;
    cusfloat fi = 0.0;
    bool pass = true;
    cusfloat X = 0.0, Y = 0.0;
    for (int j=0; j<data_table.num_points_y; j++)
    {
        for (int i=0; i<data_table.num_points_x; i++)
        {
            // Calculate numerical integral
            X = data_table.x[i];
            Y = data_table.y[j];
            fi = f(X, Y);

            // Check error
            abs_err = fi - data_table.int_value[count];
            if (std::abs(abs_err) > 1e-6)
            {
                std::cerr << "X: " << X << "- Y: " << Y << " - Int(num): " << std::setprecision(15) << fi;
                std::cerr << " - Int(ref): " << std::setprecision(15) << data_table.int_value[count];
                std::cerr << " - Diff: " << std::setprecision(15) << abs_err << " - Line: " << count+10 << std::endl;
                pass = false;
                return pass;
            }

            // Advance counter
            count++;
        }
    }

    return pass;
}


int main(int argc, char* argv[])
{
    // Check input arguments
    if (!check_num_cmd_args(argc, 4))
    {
        return 1;
    }

    std::string file_path_table_1(argv[1]);
    std::string file_path_table_2(argv[2]);
    std::string file_path_table_3(argv[3]);
    std::string file_path_table_4(argv[4]);

    // Define local variables
    bool pass = false;

    // Launch test for integral table 1 
    pass = launch_test(file_path_table_1, expint_inf_depth_num);
    if (!pass)
    {
        std::cerr << "test_pulsating_inf_depth_romberg/sub_test_table_1 failed!" << std::endl;
        return 1;
    }

    // Launch test for integral table 2
    pass = launch_test(file_path_table_2, expint_inf_depth_num_dx);
    if (!pass)
    {
        std::cerr << "test_pulsating_inf_depth_romberg/sub_test_table_2 failed!" << std::endl;
        return 1;
    }

    // Launch test for integral table 3
    pass = launch_test(file_path_table_3, wave_term_inf_depth_num);
    if (!pass)
    {
        std::cerr << "test_pulsating_inf_depth_romberg/sub_test_table_3 failed!" << std::endl;
        return 1;
    }

    // Launch test for integral table 4
    pass = launch_test(file_path_table_4, wave_term_inf_depth_num_dx);
    if (!pass)
    {
        std::cerr << "test_pulsating_inf_depth_romberg/sub_test_table_4 failed!" << std::endl;
        return 1;
    }

    return 0;
}