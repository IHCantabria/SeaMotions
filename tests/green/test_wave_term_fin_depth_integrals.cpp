
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
    bool is_loaded_1d = false;
    bool is_loaded_3d = false;
    int num_points = 0;

    ~RefData(void)
    {
        if (this->is_loaded_1d)
        {
            delete [] this->g_int;
            delete [] this->h;
        }
        else if (this->is_loaded_3d)
        {
            delete [] this->a;
            delete [] this->b;
            delete [] this->g_int;
            delete [] this->h;
        }
    }

    void load_data_1d(std::string file_path)
    {
        // Create auxiliar strings
        std::string s0, s1, s2, s3;

        // Open file unit
        std::ifstream infile(file_path);

        // Read number of points
        infile >> s0 >> this->num_points;

        // Allocate heap memory
        this->g_int = new cusfloat [this->num_points];
        this->h = new cusfloat [this->num_points];

        // Read header
        infile >> s0 >> s1;

        // Loop over points to read data
        for (int i=0; i<this->num_points; i++)
        {
            infile >> this->h[i] >> this->g_int[i];
        }

        // Activate flag to show that the data is
        // loaded
        this->is_loaded_1d = true;    

        // Close file unit
        infile.close();
    }

    void load_data_3d(std::string file_path)
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
        this->is_loaded_3d = true;    

        // Close file unit
        infile.close();
    }

    void print_data(void)
    {
        if(this->is_loaded_1d)
        {
            std::cout << "Number of points: " << this->num_points << std::endl;
            std::cout << "   H    " << "   G_int  " << std::endl;
            for (int i=0; i<this->num_points; i++)
            {
                std::cout << this->h[i] << " " << this->g_int[i] << std::endl;
            }
        }
        else if(this->is_loaded_3d)
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


bool launch_test_1d(P3* lobj, std::string file_path, std::string test_name)
{
    // Read reference data
    RefData ref_data;
    ref_data.load_data_1d(file_path);

    // Loop over data checking the subroutines
    bool flag = true;
    cusfloat diff = 0.0;
    for (int i=0; i<ref_data.num_points; i++)
    {
        // Check using direct method
        lobj->calculate_h_1D(ref_data.h[i]);
        diff = lobj->int_1d-ref_data.g_int[i];
        if (abs(diff)>EPS_INTEGRALS)
        {
            std::cout << "Test " << test_name << " - direct method: has ";
            std::cout << "an error of: " << diff;
            std::cout << " which is above the requested limits." << std::endl;
            flag = false;
            break;
        }

    }

    return flag;
}


bool launch_test_3d(P3* lobj, std::string file_path, std::string test_name)
{
    // Read reference data
    RefData ref_data;
    ref_data.load_data_3d(file_path);

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
    if (!check_num_cmd_args(argc, 14))
    {
        return 1;
    }

    std::string file_path_l1(argv[1]);
    std::string file_path_l1_da(argv[2]);
    std::string file_path_l1_db(argv[3]);
    std::string file_path_l2(argv[4]);
    std::string file_path_l3(argv[5]);
    std::string file_path_l3_da(argv[6]);
    std::string file_path_l3_db(argv[7]);
    std::string file_path_m1(argv[8]);
    std::string file_path_m1_da(argv[9]);
    std::string file_path_m1_db(argv[10]);
    std::string file_path_m2(argv[11]);
    std::string file_path_m3(argv[12]);
    std::string file_path_m3_da(argv[13]);
    std::string file_path_m3_db(argv[14]);

    // Launch test L1
    P3* l1 = new P3();
    set_data_l1(l1);
    launch_test_3d(l1, file_path_l1, "L1");

    // Launch test L1_dA
    P3* l1_da = new P3();
    set_data_l1_da(l1_da);
    launch_test_3d(l1_da, file_path_l1_da, "L1_dA");

    // Launch test L1_dB
    P3* l1_db = new P3();
    set_data_l1_db(l1_db);
    launch_test_3d(l1_db, file_path_l1_db, "L1_dB");

    // Launch test L2
    P3* l2 = new P3();
    set_data_l2(l2);
    launch_test_1d(l2, file_path_l2, "L2");

    // Launch test L3
    P3* l3 = new P3();
    set_data_l3(l3);
    launch_test_3d(l3, file_path_l3, "L3");

    // Launch test L3_dA
    P3* l3_da = new P3();
    set_data_l3_da(l3_da);
    launch_test_3d(l3_da, file_path_l3_da, "L3_dA");
    
    // Launch test L3_dB
    P3* l3_db = new P3();
    set_data_l3_db(l3_db);
    launch_test_3d(l3_db, file_path_l3_db, "L3_dB");

    // Launch test M1
    P3* m1 = new P3();
    set_data_m1(m1);
    launch_test_3d(m1, file_path_m1, "M1");

    // Launch test M1_dA
    P3* m1_da = new P3();
    set_data_m1_da(m1_da);
    launch_test_3d(m1_da, file_path_m1_da, "M1_dA");

    // Launch test M1_dB
    P3* m1_db = new P3();
    set_data_m1_db(m1_db);
    launch_test_3d(m1_db, file_path_m1_db, "M1_dB");

    // Launch test M2
    P3* m2 = new P3();
    set_data_m2(m2);
    launch_test_1d(m2, file_path_m2, "M2");

    // Launch test M3
    P3* m3 = new P3();
    set_data_m3(m3);
    launch_test_3d(m3, file_path_m3, "M3");

    // Launch test M3_dA
    P3* m3_da = new P3();
    set_data_m3_da(m3_da);
    launch_test_3d(m3_da, file_path_m3_da, "M3_dA");

    // Launch test M3_dB
    P3* m3_db = new P3();
    set_data_m3_db(m3_db);
    launch_test_3d(m3_db, file_path_m3_db, "M3_dB");

    // Delete heap memory allocation
    delete l1;
    delete l1_da;
    delete l1_db;
    delete l2;
    delete l3;
    delete l3_da;
    delete l3_db;
    delete m1;
    delete m1_da;
    delete m1_db;
    delete m2;
    delete m3;
    delete m3_da;
    delete m3_db;

    return 0;
}