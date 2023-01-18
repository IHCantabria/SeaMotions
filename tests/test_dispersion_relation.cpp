
// Include general usage libraries
#include <fstream>
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../src/config.hpp"
#include "../src/math/math_tools.hpp"
#include "../src/tools.hpp"
#include "../src/waves.hpp"


// Define precision for the calculation of the wave number
cusfloat WAVENUM_TOL = 1e-6;


struct RefData
{
    cusfloat grav_acc = 0.0;
    int num_depth = 0;
    int num_period = 0;
    cusfloat* water_depth = nullptr;
    cusfloat* wave_number = nullptr;
    cusfloat* wave_period = nullptr;

    ~RefData()
    {
        delete [] this->water_depth;
        delete [] this->wave_number;
        delete [] this->wave_period;
    }

    void load_data(std::string file_path)
    {
        // Declare auxiliar string to read header info
        std::string str_aux;

        // Open file unit
        std::ifstream infile(file_path);

        // Read gravity data
        infile >> str_aux >> this->grav_acc;

        // Read Water depth header
        infile >> str_aux >> this->num_depth;
        this->water_depth = new cusfloat[this->num_depth];
        for (int i=0; i<this->num_depth; i++)
        {
            infile >> this->water_depth[i];
        }

        // Read Wave Period header
        infile >> str_aux >> this->num_period;
        this->wave_period = new cusfloat [this->num_period];
        for (int i=0; i<this->num_period; i++)                                                                    
        {
            infile >> this->wave_period[i];
        }

        // Read Wave numbers
        infile >> str_aux;
        this->wave_number = new cusfloat[this->num_period*this->num_depth];
        for (int i=0; i<this->num_depth; i++)
        {
            for (int j=0; j<this->num_period; j++)
            {
                infile >> this->wave_number[i*this->num_period+j];
            }
        }

        // Close file unit
        infile.close();
    }
};


bool test_dispersion_real(std::string file_path)
{
    // Declare test variable
    bool pass = true;
    
    // Read reference data
    RefData ref_data;
    ref_data.load_data(file_path);

    // Loop over reference data to check
    // the calculated results are in agreement
    cusfloat diff = 0.0;
    cusfloat k = 0.0, k_ref = 0.0;
    cusfloat w = 0.0;
    for (int i=0; i<ref_data.num_depth; i++)
    {
        for (int j=0; j<ref_data.num_period; j++)
        {
            // Calculate angular frequency to feed 
            // dispersion relation function
            w = 2*PI/ref_data.wave_period[j];

            // Calculate wave number using relation dispersion
            k = w2k(w, ref_data.water_depth[i], ref_data.grav_acc);

            // Check wave number with reference data
            k_ref = ref_data.wave_number[i*ref_data.num_period+j];
            diff = std::abs(k-k_ref);
            if (diff>WAVENUM_TOL)
            {
                std::cerr << "The wave number associated to the wave period: ";
                std::cerr << ref_data.wave_period[j] << " s and depth: " << ref_data.water_depth[i]; 
                std::cerr << " m \nhas been calculated ";
                std::cerr << "out of precision." << std::endl;
                std::cerr << "   -> k_target: " << k_ref << " - ki: " << k;
                std::cerr << " - diff: " << diff << std::endl;
                std::cerr << " - k_target: " << k_ref << " - Diff" << pow2s(w)/ref_data.grav_acc-k_ref*std::tanh(k_ref*ref_data.water_depth[i]) << std::endl;
                std::cerr << " - ki: " << k << " - Diff" << pow2s(w)/ref_data.grav_acc-k*std::tanh(k*ref_data.water_depth[i]) << std::endl;
                pass = false;
            }
        }
    }

}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path_real(argv[1]);

    // Declare logic system variables
    bool pass = false;
    int return_chn = 0;

    // 
    pass = test_dispersion_real(file_path_real);
    // cusfloat T = 10.0;
    // cusfloat h = 100;
    // cusfloat w = 2*PI/T;
    // cusfloat k = w2k(w, h, 9.81);
    // int n = 10;
    // cusfloat* k = new cusfloat[n];

    // w2ki(w, h, 9.81, n, k);
    // for (int i=0; i<n; i++)
    // {
    //     std::cout << "n[" << i << "]: " << k[i] << std::endl;
    // }

    // delete [] k;

    return 0;
}