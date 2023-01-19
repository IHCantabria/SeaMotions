
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
    int num_kn = 0;
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

    void load_data_imag(std::string file_path)
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
        
        // Read number of imaginary solutions per combination
        // of water depth and wave period
        infile >> str_aux >> this->num_kn;

        // Read Wave numbers
        infile >> str_aux;
        this->wave_number = new cusfloat[this->num_period*this->num_depth*this->num_kn];
        int index = 0;
        for (int i=0; i<this->num_depth; i++)
        {
            for (int j=0; j<this->num_period; j++)
            {
                for (int k=0; k<this->num_kn; k++)
                {
                    index = (
                        i*(this->num_period*this->num_kn)
                        +j*this->num_kn
                        +k
                        );
                    infile >> this->wave_number[index];
                }
            }
        }

        // Close file unit
        infile.close();
    }

    void load_data_real(std::string file_path)
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


bool test_dispersion_imag(std::string file_path)
{
    // Declare test variable
    bool pass = true;

    // Load reference data
    RefData ref_data;
    ref_data.load_data_imag(file_path);

    // Allocate heap space for the imaginary wave numbers
    // solution
    cusfloat* kn = new cusfloat [ref_data.num_kn];

    // Loop over reference data to check
    // the calculated results are in agreement
    cusfloat diff = 0.0;
    cusfloat g = 9.81;
    int index_ref = 0;
    cusfloat k_ref = 0.0;
    cusfloat w = 0.0;
    for (int i=0; i<ref_data.num_depth; i++)
    {
        for (int j=0; j<ref_data.num_period; j++)
        {
            // Calculate imaginary wave numbers
            w = 2*PI/ref_data.wave_period[j];
            w2ki(w, ref_data.water_depth[i], g, ref_data.num_kn, kn);

            // Loop over results to compare with the reference ones
            for (int k=0; k<ref_data.num_kn; k++)
            {
                index_ref = (
                    i*(ref_data.num_period*ref_data.num_kn)
                    + j*ref_data.num_kn
                    + k
                    );
                k_ref = ref_data.wave_number[index_ref];
                diff = std::abs(kn[k]-k_ref);
                if (diff > WAVENUM_TOL)
                {
                    std::cerr << "The imaginary wave number associated to the wave period: ";
                    std::cerr << ref_data.wave_period[j] << " s, depth: " << ref_data.water_depth[i]; 
                    std::cerr << " m and branch: " << k <<" \nhas been calculated ";
                    std::cerr << "out of precision." << std::endl;
                    std::cerr << "   -> k_target: " << k_ref << " - ki: " << kn[k];
                    std::cerr << " - diff: " << diff << std::endl;
                    std::cerr << " - k_target: " << k_ref << " - Diff" << pow2s(w)/ref_data.grav_acc+k_ref*std::tan(k_ref*ref_data.water_depth[i]) << std::endl;
                    std::cerr << " - ki: " << kn[k] << " - Diff" << pow2s(w)/ref_data.grav_acc+kn[k]*std::tan(kn[k]*ref_data.water_depth[i]) << std::endl;
                    pass = false;
                    goto exit;
                }
            }
        }
    }

    exit:
        // Delete heap memory
        delete [] kn;

        return pass;
}


bool test_dispersion_real(std::string file_path)
{
    // Declare test variable
    bool pass = true;
    
    // Read reference data
    RefData ref_data;
    ref_data.load_data_real(file_path);

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
                std::cerr << "The real wave number associated to the wave period: ";
                std::cerr << ref_data.wave_period[j] << " s and depth: " << ref_data.water_depth[i]; 
                std::cerr << " m \nhas been calculated ";
                std::cerr << "out of precision." << std::endl;
                std::cerr << "   -> k_target: " << k_ref << " - ki: " << k;
                std::cerr << " - diff: " << diff << std::endl;
                std::cerr << " - k_target: " << k_ref << " - Diff" << pow2s(w)/ref_data.grav_acc-k_ref*std::tanh(k_ref*ref_data.water_depth[i]) << std::endl;
                std::cerr << " - ki: " << k << " - Diff" << pow2s(w)/ref_data.grav_acc-k*std::tanh(k*ref_data.water_depth[i]) << std::endl;
                pass = false;
                goto exit;
            }
        }
    }

    exit:
        return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 2))
    {
        return 1;
    }

    std::string file_path_real(argv[1]);
    std::string file_path_imag(argv[2]);

    // Declare logic system variables
    bool pass = false;
    int return_chn = 0;

    // Launch checking the real roots finding
    pass = test_dispersion_real(file_path_real);
    if (!pass)
    {
        return_chn = 1;
    }

    // Launch checking the imaginary roots finding
    pass = test_dispersion_imag(file_path_imag);
    if (!pass)
    {
        return_chn = 1;
    }

    return return_chn;
}