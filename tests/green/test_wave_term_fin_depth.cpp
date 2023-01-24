
// Include general usage libraries
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

// Include general usage scientific libraries
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/pulsating_fin_depth.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/tools.hpp"
#include "../../src/waves.hpp"

// Load namespace for inline complex literals
using namespace std::literals::complex_literals;


// Define test tolerances
cusfloat JOHN_TOL = 1e-7;


struct RefData
{
    cusfloat* A = nullptr;
    cuscomplex* g_series = nullptr;
    cusfloat grav_acc = 0.0;
    cusfloat* H = nullptr;
    int num_A = 0;
    int num_H = 0;
    int num_z = 0;
    cusfloat water_depth = 0;
    cusfloat* z = nullptr;
    cusfloat* zeta = nullptr;

    ~RefData()
    {
        delete [] this->A;
        delete [] this->H;
    }

    void load_data(std::string file_path)
    {
        // Declare auxiliar string to read the file line by line
        std::string str_aux;

        // Open file unit
        std::ifstream infile(file_path);

        // Read Acceleration of the Gravity
        infile >> str_aux >> this->grav_acc;

        // Read water depth
        infile >> str_aux >> this->water_depth;

        // Read A parameter space
        infile >> str_aux >> this->num_A;
        this->A = new cusfloat [this->num_A];
        for (int i=0; i<this->num_A; i++)
        {
            infile >> this->A[i];
        }

        // Read H parameter space
        infile >> str_aux >> this->num_H;
        this->H = new cusfloat [this->num_H];
        for (int i=0; i<this->num_H; i++)
        {
            infile >> this->H[i];
        }

        // Read Z parameter space
        infile >> str_aux >> this->num_z;
        this->z = new cusfloat [this->num_z];
        this->zeta = new cusfloat [this->num_z];
        for (int i=0; i<this->num_z; i++)
        {
            infile >> this->z[i];
        }
        for (int i=0; i<this->num_z; i++)
        {
            infile >> this->zeta[i];
        }

        // Read Green series representation values
        infile >> str_aux;
        int num_gs = this->num_A*this->num_H*this->num_z;
        this->g_series = new cuscomplex [num_gs];
        cusfloat a=1.0, b=0.0;
        for (int i=0; i<num_gs; i++)
        {
            // Read real (a) and imaginary (b) parts
            infile >> a >> b;

            // Storage as a complex number
            this->g_series[i] = a + b*1i;
        }

        // Close file unit
        infile.close();
    }

    void print(void)
    {
        // Print-out file content
        std::cout << "Num.A: " << this->num_A << std::endl;
        for (int i=0; i<this->num_A; i++)
        {
            std::cout << this->A[i] << std::endl;
        }

        std::cout << "Num.H: " << this->num_H << std::endl;
        for (int i=0; i<this->num_H; i++)
        {
            std::cout << this->H[i] << std::endl;
        }

        std::cout << "Num.Z: " << this->num_z << std::endl;
        for (int i=0; i<this->num_z; i++)
        {
            std::cout << this->z[i] << std::endl;
        }
        std::cout << "----------" << std::endl;
        for (int i=0; i<this->num_z; i++)
        {
            std::cout << this->zeta[i] << std::endl;
        }

        int num_gs = this->num_A*this->num_H*this->num_z;
        for (int i=0; i<num_gs; i++)
        {
            std::cout << (i+10) << " " << this->g_series[i] << std::endl;
        }
    }
};


bool launch_john(
                std::string file_path, 
                std::function <cuscomplex(
                                        cusfloat,
                                        cusfloat,
                                        cusfloat,
                                        cusfloat,
                                        WaveDispersionData&
                                        )> f_def,
                std::string function_type
                )
{
    // Define test result flag
    bool pass = true;

    // Read refernce data
    RefData ref_data;
    ref_data.load_data(file_path);
    cusfloat g = ref_data.grav_acc;
    cusfloat h = ref_data.water_depth;

    // Loop over refence data to check over all parameter
    // space
    cuscomplex jc, jr;
    cusfloat k0 = 0.0;
    cusfloat nu = 0.0;
    const int num_kn = 30;
    int rd_index = 0;
    cusfloat w = 0.0;
    for (int i=0; i<ref_data.num_A; i++)
    {
        for (int j=0; j<ref_data.num_H; j++)
        {
            // Calculate dependent variables on H parameter
            nu = ref_data.H[j]/h;
            w = std::sqrt(nu*g);
            WaveDispersionData wave_data(w, num_kn, h, g);
            wave_data.calculate_john_terms();

            for (int k=0; k<ref_data.num_z; k++)
            {
                // Calculate John series value
                jc = f_def(
                    ref_data.A[i]*ref_data.water_depth,
                    ref_data.z[k],
                    ref_data.zeta[k],
                    ref_data.water_depth,
                    wave_data
                    );

                // Get John series reference data
                rd_index = (
                            i*(ref_data.num_H*ref_data.num_z)
                            + j*ref_data.num_z
                            + k
                            );
                jr = ref_data.g_series[rd_index];

                // Compare values and check with tolerance
                if (!assert_complex_equality(jc, jr, JOHN_TOL))
                {
                    std::cerr << "John series " << function_type;
                    std::cerr <<" does not converge to the expected value " << std::endl;
                    std::cerr << "for the following input parameters: " << std::endl;
                    std::cerr << "  - R/h: " << ref_data.A[i] << std::endl;
                    std::cerr << "  - R: " << ref_data.A[i]*h << std::endl;
                    std::cerr << "  - z: " << ref_data.z[k] << std::endl;
                    std::cerr << "  - zeta: " << ref_data.zeta[k] << std::endl;
                    std::cerr << "  - h: " << h << std::endl;
                    std::cerr << "  - nu: " << nu << std::endl;
                    std::cerr << "  - k0: " << k0 << std::endl;
                    std::cerr << "  - num_kn: " << num_kn << std::endl;
                    std::cerr << "Numerical error description: " << std::endl;
                    std::cerr << " - Expected value: " << jr << std::endl;
                    std::cerr << " - Calculated value: " << jc << std::endl;
                    std::cerr << " - Difference value: " << jc-jr << std::endl;
                    pass = false;
                    goto exit;
                }
            }
        }
    }

    exit:
        return pass;
}


int main(int argc, char* argv[])
{
    // Read command line arguments
    if (!check_num_cmd_args(argc, 1))
    {
        return 1;
    }

    std::string file_path_john(argv[1]);

    // Declare variable for the logic system
    bool pass = false;
    int return_chn = 0;

    // Launch test to check the John eigenfunction
    // expansion
    pass = launch_john(file_path_john, john_series, "G");
    if (!pass)
    {
        return_chn = 1;
    }

    return return_chn;
}