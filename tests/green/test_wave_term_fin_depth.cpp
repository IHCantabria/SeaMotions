
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
#include "../../src/green/pulsating_fin_depth_cheby.hpp"
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
                    std::cerr << "  - k0: " << wave_data.k0 << std::endl;
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
    if (!check_num_cmd_args(argc, 0))
    {
        return 1;
    }

    // std::string file_path_john(argv[1]);
    // std::string file_path_john_dr(argv[2]);
    // std::string file_path_john_dz(argv[3]);

    // Declare variable for the logic system
    bool pass = false;

    // Launch test to check the John eigenfunction
    // expansion
    // pass = launch_john(file_path_john, john_series, "G");
    // if (!pass)
    // {
    //     return 1;
    // }

    // pass = launch_john(file_path_john_dr, john_series_dr, "dG_dr");
    // if (!pass)
    // {
    //     return 1;
    // }

    // pass = launch_john(file_path_john_dz, john_series_dz, "dG_dz");
    // if (!pass)
    // {
    //     return 1;
    // }
    std::cout << "here" << std::endl;
    // L1 l1 = L1();
    P3 l1 = P3();
    set_data_l1(l1);
    std::cout << "here" << std::endl;
    
    cusfloat A_max = 1.0;
    cusfloat A_min = 0.0;
    cusfloat B_max = 1.0;
    cusfloat B_min = 0.0;
    cusfloat H_max = 0.0;
    cusfloat H_min = -16;
    int num_a = 50;
    int num_b = 50;
    int num_h = 50;

    cusfloat da = (A_max-A_min)/(num_a-1);
    cusfloat db = (B_max-B_min)/(num_b-1);
    cusfloat dh = (H_max-H_min)/(num_h-1);
    cusfloat ai = 0.0;
    cusfloat bi = 0.0;
    cusfloat hi = 0.0;
    cusfloat I = 0.0;
    double elapased_time_mean = 0.0;
    int N = 100;
    for (int n=0; n<N; n++)
    {
        double t0 = get_cpu_time();
        for (int i=0; i<num_a; i++)
        {
            hi = std::pow(10.0, (H_min + i*dh));
            // std::cout << "folding: " << hi << std::endl;
            l1.fold_h(hi);
            // std::cout << " -> Done" << std::endl;

            // std::cout << "i: " << i << " - hi:" << hi << std::endl;
            for (int j=0; j<num_a; j++)
            {
                ai = A_min + j*da;
                // std::cout << "j: " << j << " - ai:" << ai << std::endl;
                for (int k=0; k<num_b; k++)
                {
                    bi = B_min + k*db;
                    // std::cout << "i: " << i << " - j: " << j << " - k: " << k;
                    // std::cout << " - ai:" << ai << " - bi:" << bi << " - hi:" << hi << std::endl;
                    I = l1.get_value_ab(ai, bi);
                    // std::cout << "-> Done: " << I << std::endl;
                }
            }

        }
        double t1 = get_cpu_time();
        elapased_time_mean += (t1-t0);
        // std::cout << "Elapsed Time [s]: " << t1-t0 << std::endl;
    }

    std::cout << "Elapsed Time Mean [s]: " << elapased_time_mean/N << std::endl;
    // l1.fold_h(0.184207);
    // std::cout << "l1: " << l1.get_value_ab(1.0, 1.0) << std::endl;
    // std::cout << "l1: " << l1.get_value_abh(0.1, 0.1, 20.0) << std::endl;
    

    return 0;
}