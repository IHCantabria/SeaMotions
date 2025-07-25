
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
#include "../../src/green/pulsating_fin_depth_v2.hpp"
#include "../../src/green/integrals_database.hpp"
#include "../../src/math/bessel_factory.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/tools.hpp"
#include "../../src/waves/wave_dispersion_fo.hpp"

// Load namespace for inline complex literals
using namespace std::literals::complex_literals;


// Define test tolerances
cusfloat JOHN_TOL = 1e-3;


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
            this->H[i] += 1e-9;
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


template<int fnc_type, int mode_f, int mode_dfdr, int mode_dfdz>
bool launch_integral(
                        std::string file_path, 
                        std::string function_type,
                        int         function_pos
                    )
{
    // Define test result flag
    bool pass = true;

    // Read refernce data
    RefData ref_data;
    ref_data.load_data(file_path);
    cusfloat g = ref_data.grav_acc;
    cusfloat h = ref_data.water_depth;

    // Create an instance of Bessel Factory class
    constexpr std::size_t N = 1;
    BesselFactoryVecUpTo<N> bessel_factory;

    // Loop over refence data to check over all parameter
    // space
    cusfloat            Avec[N];
    cuscomplex          Gp;
    cuscomplex          G[N];
    cuscomplex          G_dr[N];
    cuscomplex          G_dz[N];
    cusfloat            zvec[N];
    cusfloat            zetavec[N];
    cuscomplex          jc, jr;
    cusfloat            nu = 0.0;
    constexpr int       num_kn = 30;
    int                 rd_index = 0;
    cusfloat            w = 0.1;
    WaveDispersionFONK  wave_data(w, h, g);

    for (int i=0; i<ref_data.num_H; i++)
    {
        // Calculate dependent variables on H parameter
        nu = ref_data.H[i]/h;
        w = std::sqrt(nu*g);
        wave_data.update_full( w, h, g );

        // Fold integrals coefficients
        if constexpr( fnc_type == 1 )
        {
            fold_database( ref_data.H[i] );
            std::cout << "Database folded!" << std::endl;
        }

        std::cout << "I: " << i << std::endl;
        for (int j=0; j<ref_data.num_A; j++)
        {
            for (int k=0; k<ref_data.num_z; k++)
            {
                // Storage scalar values into vector ones
                Avec[0]     = ref_data.A[j] * ref_data.water_depth;
                zvec[0]     = ref_data.z[k];
                zetavec[0]  = ref_data.zeta[k];

                // Calculate index data
                rd_index = (
                                i*(ref_data.num_H*ref_data.num_A)
                                + j*ref_data.num_A
                                + k
                            );

                if ( rd_index == 300 )
                {
                    double aa = 0.0;
                }
                std::cout << "I: " << i << " - J: " << j << " - K: " << k << std::endl;
                // Calculate Green function integral value
                // wave_term_fin_depth_integral<N>( N, Avec, zvec, zetavec, h, bessel_factory, wave_data, G, G_dr, G_dz );
                if constexpr( fnc_type == 0)
                {
                    Gp = john_series(Avec[0], zvec[0], zetavec[0], h , wave_data);
                    john_series<N, mode_f, mode_dfdr, mode_dfdz>( Avec, zvec, zetavec, h, bessel_factory, wave_data, G, G_dr, G_dz );

                    // std::cout << "G: " << G[0] << " - Gp: " << Gp << std::endl;
                }
                else
                {
                    wave_term_integral<N, mode_f, mode_dfdr, mode_dfdz>( Avec, zvec, zetavec, h, bessel_factory, wave_data, G, G_dr, G_dz );
                }
                // custom_template<N>( Avec );

                if ( function_pos == 0 )
                {
                    jc = G[0];
                }
                else if ( function_pos == 1 )
                {
                    jc = G_dr[0];
                }
                else
                {
                    jc = G_dz[0];
                }

                // We compare against the conjugate
                // vector as the reference values are ouputed with time sign exp(iwt) while 
                // Seamotions works with exp(-iwt) (waves goes from -X to +X)
                jr = std::conj( ref_data.g_series[rd_index] );

                // Compare values and check with tolerance
                if (!assert_complex_equality(jc, jr, JOHN_TOL))
                {
                    std::cerr << "Integral method " << function_type;
                    std::cerr <<" does not converge to the expected value " << std::endl;
                    std::cerr << "for the following input parameters: " << std::endl;
                    std::cerr << "  - R/h: " << ref_data.A[j] << std::endl;
                    std::cerr << "  - R: " << ref_data.A[j]*h << std::endl;
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
                    std::cerr << " - Index: " << rd_index << std::endl;
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
    if (!check_num_cmd_args(argc, 6))
    {
        return 1;
    }

    // std::string file_path_john(argv[1]);
    // std::string file_path_john_dr(argv[2]);
    // std::string file_path_john_dz(argv[3]);
    // std::string file_path_Gint(argv[4]);
    // std::string file_path_Gint_dr(argv[5]);
    // std::string file_path_Gint_dz(argv[6]);

    std::string file_path_john( "D:/sergio/developments/SeaMotions/tests/tests_data/green_fin_depth_tables/G_Agt1.dat" );
    std::string file_path_john_dr( "D:/sergio/developments/SeaMotions/tests/tests_data/green_fin_depth_tables/G_Agt1_dR.dat" );
    std::string file_path_john_dz( "D:/sergio/developments/SeaMotions/tests/tests_data/green_fin_depth_tables/G_Agt1_dz.dat" );
    std::string file_path_Gint( "D:/sergio/developments/SeaMotions/tests/tests_data/green_fin_depth_tables/G_Alt1.dat" );
    std::string file_path_Gint_dr( "D:/sergio/developments/SeaMotions/tests/tests_data/green_fin_depth_tables/G_Alt1_dR.dat" );
    std::string file_path_Gint_dz( "D:/sergio/developments/SeaMotions/tests/tests_data/green_fin_depth_tables/G_Alt1_dz.dat" );
    

    // Declare variable for the logic system
    bool pass = false;

    // Launch test to check the John eigenfunction
    // expansion
    std::cout << "here" << std::endl;
    pass = launch_integral<0, G_ON, DGDR_OFF, DGDZ_OFF>(
                                                            file_path_john, 
                                                            "G",
                                                            0
                                                        );
    if (!pass)
    {
        return 1;
    }

    pass = launch_integral<0, G_OFF, DGDR_ON, DGDZ_OFF>(
                                                            file_path_john_dr, 
                                                            "dG_dr",
                                                            1
                                                        );
    if (!pass)
    {
        return 1;
    }

    pass = launch_integral<0, G_OFF, DGDR_OFF, DGDZ_ON>(
                                                            file_path_john_dz,
                                                            "dG_dz",
                                                            2
                                                        );
    if (!pass)
    {
        return 1;
    }

    // Launch test to check the Green function
    // integral approximation method
    pass = launch_integral<1, G_ON, DGDR_OFF, DGDZ_OFF>(
                                                            file_path_Gint, 
                                                            "G_integral",
                                                            0
                                                        );
    if (!pass)
    {
        return 1;
    }

    pass = launch_integral<1, G_OFF, DGDR_ON, DGDZ_OFF>(
                                                            file_path_Gint_dr, 
                                                            "G_integral_dr",
                                                            1
                                                        );
    if (!pass)
    {
        return 1;
    }

    pass = launch_integral<1, G_OFF, DGDR_OFF, DGDZ_ON>(
                                                            file_path_Gint_dz, 
                                                            "G_integral_dz",
                                                            2
                                                        );
    if (!pass)
    {
        return 1;
    }

    return 0;
}