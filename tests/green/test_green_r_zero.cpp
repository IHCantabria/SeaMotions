
// Include local modules
#include "../../src/config.hpp"
#include "../../src/green/integrals_db.hpp"
#include "../../src/green/pulsating_fin_depth.hpp"
#include "../../src/waves.hpp"


void    analytical_derivative_steady(
                                        cusfloat*       field_point,
                                        cusfloat        x,
                                        cusfloat        y,
                                        cusfloat        z,
                                        cusfloat        water_depth,
                                        cuscomplex&     dG_dX,
                                        cuscomplex&     dG_dY
                                    )
{
    cusfloat    R       =   std::sqrt(
                                        pow2s( field_point[0] - x )
                                        +
                                        pow2s( field_point[1] - y )
                                    );
    
    // Calculate derivative with respect to R
    cuscomplex  dG_dR   = G_integral_steady_dr(
                                                        R,
                                                        field_point[2],
                                                        z,
                                                        water_depth
                                                    );
    
    // Calculate X and Y cartesian coordinates derivatives
    cusfloat    dX      =   field_point[0] - x;
                dG_dX   =   dG_dR * dX / R;
    cusfloat    dY      =   field_point[1] - y;
                dG_dY   =   dG_dR * dY / R;
}


void    analytical_derivative_wave(
                                        cusfloat*           field_point,
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        cusfloat            water_depth,
                                        WaveDispersionFO* wave_data,
                                        IntegralsDb*        integrals_db,
                                        cuscomplex&         dG_dX,
                                        cuscomplex&         dG_dY
                                    )
{
    cusfloat    R       =   std::sqrt(
                                        pow2s( field_point[0] - x )
                                        +
                                        pow2s( field_point[1] - y )
                                    );
    
    // Calculate derivative with respect to R
    cuscomplex  dG_dR   = G_integral_wave_dr(
                                                        R,
                                                        field_point[2],
                                                        z,
                                                        water_depth,
                                                        *wave_data,
                                                        *integrals_db
                                                    );
    
    // Calculate X and Y cartesian coordinates derivatives
    cusfloat    dX      =   field_point[0] - x;
                dG_dX   =   dG_dR * dX / R;
    cusfloat    dY      =   field_point[1] - y;
                dG_dY   =   dG_dR * dY / R;
}


void    numerical_derivative_steady(
                                        cusfloat*       field_point,
                                        cusfloat        x,
                                        cusfloat        y,
                                        cusfloat        z,
                                        cusfloat        water_depth,
                                        cuscomplex&     dG_dX,
                                        cuscomplex&     dG_dY
                                    )
{
    // Define derivation interval
    cusfloat    eps     = 1e-6;

    // Calculate dG_dX
    cusfloat    x0      = x - eps / 2.0;
    cusfloat    x1      = x + eps / 2.0;

    cuscomplex  G0x     = G_integral_steady(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x0,
                                                y,
                                                z,
                                                water_depth
                                            );

    cuscomplex  G1x     = G_integral_steady(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x1,
                                                y,
                                                z,
                                                water_depth
                                            );

                dG_dX   = ( G1x - G0x ) / eps;

    // Calculate dG_dY
    cusfloat    y0      = y - eps / 2.0;
    cusfloat    y1      = y + eps / 2.0;

    cuscomplex  G0y     = G_integral_steady(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x,
                                                y0,
                                                z,
                                                water_depth
                                            );

    cuscomplex  G1y     = G_integral_steady(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x,
                                                y1,
                                                z,
                                                water_depth
                                            );

                dG_dY   = ( G1y - G0y ) / eps;
}


void    numerical_derivative_wave(
                                        cusfloat*           field_point,
                                        cusfloat            x,
                                        cusfloat            y,
                                        cusfloat            z,
                                        cusfloat            water_depth,
                                        WaveDispersionFO* wave_data,
                                        IntegralsDb*        integrals_db,
                                        cuscomplex&         dG_dX,
                                        cuscomplex&         dG_dY
                                    )
{
    // Define derivation interval
    cusfloat    eps     = 1e-6;

    // Calculate dG_dX
    cusfloat    x0      = x - eps / 2.0;
    cusfloat    x1      = x + eps / 2.0;

    cuscomplex  G0x     = G_integral_wave(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x0,
                                                y,
                                                z,
                                                water_depth,
                                                *wave_data,
                                                *integrals_db
                                            );

    cuscomplex  G1x     = G_integral_wave(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x1,
                                                y,
                                                z,
                                                water_depth,
                                                *wave_data,
                                                *integrals_db
                                            );

                dG_dX   = ( G1x - G0x ) / eps;

    // Calculate dG_dY
    cusfloat    y0      = y - eps / 2.0;
    cusfloat    y1      = y + eps / 2.0;

    cuscomplex  G0y     =   G_integral_wave(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x,
                                                y0,
                                                z,
                                                water_depth,
                                                *wave_data,
                                                *integrals_db
                                            );

    cuscomplex  G1y     =   G_integral_wave(
                                                field_point[0],
                                                field_point[1],
                                                field_point[2],
                                                x,
                                                y1,
                                                z,
                                                water_depth,
                                                *wave_data,
                                                *integrals_db
                                            );

                dG_dY   = ( G1y - G0y ) / eps;
}


void    test_steady( 
                                        cusfloat        water_depth
                    )
{
    /***********************************************/
    /********* Perform tests for R=0 ***************/
    /***********************************************/

    // Define field and source point positions
    cusfloat    field_point[3]  =   { 1e-3, 0.0, -9.0 };
    cusfloat    source_point[3] =   { 0.0, 0.0, -10.0 };

    // Calculate analytical derivatives
    cuscomplex  dG_dX_ana( 0.0, 0.0 ) ;
    cuscomplex  dG_dY_ana( 0.0, 0.0 ) ;

    analytical_derivative_steady(
                                    field_point,
                                    source_point[0],
                                    source_point[1],
                                    source_point[2],
                                    water_depth,
                                    dG_dX_ana,
                                    dG_dY_ana
                                );

    // Calculate analytical derivatives
    cuscomplex  dG_dX_num( 0.0, 0.0 ) ;
    cuscomplex  dG_dY_num( 0.0, 0.0 ) ;

    numerical_derivative_steady(
                                    field_point,
                                    source_point[0],
                                    source_point[1],
                                    source_point[2],
                                    water_depth,
                                    dG_dX_num,
                                    dG_dY_num
                                );

    // Show comparative
    std::cout << "GREEN FUNCTION STEADY TERM: " << std::endl;
    std::cout << "--> dG_dX_ana: " << dG_dX_ana << " - dG_dX_num: " << dG_dX_num << std::endl;
    std::cout << "--> dG_dY_ana: " << dG_dY_ana << " - dG_dY_num: " << dG_dY_num << std::endl;

}


void    test_wave( 
                                        cusfloat        ang_freq,
                                        cusfloat        water_depth,
                                        cusfloat        grav_acc
                )
{
    /***********************************************/
    /*************** Load Wave Data ****************/
    /***********************************************/
    // Calculate wave numbers
    WaveDispersionFO* wave_data   = new WaveDispersionFO( 
                                                                ang_freq,
                                                                30,
                                                                water_depth,
                                                                grav_acc
                                                            );
    wave_data->calculate_john_terms( );

    // Load integrals database
    IntegralsDb*    integrals_db    = new IntegralsDb( );
    
    // Fold for current frequency and water depth
    cusfloat H = pow2s( ang_freq ) * water_depth / grav_acc;
    integrals_db->fold_h( H );

    /***********************************************/
    /********* Perform tests for R=0 ***************/
    /***********************************************/

    // Define field and source point positions
    cusfloat    field_point[3]  =   { 1e-3, 0.0, -9.0 };
    cusfloat    source_point[3] =   { 0.0, 0.0, -10.0 };

    // Calculate analytical derivatives
    cuscomplex  dG_dX_ana( 0.0, 0.0 ) ;
    cuscomplex  dG_dY_ana( 0.0, 0.0 ) ;

    analytical_derivative_wave(
                                    field_point,
                                    source_point[0],
                                    source_point[1],
                                    source_point[2],
                                    water_depth,
                                    wave_data,
                                    integrals_db,
                                    dG_dX_ana,
                                    dG_dY_ana
                                );

    // Calculate analytical derivatives
    cuscomplex  dG_dX_num( 0.0, 0.0 ) ;
    cuscomplex  dG_dY_num( 0.0, 0.0 ) ;

    numerical_derivative_wave(
                                    field_point,
                                    source_point[0],
                                    source_point[1],
                                    source_point[2],
                                    water_depth,
                                    wave_data,
                                    integrals_db,
                                    dG_dX_num,
                                    dG_dY_num
                            );

    // Show comparative
    std::cout << "GREEN FUNCTION WAVE TERM: " << std::endl;
    std::cout << "dG_dX_ana: " << dG_dX_ana << " - dG_dX_num: " << dG_dX_num << std::endl;
    std::cout << "dG_dY_ana: " << dG_dY_ana << " - dG_dY_num: " << dG_dY_num << std::endl;

}


int main( void )
{
    // Define problem properties
    cusfloat    ang_freq    = 0.1;
    cusfloat    grav_acc    = 9.81;
    cusfloat    water_depth = 1000.0;

    // Test derivatives over Steady term
    test_steady(
                    grav_acc
                );

    // Test derivatives over Wave term
    test_wave(
                    ang_freq,
                    water_depth,
                    grav_acc
            );

    return 0;
}