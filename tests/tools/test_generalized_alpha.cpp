/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * Test: test_generalized_alpha
 * ----------------------------
 * Validates the Generalized-Alpha time solver against the analytical solution
 * of a 1-DOF damped harmonic oscillator:
 *
 *   m * u'' + c * u' + k * u = 0
 *
 * System parameters:
 *   m       = 1.0 kg
 *   omega_n = 2*pi rad/s  (natural frequency, T_n = 1 s)
 *   zeta    = 0.10        (10% critical damping)
 *   k       = m * omega_n^2 = 4*pi^2 N/m
 *   c       = 2 * m * omega_n * zeta = 0.4*pi N·s/m
 *
 * Initial conditions:
 *   u(0) = 1.0 m,  u'(0) = 0.0 m/s
 *
 * Analytical solution (free damped vibration):
 *   omega_d = omega_n * sqrt(1 - zeta^2)
 *   u(t)    = exp(-zeta * omega_n * t) * [cos(omega_d*t) + (zeta/sqrt(1-zeta^2))*sin(omega_d*t)]
 *
 * Integration settings:
 *   dt      = 0.001 s
 *   t_end   = 3.0 s  (3 natural periods)
 *   rho_inf tested: 1.0 (no dissipation, equals avg-acceleration) and 0.5 (moderate dissipation)
 *
 * Pass criterion: |u_numerical - u_analytical| < ABS_TOL at check times.
 *
 * Note: with rho_inf = 1 the Generalized-Alpha scheme reduces exactly to the
 * constant-average-acceleration Newmark method, so results must match. With
 * rho_inf < 1 high-frequency modes are damped and amplitude errors increase
 * slightly, so ABS_TOL is relaxed to 5e-3 for those cases.
 */

// Include general usage libraries
#include <iostream>
#include <cmath>
#include <string>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/generalized_alpha.hpp"

//--------------------------------------------------------------------
//-- Test tolerances
//--------------------------------------------------------------------
// rho_inf = 1.0 → identical to Newmark avg-acceleration: tight tolerance
constexpr cusfloat ABS_TOL_NO_DISS = 1e-3;
// rho_inf = 0.5 → moderate numerical dissipation: slightly relaxed tolerance
constexpr cusfloat ABS_TOL_DISS    = 5e-3;

//--------------------------------------------------------------------
//-- External force functor (zero force - free vibration)
//--------------------------------------------------------------------
struct ZeroFext
{
    void operator()(
                        cusfloat    /* t   */,
                        cusfloat    /* dt  */,
                        cusfloat*   /* pos */,
                        cusfloat*   /* vel */,
                        cusfloat*   /* acc */,
                        cusfloat*   fext
                    )
    {
        fext[0] = 0.0;
    }
};

//--------------------------------------------------------------------
//-- Analytical solution
//--------------------------------------------------------------------
/**
 * @brief Analytical position for a free damped 1-DOF oscillator.
 *
 * u(t) = exp(-zeta*omega_n*t) * [cos(omega_d*t) + (zeta/sqrt(1-zeta^2))*sin(omega_d*t)]
 * with initial conditions u(0)=1, v(0)=0.
 *
 * @param t       Evaluation time [s]
 * @param omega_n Natural frequency [rad/s]
 * @param zeta    Damping ratio [-]
 */
cusfloat sdof_analytical_pos( cusfloat t, cusfloat omega_n, cusfloat zeta )
{
    cusfloat omega_d = omega_n * std::sqrt( 1.0 - zeta*zeta );
    cusfloat decay   = std::exp( -zeta * omega_n * t );
    return decay * ( std::cos( omega_d * t ) + ( zeta / std::sqrt( 1.0 - zeta*zeta ) ) * std::sin( omega_d * t ) );
}

//--------------------------------------------------------------------
//-- Helper: build a 1x1 CSRMatrix from a scalar value
//--------------------------------------------------------------------
CSRMatrix* build_csr_1x1( cusfloat value )
{
    cusfloat dense[1] = { value };
    return new CSRMatrix( 1, dense );
}

//--------------------------------------------------------------------
//-- Core test function for a given rho_inf value
//--------------------------------------------------------------------
/**
 * @brief Run the Generalized-Alpha solver for a given rho_inf and check
 *        the solution against the analytical SDOF response.
 *
 * @param rho_inf   Spectral radius at infinite frequency ∈ [0, 1]
 * @param abs_tol   Absolute error tolerance for the position checks
 */
bool test_generalized_alpha_sdof( cusfloat rho_inf, cusfloat abs_tol )
{
    bool pass = true;

    //----------------------------------------------------------------
    // System parameters
    //----------------------------------------------------------------
    constexpr cusfloat m       = 1.0;
    constexpr cusfloat omega_n = 2.0 * 3.141592653589793;   // 2*pi rad/s -> T_n = 1 s
    constexpr cusfloat zeta    = 0.10;
    const     cusfloat k       = m * omega_n * omega_n;
    const     cusfloat c       = 2.0 * m * omega_n * zeta;

    //----------------------------------------------------------------
    // Build 1×1 sparse system matrices
    //----------------------------------------------------------------
    CSRMatrix* M_csr = build_csr_1x1( m );
    CSRMatrix* K_csr = build_csr_1x1( k );
    CSRMatrix* C_csr = build_csr_1x1( c );

    //----------------------------------------------------------------
    // Integration settings
    //----------------------------------------------------------------
    constexpr cusfloat dt    = 0.001;
    constexpr cusfloat t0    = 0.0;
    constexpr cusfloat t_end = 3.0;

    //----------------------------------------------------------------
    // Initial conditions (acceleration is corrected by the solver)
    //----------------------------------------------------------------
    cusfloat u0 = 1.0;
    cusfloat v0 = 0.0;
    cusfloat a0 = 0.0;

    //----------------------------------------------------------------
    // Restrictions (none imposed: all DOFs are free)
    //----------------------------------------------------------------
    int restrictions[1] = { 0 };

    //----------------------------------------------------------------
    // External force functor
    //----------------------------------------------------------------
    ZeroFext fcn;

    //----------------------------------------------------------------
    // Build the Generalized-Alpha solver
    // The solver corrects a0 internally via M*a0 = F - K*u0 - C*v0
    //----------------------------------------------------------------
    GeneralizedAlpha<ZeroFext*> solver(
                                            &fcn,
                                            M_csr,
                                            K_csr,
                                            C_csr,
                                            dt,
                                            t0,
                                            rho_inf,
                                            &u0,
                                            &v0,
                                            &a0,
                                            restrictions
                                        );

    //----------------------------------------------------------------
    // Report computed algorithm parameters
    //----------------------------------------------------------------
    std::cout << "  rho_inf = " << rho_inf
              << " | alpha_m = " << solver.alpha_m
              << " | alpha_f = " << solver.alpha_f
              << " | gamma = "   << solver.gamma
              << " | beta = "    << solver.beta << "\n";

    //----------------------------------------------------------------
    // Define check times and their step indices
    //----------------------------------------------------------------
    constexpr int n_checks = 6;
    const cusfloat check_times[n_checks] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };

    int check_steps[n_checks];
    for ( int i=0; i<n_checks; i++ )
    {
        check_steps[i] = static_cast<int>( std::round( check_times[i] / dt ) );
    }

    //----------------------------------------------------------------
    // Time integration loop
    //----------------------------------------------------------------
    int total_steps = static_cast<int>( std::round( t_end / dt ) );
    int next_check  = 0;

    for ( int step=0; step<total_steps; step++ )
    {
        solver.step();

        // Check at requested times
        if ( next_check < n_checks && (step+1) == check_steps[next_check] )
        {
            cusfloat t_check  = check_times[next_check];
            cusfloat u_num    = solver.y_pos[0];
            cusfloat u_ref    = sdof_analytical_pos( t_check, omega_n, zeta );
            cusfloat abs_err  = std::abs( u_num - u_ref );

            std::cout << "  t = " << t_check << " s"
                      << " | u_num = " << u_num
                      << " | u_ref = " << u_ref
                      << " | |error| = " << abs_err;

            if ( abs_err > abs_tol )
            {
                std::cout << " -> FAIL (tol = " << abs_tol << ")\n";
                pass = false;
            }
            else
            {
                std::cout << " -> PASS\n";
            }

            next_check++;
        }
    }

    //----------------------------------------------------------------
    // Cleanup
    //----------------------------------------------------------------
    delete M_csr;
    delete K_csr;
    delete C_csr;

    return pass;
}


int main( void )
{
    std::cout << "\n====================================================\n";
    std::cout << "  Test: GeneralizedAlpha - 1-DOF damped oscillator\n";
    std::cout << "====================================================\n\n";

    bool all_pass = true;

    //----------------------------------------------------------------
    // Case 1: rho_inf = 1.0
    // No numerical dissipation. Reduces exactly to the Newmark
    // constant-average-acceleration method (beta=1/4, gamma=1/2).
    //----------------------------------------------------------------
    std::cout << "--- Case 1: rho_inf = 1.0 (no numerical dissipation) ---\n";
    bool pass1 = test_generalized_alpha_sdof( 1.0, ABS_TOL_NO_DISS );
    all_pass = all_pass && pass1;
    std::cout << "  Case 1 result: " << (pass1 ? "PASS" : "FAIL") << "\n\n";

    //----------------------------------------------------------------
    // Case 2: rho_inf = 0.5
    // Moderate high-frequency dissipation. Errors are slightly larger
    // due to algorithmic damping of the response frequency itself.
    //----------------------------------------------------------------
    std::cout << "--- Case 2: rho_inf = 0.5 (moderate dissipation) ---\n";
    bool pass2 = test_generalized_alpha_sdof( 0.5, ABS_TOL_DISS );
    all_pass = all_pass && pass2;
    std::cout << "  Case 2 result: " << (pass2 ? "PASS" : "FAIL") << "\n\n";

    //----------------------------------------------------------------
    // Case 3: rho_inf = 0.0
    // Maximum high-frequency dissipation (maximally dissipative scheme).
    // Errors on this well-resolved system (omega_n*dt << 1) are still small.
    //----------------------------------------------------------------
    std::cout << "--- Case 3: rho_inf = 0.0 (maximum dissipation) ---\n";
    bool pass3 = test_generalized_alpha_sdof( 0.0, ABS_TOL_DISS );
    all_pass = all_pass && pass3;
    std::cout << "  Case 3 result: " << (pass3 ? "PASS" : "FAIL") << "\n\n";

    std::cout << "----------------------------------------------------\n";
    std::cout << "  Overall result: " << (all_pass ? "PASS" : "FAIL") << "\n";
    std::cout << "----------------------------------------------------\n\n";

    return all_pass ? 0 : 1;
}
