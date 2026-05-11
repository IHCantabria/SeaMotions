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
 * Test: test_generalized_alpha_nl
 * --------------------------------
 * Validates the nonlinear Generalized-Alpha solver (GeneralizedAlphaNL) through
 * two test cases.
 *
 * -------------------------------------------------------------------------
 * Case 1 — Linear SDOF through the nonlinear interface
 * -------------------------------------------------------------------------
 * System:
 *   m * u'' + c * u' + k * u = 0      (free damped oscillation)
 *
 * Parameters:
 *   m = 1.0 kg,  omega_n = 2*pi rad/s,  zeta = 0.10
 *   k = m*omega_n^2,  c = 2*m*omega_n*zeta
 *   u(0) = 1.0 m,  v(0) = 0 m/s
 *
 * The internal force residual is expressed as:
 *   R_fcn = (1-alpha_f)*(k*u_{n+1} + c*v_{n+1})
 *          + alpha_f   *(k*u_n     + c*v_n    )
 *          - F_ext     (= 0, free vibration)
 *
 * The mass inertia contribution is added by the solver itself.
 *
 * Tangent functors return the constant linear operators:
 *   K_T = k,   C_T = c
 *
 * Verification:
 *   - Numerical solution matches the analytical response within ABS_TOL = 1e-3.
 *   - Newton-Raphson converges in exactly 1 iteration at every step
 *     (quadratic convergence collapses to one step for a linear problem).
 *
 * Analytical solution (free damped oscillator, u(0)=1, v(0)=0):
 *   omega_d = omega_n * sqrt(1 - zeta^2)
 *   u(t)    = exp(-zeta*omega_n*t) * [cos(omega_d*t) + (zeta/sqrt(1-zeta^2))*sin(omega_d*t)]
 *
 * -------------------------------------------------------------------------
 * Case 2 — Duffing oscillator (nonlinear stiffness)
 * -------------------------------------------------------------------------
 * System:
 *   m * u'' + c * u' + k * u + k_nl * u^3 = 0
 *
 * Parameters:
 *   m = 1.0 kg,  k = (2*pi)^2 N/m,  c = 2*m*omega_n*0.05 (5% damping)
 *   k_nl = 10 * k  (strong nonlinearity)
 *   u(0) = 1.0 m,  v(0) = 0 m/s
 *
 * The residual at the alpha-f level is:
 *   R_fcn = (1-alpha_f)*[k*u_{n+1} + k_nl*u_{n+1}^3 + c*v_{n+1}]
 *          + alpha_f   *[k*u_n     + k_nl*u_n^3     + c*v_n    ]
 *
 * Tangent operators:
 *   K_T = k + 3*k_nl*u^2   (secant of the cubic spring)
 *   C_T = c
 *
 * Verification:
 *   - NR converges at every step (nr_converged = true throughout).
 *   - Total mechanical energy E = 0.5*m*v^2 + 0.5*k*u^2 + 0.25*k_nl*u^4
 *     is monotonically decreasing (damping is positive).
 *
 * Integration settings (both cases):
 *   dt = 0.001 s,  t_end = 3.0 s,  rho_inf = 0.5
 */

// Include general usage libraries
#include <iostream>
#include <cmath>
#include <string>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_tools.hpp"
#include "../../src/math/generalized_alpha_nl.hpp"

//--------------------------------------------------------------------
//-- Test tolerances
//--------------------------------------------------------------------
constexpr cusfloat ABS_TOL = 1e-3;     // Position error tolerance

//--------------------------------------------------------------------
//-- Helper: build a 1x1 CSRMatrix from a scalar value
//--------------------------------------------------------------------
CSRMatrix* build_csr_1x1( cusfloat value )
{
    cusfloat dense[1] = { value };
    return new CSRMatrix( 1, dense );
}

//--------------------------------------------------------------------
//-- Analytical solution (free damped SDOF, u0=1, v0=0)
//--------------------------------------------------------------------
cusfloat sdof_analytical_pos( cusfloat t, cusfloat omega_n, cusfloat zeta )
{
    cusfloat omega_d = omega_n * std::sqrt( 1.0 - zeta*zeta );
    cusfloat decay   = std::exp( -zeta * omega_n * t );
    return decay * ( std::cos( omega_d * t )
                   + ( zeta / std::sqrt( 1.0 - zeta*zeta ) ) * std::sin( omega_d * t ) );
}

//--------------------------------------------------------------------
//-- Dummy linear force functor (unused by GeneralizedAlphaNL but
//-- required as TLinear template parameter)
//--------------------------------------------------------------------
struct DummyLinear
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
//-- Case 1: Linear SDOF residual functor
//--------------------------------------------------------------------
/**
 * @brief Residual contribution (internal forces minus external forces),
 *        excluding the inertia term that GeneralizedAlphaNL adds itself.
 *
 * R_fcn = (1-alpha_f)*(k*u_{n+1} + c*v_{n+1})
 *        + alpha_f   *(k*u_n     + c*v_n    )
 */
struct LinearSDOFResidual
{
    cusfloat k;
    cusfloat c;
    cusfloat alpha_f;

    void operator()(
                        cusfloat    /* t        */,
                        cusfloat    /* dt       */,
                        cusfloat*   u_pos,
                        cusfloat*   u_vel,
                        cusfloat*   /* u_acc   */,
                        cusfloat*   u_pos_old,
                        cusfloat*   u_vel_old,
                        cusfloat*   /* u_acc_old */,
                        cusfloat*   residual
                    )
    {
        residual[0] = ( 1.0 - alpha_f ) * ( k * u_pos[0] + c * u_vel[0] )
                    +          alpha_f   * ( k * u_pos_old[0] + c * u_vel_old[0] );
    }
};

//--------------------------------------------------------------------
//-- Case 1: Linear SDOF tangent functor
//--------------------------------------------------------------------
struct LinearSDOFTangent
{
    cusfloat k;
    cusfloat c;

    void operator()(
                        cusfloat    /* t    */,
                        cusfloat    /* dt   */,
                        cusfloat*   /* pos  */,
                        cusfloat*   /* vel  */,
                        cusfloat*   /* acc  */,
                        CSRMatrix*& K_T,
                        CSRMatrix*& C_T
                    )
    {
        cusfloat k_dense[1] = { k };
        cusfloat c_dense[1] = { c };
        K_T = new CSRMatrix( 1, k_dense );
        C_T = new CSRMatrix( 1, c_dense );
    }
};

//--------------------------------------------------------------------
//-- Case 1: Test
//--------------------------------------------------------------------
bool test_linear_sdof_through_nl_interface( cusfloat rho_inf )
{
    bool pass = true;

    //----------------------------------------------------------------
    // System parameters
    //----------------------------------------------------------------
    constexpr cusfloat m       = 1.0;
    constexpr cusfloat omega_n = 2.0 * 3.141592653589793;   // 2*pi rad/s
    constexpr cusfloat zeta    = 0.10;
    const     cusfloat k       = m * omega_n * omega_n;
    const     cusfloat c       = 2.0 * m * omega_n * zeta;

    //----------------------------------------------------------------
    // Build 1×1 mass matrix
    //----------------------------------------------------------------
    CSRMatrix* M_csr = build_csr_1x1( m );

    //----------------------------------------------------------------
    // Integration settings
    //----------------------------------------------------------------
    constexpr cusfloat dt    = 0.001;
    constexpr cusfloat t0    = 0.0;
    constexpr cusfloat t_end = 3.0;

    //----------------------------------------------------------------
    // Initial conditions
    //----------------------------------------------------------------
    cusfloat u0 = 1.0;
    cusfloat v0 = 0.0;
    cusfloat a0 = 0.0;

    //----------------------------------------------------------------
    // Compute alpha_f for this rho_inf so that the residual functor
    // can blend the forces at the correct interpolation level
    //----------------------------------------------------------------
    cusfloat alpha_f_val = rho_inf / ( rho_inf + 1.0 );

    //----------------------------------------------------------------
    // Build functors
    //----------------------------------------------------------------
    LinearSDOFResidual  res_fcn = { k, c, alpha_f_val };
    LinearSDOFTangent   tan_fcn = { k, c };
    DummyLinear         lin_fcn;

    //----------------------------------------------------------------
    // Build the solver (no DOF restrictions)
    //----------------------------------------------------------------
    GeneralizedAlphaNL<LinearSDOFResidual, LinearSDOFTangent, DummyLinear> solver(
                                                                                        &res_fcn,
                                                                                        &tan_fcn,
                                                                                        &lin_fcn,
                                                                                        M_csr,
                                                                                        dt,
                                                                                        t0,
                                                                                        rho_inf,
                                                                                        50,
                                                                                        1e-10,
                                                                                        1e-8,
                                                                                        &u0,
                                                                                        &v0,
                                                                                        &a0
                                                                                    );

    std::cout << "  rho_inf = " << rho_inf
              << " | alpha_m = " << solver.alpha_m
              << " | alpha_f = " << solver.alpha_f
              << " | gamma = "   << solver.gamma
              << " | beta = "    << solver.beta << "\n";

    //----------------------------------------------------------------
    // Define check times
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
    bool nr_always_1iter = true;

    for ( int step=0; step<total_steps; step++ )
    {
        solver.step();

        // For a linear problem, NR should converge in exactly 1 iteration
        if ( solver.nr_iter > 1 )
        {
            nr_always_1iter = false;
        }
        if ( !solver.nr_converged )
        {
            std::cout << "  [WARN] NR did not converge at step " << step << "\n";
            pass = false;
        }

        // Check solution at requested times
        if ( next_check < n_checks && (step+1) == check_steps[next_check] )
        {
            cusfloat t_check = check_times[next_check];
            cusfloat u_num   = solver.y_pos[0];
            cusfloat u_ref   = sdof_analytical_pos( t_check, omega_n, zeta );
            cusfloat abs_err = std::abs( u_num - u_ref );

            std::cout << "  t = " << t_check << " s"
                      << " | u_num = " << u_num
                      << " | u_ref = " << u_ref
                      << " | |error| = " << abs_err
                      << " | NR iter = " << solver.nr_iter;

            if ( abs_err > ABS_TOL )
            {
                std::cout << " -> FAIL (tol = " << ABS_TOL << ")\n";
                pass = false;
            }
            else
            {
                std::cout << " -> PASS\n";
            }

            next_check++;
        }
    }

    if ( nr_always_1iter )
    {
        std::cout << "  NR iter = 1 at every step (linear problem) -> PASS\n";
    }
    else
    {
        std::cout << "  NR iter > 1 detected (unexpected for linear problem) -> FAIL\n";
        pass = false;
    }

    delete M_csr;

    return pass;
}


//--------------------------------------------------------------------
//-- Case 2: Duffing oscillator residual functor
//--------------------------------------------------------------------
/**
 * @brief Residual contribution for the Duffing oscillator:
 *
 *   m * u'' + c * u' + k * u + k_nl * u^3 = 0
 *
 * R_fcn = (1-alpha_f)*[k*u_{n+1} + k_nl*u_{n+1}^3 + c*v_{n+1}]
 *        + alpha_f   *[k*u_n     + k_nl*u_n^3     + c*v_n    ]
 */
struct DuffingResidual
{
    cusfloat k;
    cusfloat k_nl;
    cusfloat c;
    cusfloat alpha_f;

    void operator()(
                        cusfloat    /* t        */,
                        cusfloat    /* dt       */,
                        cusfloat*   u_pos,
                        cusfloat*   u_vel,
                        cusfloat*   /* u_acc   */,
                        cusfloat*   u_pos_old,
                        cusfloat*   u_vel_old,
                        cusfloat*   /* u_acc_old */,
                        cusfloat*   residual
                    )
    {
        cusfloat f_int_new = k * u_pos[0]     + k_nl * u_pos[0]*u_pos[0]*u_pos[0]     + c * u_vel[0];
        cusfloat f_int_old = k * u_pos_old[0] + k_nl * u_pos_old[0]*u_pos_old[0]*u_pos_old[0] + c * u_vel_old[0];
        residual[0] = ( 1.0 - alpha_f ) * f_int_new + alpha_f * f_int_old;
    }
};

//--------------------------------------------------------------------
//-- Case 2: Duffing oscillator tangent functor
//--------------------------------------------------------------------
struct DuffingTangent
{
    cusfloat k;
    cusfloat k_nl;
    cusfloat c;

    void operator()(
                        cusfloat    /* t    */,
                        cusfloat    /* dt   */,
                        cusfloat*   u_pos,
                        cusfloat*   /* vel  */,
                        cusfloat*   /* acc  */,
                        CSRMatrix*& K_T,
                        CSRMatrix*& C_T
                    )
    {
        cusfloat k_t_val    = k + 3.0 * k_nl * u_pos[0]*u_pos[0];
        cusfloat k_t[1]     = { k_t_val };
        cusfloat c_t[1]     = { c };
        K_T = new CSRMatrix( 1, k_t );
        C_T = new CSRMatrix( 1, c_t );
    }
};

//--------------------------------------------------------------------
//-- Case 2: Test
//--------------------------------------------------------------------
bool test_duffing_oscillator( cusfloat rho_inf )
{
    bool pass = true;

    //----------------------------------------------------------------
    // System parameters
    //----------------------------------------------------------------
    constexpr cusfloat m       = 1.0;
    constexpr cusfloat omega_n = 2.0 * 3.141592653589793;   // 2*pi rad/s
    constexpr cusfloat zeta    = 0.05;
    const     cusfloat k       = m * omega_n * omega_n;
    const     cusfloat k_nl    = 10.0 * k;   // strong cubic nonlinearity
    const     cusfloat c       = 2.0 * m * omega_n * zeta;

    //----------------------------------------------------------------
    // Build 1×1 mass matrix
    //----------------------------------------------------------------
    CSRMatrix* M_csr = build_csr_1x1( m );

    //----------------------------------------------------------------
    // Integration settings
    //----------------------------------------------------------------
    constexpr cusfloat dt    = 0.001;
    constexpr cusfloat t0    = 0.0;
    constexpr cusfloat t_end = 3.0;

    //----------------------------------------------------------------
    // Initial conditions
    //----------------------------------------------------------------
    cusfloat u0 = 1.0;
    cusfloat v0 = 0.0;
    cusfloat a0 = 0.0;

    cusfloat alpha_f_val = rho_inf / ( rho_inf + 1.0 );

    //----------------------------------------------------------------
    // Build functors
    //----------------------------------------------------------------
    DuffingResidual res_fcn = { k, k_nl, c, alpha_f_val };
    DuffingTangent  tan_fcn = { k, k_nl, c };
    DummyLinear     lin_fcn;

    //----------------------------------------------------------------
    // Build the solver
    //----------------------------------------------------------------
    GeneralizedAlphaNL<DuffingResidual, DuffingTangent, DummyLinear> solver(
                                                                                    &res_fcn,
                                                                                    &tan_fcn,
                                                                                    &lin_fcn,
                                                                                    M_csr,
                                                                                    dt,
                                                                                    t0,
                                                                                    rho_inf,
                                                                                    50,
                                                                                    1e-10,
                                                                                    1e-8,
                                                                                    &u0,
                                                                                    &v0,
                                                                                    &a0
                                                                                );

    std::cout << "  rho_inf = " << rho_inf
              << " | alpha_m = " << solver.alpha_m
              << " | alpha_f = " << solver.alpha_f << "\n";

    //----------------------------------------------------------------
    // Total mechanical energy: E = 0.5*m*v^2 + 0.5*k*u^2 + 0.25*k_nl*u^4
    //----------------------------------------------------------------
    auto energy = [&]( cusfloat u, cusfloat v ) -> cusfloat {
        return 0.5*m*v*v + 0.5*k*u*u + 0.25*k_nl*u*u*u*u;
    };

    cusfloat E_prev = energy( solver.y_pos[0], solver.y_vel[0] );
    cusfloat E_max  = E_prev;

    //----------------------------------------------------------------
    // Time integration loop
    //----------------------------------------------------------------
    int total_steps = static_cast<int>( std::round( t_end / dt ) );

    // Report energy and NR status at coarse intervals
    constexpr int n_reports = 6;
    const cusfloat report_times[n_reports] = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 };
    int report_steps[n_reports];
    for ( int i=0; i<n_reports; i++ )
    {
        report_steps[i] = static_cast<int>( std::round( report_times[i] / dt ) );
    }
    int next_report = 0;

    for ( int step=0; step<total_steps; step++ )
    {
        solver.step();

        if ( !solver.nr_converged )
        {
            std::cout << "  [FAIL] NR did not converge at step " << step
                      << " (t = " << solver.time << " s)"
                      << " | ||R|| = " << solver.nr_res_abs << "\n";
            pass = false;
        }

        // Report at coarse output times
        if ( next_report < n_reports && (step+1) == report_steps[next_report] )
        {
            cusfloat E_cur = energy( solver.y_pos[0], solver.y_vel[0] );

            std::cout << "  t = " << report_times[next_report] << " s"
                      << " | u = "    << solver.y_pos[0]
                      << " | E = "    << E_cur
                      << " | NR iter = " << solver.nr_iter
                      << " | ||R|| = " << solver.nr_res_abs;

            // Energy must not increase (positive damping)
            if ( E_cur > E_max * ( 1.0 + 1e-6 ) )   // small tolerance for rounding
            {
                std::cout << " -> FAIL (energy increased)\n";
                pass = false;
            }
            else
            {
                std::cout << " -> PASS\n";
            }

            E_max       = std::max( E_max, E_cur );
            next_report++;
        }
    }

    // Final check: solution has decayed significantly from initial energy
    cusfloat E_final  = energy( solver.y_pos[0], solver.y_vel[0] );
    cusfloat E_init   = energy( u0, v0 );
    cusfloat E_ratio  = E_final / E_init;

    std::cout << "  Energy ratio E_final/E_init = " << E_ratio;
    if ( E_ratio < 0.99 )
    {
        std::cout << " (energy decayed as expected) -> PASS\n";
    }
    else
    {
        std::cout << " (energy did not decay as expected) -> FAIL\n";
        pass = false;
    }

    delete M_csr;

    return pass;
}


//--------------------------------------------------------------------
//-- main
//--------------------------------------------------------------------
int main( void )
{
    std::cout << "\n==========================================================\n";
    std::cout << "  Test: GeneralizedAlphaNL - 1-DOF nonlinear solver\n";
    std::cout << "==========================================================\n\n";

    bool all_pass = true;

    //----------------------------------------------------------------
    // Case 1: Linear SDOF expressed through the nonlinear interface.
    // NR should collapse to 1 iteration at every step, and the solution
    // must match the analytical response within ABS_TOL.
    //----------------------------------------------------------------
    std::cout << "--- Case 1: Linear SDOF through NL interface (rho_inf = 0.5) ---\n";
    bool pass1 = test_linear_sdof_through_nl_interface( 0.5 );
    all_pass = all_pass && pass1;
    std::cout << "  Case 1 result: " << (pass1 ? "PASS" : "FAIL") << "\n\n";

    //----------------------------------------------------------------
    // Case 2: Duffing oscillator with strong cubic nonlinearity.
    // NR must converge at every step and energy must monotonically decay.
    //----------------------------------------------------------------
    std::cout << "--- Case 2: Duffing oscillator (rho_inf = 0.5) ---\n";
    bool pass2 = test_duffing_oscillator( 0.5 );
    all_pass = all_pass && pass2;
    std::cout << "  Case 2 result: " << (pass2 ? "PASS" : "FAIL") << "\n\n";

    std::cout << "----------------------------------------------------------\n";
    std::cout << "  Overall result: " << (all_pass ? "PASS" : "FAIL") << "\n";
    std::cout << "==========================================================\n\n";

    return all_pass ? 0 : 1;
}
