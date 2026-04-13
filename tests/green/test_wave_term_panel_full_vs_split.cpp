/*
 * @file    test_wave_term_panel_full_vs_split.cpp
 * @brief   Validates the full wave-term integration against split wave+steady assembly.
 *
 * This test integrates the finite-depth Green function over a single quadrilateral panel
 * using Gauss points and compares two equivalent paths:
 *  - Full path: wave_term_integral(..., is_full=true) integrated over the panel.
 *  - Split path: wave_term_integral(..., is_full=false) integrated over the panel plus
 *    the exact steady contribution from Newman's source formulation.
 *
 * It also checks that the steady increment recovered from (full - wave-only) matches the
 * exact steady term assembled with the same regular-frequency image-source convention.
 */

// Include general usage libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/containers.hpp"
#include "../../src/green/integrals_database.hpp"
#include "../../src/green/pulsating_fin_depth.hpp"
#include "../../src/green/source.hpp"
#include "../../src/math/integration.hpp"
#include "../../src/math/math_tools.hpp"


struct PanelInfluenceResult
{
    cuscomplex phi         = cuscomplex(0.0, 0.0);
    cuscomplex velocity[3] = { cuscomplex(0.0, 0.0), cuscomplex(0.0, 0.0), cuscomplex(0.0, 0.0) };
};


void define_submerged_square_panel(
                                    PanelGeom& panel,
                                    cusfloat   scale,
                                    cusfloat   z0
                                )
{
    panel.num_nodes = 4;

    panel.x[0] = -0.5 * scale;
    panel.x[1] = -0.5 * scale;
    panel.x[2] =  0.5 * scale;
    panel.x[3] =  0.5 * scale;

    panel.y[0] = -0.5 * scale;
    panel.y[1] =  0.5 * scale;
    panel.y[2] =  0.5 * scale;
    panel.y[3] = -0.5 * scale;

    panel.z[0] = z0;
    panel.z[1] = z0;
    panel.z[2] = z0;
    panel.z[3] = z0;
}


void initialize_panel(
                        PanelGeom& panel,
                        cusfloat*  cog
                    )
{
    panel.calculate_properties(cog);
    panel.calculate_integration_properties<NUM_GP>();
}


PanelInfluenceResult calculate_exact_steady_regular(
                                                        PanelGeom&      panel,
                                                        PanelGeom&      panel_mirror,
                                                        const cusfloat* field_point,
                                                        cusfloat        water_depth
                                                    )
{
    PanelInfluenceResult result;

    cusfloat vel_0[3] = {0.0, 0.0, 0.0};
    cusfloat vel_1[3] = {0.0, 0.0, 0.0};
    cusfloat vel_2[3] = {0.0, 0.0, 0.0};
    cusfloat vel_3[3] = {0.0, 0.0, 0.0};
    cusfloat vel_4[3] = {0.0, 0.0, 0.0};
    cusfloat vel_5[3] = {0.0, 0.0, 0.0};

    cusfloat pot_0 = 0.0;
    cusfloat pot_1 = 0.0;
    cusfloat pot_2 = 0.0;
    cusfloat pot_3 = 0.0;
    cusfloat pot_4 = 0.0;
    cusfloat pot_5 = 0.0;

    cusfloat shifted_field_point[3] = {0.0, 0.0, 0.0};

    calculate_source_newman(&panel, const_cast<cusfloat*>(field_point), 0, 0, vel_0, pot_0);

    shifted_field_point[0] = field_point[0];
    shifted_field_point[1] = field_point[1];
    shifted_field_point[2] = field_point[2] + 2.0 * water_depth;
    calculate_source_newman(&panel_mirror, shifted_field_point, 0, 0, vel_1, pot_1);

    shifted_field_point[0] = field_point[0];
    shifted_field_point[1] = field_point[1];
    shifted_field_point[2] = field_point[2];
    calculate_source_newman(&panel_mirror, shifted_field_point, 0, 0, vel_2, pot_2);

    shifted_field_point[0] = field_point[0];
    shifted_field_point[1] = field_point[1];
    shifted_field_point[2] = field_point[2] + 2.0 * water_depth;
    calculate_source_newman(&panel, shifted_field_point, 0, 0, vel_3, pot_3);

    shifted_field_point[0] = field_point[0];
    shifted_field_point[1] = field_point[1];
    shifted_field_point[2] = -field_point[2] + 2.0 * water_depth;
    calculate_source_newman(&panel_mirror, shifted_field_point, 0, 0, vel_4, pot_4);

    shifted_field_point[0] = field_point[0];
    shifted_field_point[1] = field_point[1];
    shifted_field_point[2] = field_point[2] + 4.0 * water_depth;
    calculate_source_newman(&panel_mirror, shifted_field_point, 0, 0, vel_5, pot_5);

    result.phi         = cuscomplex((pot_0 + pot_1 + pot_2 + pot_3 + pot_4 + pot_5) / (4.0 * PI), 0.0);

    result.velocity[0] = cuscomplex(-(vel_0[0] + vel_1[0] + vel_2[0] + vel_3[0] + vel_4[0] + vel_5[0]) / (4.0 * PI), 0.0);
    result.velocity[1] = cuscomplex(-(vel_0[1] + vel_1[1] + vel_2[1] + vel_3[1] + vel_4[1] + vel_5[1]) / (4.0 * PI), 0.0);
    result.velocity[2] = cuscomplex(-(vel_0[2] + vel_1[2] + vel_2[2] + vel_3[2] - vel_4[2] + vel_5[2]) / (4.0 * PI), 0.0);

    return result;
}


template<int N>
PanelInfluenceResult integrate_wave_panel(
                                            PanelGeom&               panel,
                                            const cusfloat*          field_point,
                                            cusfloat                 water_depth,
                                            BesselFactoryVecUpTo<N>& bessel_factory,
                                            WaveDispersionFONK&      wave_data,
                                            bool                    is_full
                                        )
{
    cusfloat dX[N];
    cusfloat dY[N];
    cusfloat R[N];
    cusfloat z[N];
    cusfloat zeta[N];

    cuscomplex G[N];
    cuscomplex G_dr[N];
    cuscomplex G_dz[N];
    cuscomplex G_dzeta[N];

    for (int i = 0; i < N; i++)
    {
        dX[i]   = field_point[0] - panel.gauss_points_global_x[i];
        dY[i]   = field_point[1] - panel.gauss_points_global_y[i];
        R[i]    = std::sqrt(pow2s(dX[i]) + pow2s(dY[i]));
        z[i]    = field_point[2];
        zeta[i] = panel.gauss_points_global_z[i];
    }

    wave_term_integral<N, G_ON, DGDR_ON, DGDZ_ON, FSLID_OFF>(
                                                                R,
                                                                z,
                                                                zeta,
                                                                water_depth,
                                                                bessel_factory,
                                                                wave_data,
                                                                G,
                                                                G_dr,
                                                                G_dz,
                                                                G_dzeta,
                                                                is_full
                                                            );

    cuscomplex dG_dx[N];
    cuscomplex dG_dy[N];
    for (int i = 0; i < N; i++)
    {
        cusfloat safe_R = std::max(R[i], cusfloat(1e-12));
        dG_dx[i]        = G_dr[i] * dX[i] / safe_R;
        dG_dy[i]        = G_dr[i] * dY[i] / safe_R;
    }

    PanelInfluenceResult result;
    gauss2d_loop<NUM_GP>(result.phi,         G,     &panel);
    gauss2d_loop<NUM_GP>(result.velocity[0], dG_dx, &panel);
    gauss2d_loop<NUM_GP>(result.velocity[1], dG_dy, &panel);
    gauss2d_loop<NUM_GP>(result.velocity[2], G_dz, &panel);

    result.phi         /= (4.0 * PI);
    result.velocity[0] /= (4.0 * PI);
    result.velocity[1] /= (4.0 * PI);
    result.velocity[2] /= (4.0 * PI);

    return result;
}


cusfloat complex_abs(const cuscomplex& value)
{
    return std::abs(value);
}


cusfloat complex_vector_norm(const cuscomplex* values)
{
    cusfloat norm_sq = 0.0;
    for (int i = 0; i < 3; i++)
    {
        norm_sq += pow2s(complex_abs(values[i]));
    }
    return std::sqrt(norm_sq);
}


bool compare_case(
                    const char*                 label,
                    const PanelInfluenceResult& calculated,
                    const PanelInfluenceResult& reference,
                    cusfloat                    max_phi_rel_error,
                    cusfloat                    max_vel_rel_error,
                    bool                        verbose = false
                )
{
    cusfloat phi_rel_error   = complex_abs(calculated.phi - reference.phi) /
                               std::max(complex_abs(reference.phi), cusfloat(1e-12));
    cuscomplex vel_diff[3]   = {
                                  calculated.velocity[0] - reference.velocity[0],
                                  calculated.velocity[1] - reference.velocity[1],
                                  calculated.velocity[2] - reference.velocity[2]
                              };
    cusfloat vel_rel_error   = complex_vector_norm(vel_diff) /
                               std::max(complex_vector_norm(reference.velocity), cusfloat(1e-12));

    if ( phi_rel_error > max_phi_rel_error || vel_rel_error > max_vel_rel_error || verbose )
    {
        std::cout << label << " comparison:" << std::endl;
        std::cout << "    phi calc = " << calculated.phi << " ref = " << reference.phi << std::endl;
        std::cout << "    vel calc = (" << calculated.velocity[0] << ", " << calculated.velocity[1] << ", " << calculated.velocity[2] << ")" << std::endl;
        std::cout << "    vel ref  = (" << reference.velocity[0] << ", " << reference.velocity[1] << ", " << reference.velocity[2] << ")" << std::endl;
        std::cout << "    phi rel err = " << phi_rel_error << " max = " << max_phi_rel_error << std::endl;
        std::cout << "    vel rel err = " << vel_rel_error << " max = " << max_vel_rel_error << std::endl;
    }

    if (phi_rel_error > max_phi_rel_error || vel_rel_error > max_vel_rel_error)
    {
        std::cerr << label << " failed" << std::endl;
        return false;
    }

    return true;
}


bool run_field_point_case(
                            PanelGeom&                    panel,
                            PanelGeom&                    panel_mirror,
                            const cusfloat*               field_point,
                            cusfloat                      water_depth,
                            BesselFactoryVecUpTo<NUM_GP2>& bessel_factory,
                            WaveDispersionFONK&           wave_data,
                            cusfloat                      max_total_phi_rel_error,
                            cusfloat                      max_total_vel_rel_error,
                            cusfloat                      max_steady_phi_rel_error,
                            cusfloat                      max_steady_vel_rel_error
                        )
{
    PanelInfluenceResult steady_exact = calculate_exact_steady_regular(panel, panel_mirror, field_point, water_depth);
    PanelInfluenceResult wave_only    = integrate_wave_panel<NUM_GP2>(panel, field_point, water_depth, bessel_factory, wave_data, false);
    PanelInfluenceResult full         = integrate_wave_panel<NUM_GP2>(panel, field_point, water_depth, bessel_factory, wave_data, true);

    PanelInfluenceResult total_split;
    total_split.phi        = steady_exact.phi + wave_only.phi;
    for (int i = 0; i < 3; i++)
    {
        total_split.velocity[i] = steady_exact.velocity[i] + wave_only.velocity[i];
    }

    PanelInfluenceResult steady_from_full;
    steady_from_full.phi = full.phi - wave_only.phi;
    for (int i = 0; i < 3; i++)
    {
        steady_from_full.velocity[i] = full.velocity[i] - wave_only.velocity[i];
    }

    if (!compare_case("Full vs split total", full, total_split, max_total_phi_rel_error, max_total_vel_rel_error, true))
    {
        return false;
    }

    if (!compare_case("Recovered steady vs Newman steady", steady_from_full, steady_exact, max_steady_phi_rel_error, max_steady_vel_rel_error, true))
    {
        return false;
    }

    return true;
}


int main(void)
{
    cusfloat cog[3] = {0.0, 0.0, 0.0};
    PanelGeom panel;
    PanelGeom panel_mirror;

    define_submerged_square_panel(panel, 2.0, -1.0);
    define_submerged_square_panel(panel_mirror, 2.0, 1.0);
    initialize_panel(panel, cog);
    initialize_panel(panel_mirror, cog);

    const cusfloat water_depth = 50.0;
    const cusfloat grav_acc    = 9.81;
    const cusfloat period      = 16.129032258064516;
    const cusfloat ang_freq    = 2.0 * PI / period;
    const cusfloat H           = pow2s(ang_freq) / grav_acc * water_depth;

    fold_database(H);

    WaveDispersionFONK wave_data(ang_freq, water_depth, grav_acc);
    wave_data.update_full(ang_freq, water_depth, grav_acc);
    BesselFactoryVecUpTo<NUM_GP2> bessel_factory;

    struct FieldPointCase
    {
        cusfloat point[3];
        cusfloat total_phi_tol;
        cusfloat total_vel_tol;
        cusfloat steady_phi_tol;
        cusfloat steady_vel_tol;
    };

    const FieldPointCase cases[] = {
        {{ 4.0,   0.0,  -1.0 }, 0.08, 0.10, 0.08, 0.10},
        {{ 0.0,  -4.0,  -1.0 }, 0.08, 0.10, 0.08, 0.10},
        {{ 1.5,   0.5,  -0.25}, 0.12, 0.15, 0.12, 0.15},
        {{-1.5,  -0.75, -2.5 }, 0.08, 0.10, 0.08, 0.10}
    };

    for (const auto& item : cases)
    {
        if (!run_field_point_case(
                                    panel,
                                    panel_mirror,
                                    item.point,
                                    water_depth,
                                    bessel_factory,
                                    wave_data,
                                    item.total_phi_tol,
                                    item.total_vel_tol,
                                    item.steady_phi_tol,
                                    item.steady_vel_tol
                                ))
        {
            return 1;
        }
    }

    return 0;
}