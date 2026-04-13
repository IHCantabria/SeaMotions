// Include general usage libraries
#include <algorithm>
#include <cmath>
#include <iostream>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/containers/containers.hpp"
#include "../../src/green/source.hpp"
#include "../../src/math/gauss_t.hpp"
#include "../../src/math/math_tools.hpp"


struct GaussSourceResult
{
    cusfloat phi = 0.0;
    cusfloat velocity[3] = {0.0, 0.0, 0.0};
};

void define_square_panel(PanelGeom &panel, cusfloat scale)
{
    // Define number of nodes for the panel
    panel.num_nodes = 4;

    // Define X position
    panel.x[0] = -0.5*scale;
    panel.x[1] = -0.5*scale;
    panel.x[2] = 0.5*scale;
    panel.x[3] = 0.5*scale;

    // Define Y position
    panel.y[0] = -0.5*scale;
    panel.y[1] = 0.5*scale;
    panel.y[2] = 0.5*scale;
    panel.y[3] = -0.5*scale;

    // Define Z position
    panel.z[0] = 0.0;
    panel.z[1] = 0.0;
    panel.z[2] = 0.0;
    panel.z[3] = 0.0;
}


template<int NGP>
GaussSourceResult integrate_source_gauss(const PanelGeom& panel, const cusfloat* field_point)
{
    GaussSourceResult result;

    for (int i = 0; i < NGP * NGP; i++)
    {
        cusfloat r_vec[3];
        r_vec[0] = field_point[0] - panel.gauss_points_global_x[i];
        r_vec[1] = field_point[1] - panel.gauss_points_global_y[i];
        r_vec[2] = field_point[2] - panel.gauss_points_global_z[i];

        cusfloat r = std::sqrt(pow2s(r_vec[0]) + pow2s(r_vec[1]) + pow2s(r_vec[2]));

        cusfloat weight = GaussPointsT<2, NGP>::weights_x[i] *
                          GaussPointsT<2, NGP>::weights_y[i] *
                          panel.jac_det_gauss_points[i];

        result.phi += weight / r;

        cusfloat inv_r3 = 1.0 / (r * r * r);
        result.velocity[0] += weight * r_vec[0] * inv_r3;
        result.velocity[1] += weight * r_vec[1] * inv_r3;
        result.velocity[2] += weight * r_vec[2] * inv_r3;
    }

    return result;
}


cusfloat vector_norm(const cusfloat* v)
{
    return std::sqrt(pow2s(v[0]) + pow2s(v[1]) + pow2s(v[2]));
}


bool run_distance_case(
                        PanelGeom& panel,
                        const cusfloat* normal,
                        cusfloat ratio,
                        cusfloat max_rel_phi,
                        cusfloat max_rel_vel
                    )
{
    cusfloat field_point[3];
    for (int i = 0; i < 3; i++)
    {
        field_point[i] = panel.center[i] + normal[i] * ratio * panel.length;
    }

    cusfloat phi_exact = 0.0;
    cusfloat vel_exact[3] = {0.0, 0.0, 0.0};
    calculate_source_newman(
                                &panel,
                                field_point,
                                0,
                                0,
                                vel_exact,
                                phi_exact
                            );

    auto gauss = integrate_source_gauss<NUM_GP>(panel, field_point);

    cusfloat phi_abs_err = std::abs(phi_exact - gauss.phi);
    cusfloat phi_ref = std::max(std::abs(phi_exact), cusfloat(1e-12));
    cusfloat phi_rel_err = phi_abs_err / phi_ref;

    cusfloat vel_abs_err_vec[3] = {
                                    std::abs(vel_exact[0] - gauss.velocity[0]),
                                    std::abs(vel_exact[1] - gauss.velocity[1]),
                                    std::abs(vel_exact[2] - gauss.velocity[2])
                                };
    cusfloat vel_abs_err = vector_norm(vel_abs_err_vec);
    cusfloat vel_ref = std::max(vector_norm(vel_exact), cusfloat(1e-12));
    cusfloat vel_rel_err = vel_abs_err / vel_ref;

    if (phi_rel_err > max_rel_phi || vel_rel_err > max_rel_vel)
    {
        std::cerr << "Gauss integration failed for r/length = " << ratio << std::endl;
        std::cerr << "    phi_exact = " << phi_exact << ", phi_gauss = " << gauss.phi << std::endl;
        std::cerr << "    rel_err_phi = " << phi_rel_err << ", max = " << max_rel_phi << std::endl;
        std::cerr << "    rel_err_vel = " << vel_rel_err << ", max = " << max_rel_vel << std::endl;
        return false;
    }

    return true;
}


bool run_point_case(
                        PanelGeom& panel,
                        const cusfloat* field_point,
                        const char* label,
                        cusfloat max_rel_phi,
                        cusfloat max_rel_vel
                   )
{
    cusfloat phi_exact = 0.0;
    cusfloat vel_exact[3] = {0.0, 0.0, 0.0};
    calculate_source_newman(
                                &panel,
                                const_cast<cusfloat*>(field_point),
                                0,
                                0,
                                vel_exact,
                                phi_exact
                            );

    auto gauss = integrate_source_gauss<NUM_GP>(panel, field_point);

    cusfloat phi_abs_err = std::abs(phi_exact - gauss.phi);
    cusfloat phi_ref = std::max(std::abs(phi_exact), cusfloat(1e-12));
    cusfloat phi_rel_err = phi_abs_err / phi_ref;

    cusfloat vel_abs_err_vec[3] = {
                                    std::abs(vel_exact[0] - gauss.velocity[0]),
                                    std::abs(vel_exact[1] - gauss.velocity[1]),
                                    std::abs(vel_exact[2] - gauss.velocity[2])
                                };
    cusfloat vel_abs_err = vector_norm(vel_abs_err_vec);
    cusfloat vel_ref = std::max(vector_norm(vel_exact), cusfloat(1e-12));
    cusfloat vel_rel_err = vel_abs_err / vel_ref;

    if (phi_rel_err > max_rel_phi || vel_rel_err > max_rel_vel)
    {
        std::cerr << "Gauss integration failed for " << label << std::endl;
        std::cerr << "    point coords = (" << field_point[0] << ", " << field_point[1] << ", " << field_point[2] << ")" << std::endl;
        std::cerr << "    phi_exact = " << phi_exact << ", phi_gauss = " << gauss.phi << std::endl;
        std::cerr << "    vel_exact = (" << vel_exact[0] << ", " << vel_exact[1] << ", " << vel_exact[2] << ")" << std::endl;
        std::cerr << "    vel_gauss = (" << gauss.velocity[0] << ", " << gauss.velocity[1] << ", " << gauss.velocity[2] << ")" << std::endl;
        std::cerr << "    rel_err_phi = " << phi_rel_err << ", max = " << max_rel_phi << std::endl;
        std::cerr << "    rel_err_vel = " << vel_rel_err << ", max = " << max_rel_vel << std::endl;
        return false;
    }

    return true;
}


int main(void)
{
    PanelGeom panel;
    define_square_panel(panel, 2.0);

    cusfloat cog[3] = {0.0, 0.0, 0.0};
    panel.calculate_properties(cog);
    panel.calculate_integration_properties<2>();

    cusfloat normal[3] = { panel.normal_vec[0], panel.normal_vec[1], panel.normal_vec[2] };

    struct DistanceCase
    {
        cusfloat ratio;
        cusfloat max_rel_phi;
        cusfloat max_rel_vel;
    };

    const DistanceCase cases[] = {
        {0.5, 0.03, 0.06},
        {1.00, 0.03, 0.06},
        {2.00, 0.015, 0.03},
        {4.00, 0.008, 0.02}
    };

    for (const auto& item : cases)
    {
        if (!run_distance_case(panel, normal, item.ratio, item.max_rel_phi, item.max_rel_vel))
        {
            return 1;
        }
    }

    const cusfloat in_plane_far_1[3] = { 1.6, 0.0, 0.0 };
    const cusfloat in_plane_far_2[3] = { 0.0, -1.6, 0.0 };
    const cusfloat out_plane_1[3] = { 1.6, 0.4, 0.5 };
    const cusfloat out_plane_2[3] = { -1.6, -0.4, -0.5 };

    if (!run_point_case(panel, in_plane_far_1, "In-plane point (x=1.6, y=0, z=0)", 0.12, 0.25))
    {
        return 1;
    }
    if (!run_point_case(panel, in_plane_far_2, "In-plane point (x=0, y=-1.6, z=0)", 0.12, 0.25))
    {
        return 1;
    }
    if (!run_point_case(panel, out_plane_1, "Out-of-plane point (0.6, 0.4, 0.5)", 0.08, 0.18))
    {
        return 1;
    }
    if (!run_point_case(panel, out_plane_2, "Out-of-plane point (-0.6, -0.4, -0.5)", 0.08, 0.18))
    {
        return 1;
    }

    return 0;
}
