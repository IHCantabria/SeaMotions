// Include general usage libraries
#include <cmath>
#include <iostream>
#include <functional>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/math_constants.hpp"
#include "../../src/mesh/annulus_nodes.hpp"

int main()
{
    const cusfloat r0 = static_cast<cusfloat>(0.5);
    const cusfloat r1 = static_cast<cusfloat>(2.0);
    const int theta_np = 64;

    const cusfloat area = PI * (r1 * r1 - r0 * r0);

    auto run_check = [&](AngularRule rule, const char* label, cusfloat rel_tol) -> bool {
        auto cloud = create_annulus_nodes<6>(r0, r1, theta_np, rule);

        cusfloat sum = 0.0;
        for (int i = 0; i < cloud.nodes_np; ++i)
        {
            sum += cloud.w[static_cast<size_t>(i)];
        }

        const cusfloat abs_diff = std::abs(sum - area);
        const cusfloat rel_diff = abs_diff / std::max(area, static_cast<cusfloat>(1e-12));

        if (rel_diff > rel_tol)
        {
            std::cerr << "Annulus quadrature test failed (" << label << ").\n";
            std::cerr << "Expected area: " << area << "\n";
            std::cerr << "Computed sum: " << sum << "\n";
            std::cerr << "Abs diff    : " << abs_diff << "\n";
            std::cerr << "Rel diff    : " << rel_diff << "\n";
            std::cerr << "Rel tol     : " << rel_tol << "\n";
            return false;
        }

        std::cout << "Annulus quadrature test passed (" << label << ").\n";
        std::cout << "Area        : " << area << "\n";
        std::cout << "Sum weights : " << sum << "\n";
        return true;
    };

    const bool ok_trap = run_check(AngularRule::Trapezoidal, "trapezoidal", static_cast<cusfloat>(1e-4));
    const bool ok_cheb = run_check(AngularRule::GaussChebyshev, "gauss-chebyshev", static_cast<cusfloat>(5e-4));

    // Additional annulus integration checks with analytic references.
    // The reference values correspond to integrals in (r, theta) without the
    // polar Jacobian r, so we divide by r when using polar-area weights.
    auto run_integral_check = [&](
        const char* label,
        cusfloat r0_in,
        cusfloat r1_in,
        int theta_np_in,
        cusfloat reference,
        cusfloat rel_tol,
        const std::function<cusfloat(cusfloat, cusfloat)>& integrand
    ) -> bool {
        auto cloud = create_annulus_nodes<8>(r0_in, r1_in, theta_np_in, AngularRule::Trapezoidal);

        cusfloat sum = 0.0;
        for (int i = 0; i < cloud.nodes_np; ++i)
        {
            const cusfloat x = cloud.x[static_cast<size_t>(i)];
            const cusfloat y = cloud.y[static_cast<size_t>(i)];
            const cusfloat r = std::sqrt(x * x + y * y);
            cusfloat theta = std::atan2(y, x);
            if (theta < static_cast<cusfloat>(0.0))
            {
                theta += static_cast<cusfloat>(2.0) * PI;
            }
            sum += cloud.w[static_cast<size_t>(i)] * (integrand(r, theta));
        }

        const cusfloat abs_diff = std::abs(sum - reference);
        const cusfloat rel_diff = abs_diff / std::max(reference, static_cast<cusfloat>(1e-12));

        if (rel_diff > rel_tol)
        {
            std::cerr << "Annulus integral test failed (" << label << ").\n";
            std::cerr << "Expected   : " << reference << "\n";
            std::cerr << "Computed   : " << sum << "\n";
            std::cerr << "Abs diff   : " << abs_diff << "\n";
            std::cerr << "Rel diff   : " << rel_diff << "\n";
            std::cerr << "Rel tol    : " << rel_tol << "\n";
            return false;
        }

        std::cout << "Annulus integral test passed (" << label << ").\n";
        std::cout << "Reference : " << reference << "\n";
        std::cout << "Computed  : " << sum << "\n";
        return true;
    };

    const cusfloat ref1 = static_cast<cusfloat>(0.576884156739537);
    const cusfloat ref2 = static_cast<cusfloat>(7.409890356169836);
    const cusfloat ref3 = static_cast<cusfloat>(0.453887272039892);

    const bool ok_int_1 = run_integral_check(
        "exp(-r^2) * (1 + cos(3*theta))",
        static_cast<cusfloat>(1.3),
        static_cast<cusfloat>(2.65),
        256,
        ref1,
        static_cast<cusfloat>(8e-4),
        [](cusfloat r, cusfloat theta) {
            return std::exp(-r * r) * (static_cast<cusfloat>(1.0) + std::cos(static_cast<cusfloat>(3.0) * theta));
        }
    );

    const bool ok_int_2 = run_integral_check(
        "exp(-r) * (1 + theta / (2*pi))",
        static_cast<cusfloat>(0.75),
        static_cast<cusfloat>(5.0),
        256,
        ref2,
        static_cast<cusfloat>(1e-2),
        [](cusfloat r, cusfloat theta) {
            return std::exp(-r) * (static_cast<cusfloat>(1.0) + theta / (static_cast<cusfloat>(2.0) * PI));
        }
    );

    const bool ok_int_3 = run_integral_check(
        "exp(-r) * abs(sin(theta))",
        static_cast<cusfloat>(3.66),
        static_cast<cusfloat>(7.142),
        256,
        ref3,
        static_cast<cusfloat>(1e-3),
        [](cusfloat r, cusfloat theta) {
            return std::exp(-r) * std::abs(std::sin(theta));
        }
    );

    return (ok_trap && ok_cheb && ok_int_1 && ok_int_2 && ok_int_3) ? 0 : 1;
}
