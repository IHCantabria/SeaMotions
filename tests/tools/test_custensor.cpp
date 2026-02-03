// Include general usage libraries
#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>

// Include local modules
#include "../../src/config.hpp"
#include "../../src/math/custensor/custensor.hpp"

using namespace cut;

static inline bool almost_equal(cusfloat a, cusfloat b, cusfloat tol = EPS_PRECISION * 10)
{
    return std::abs(a - b) <= tol;
}

int main()
{
    // 1D constructors and element access
    size_t n = 3;
    CusTensor<cusfloat> a(n);
    CusTensor<cusfloat> b(n);
    a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
    b[0] = 4.0; b[1] = 8.0; b[2] = 12.0;

    assert(a.size() == n);
    assert(b.shape().size() == 1 && b.shape()[0] == n);

    // Binary ops and assignment from expressions
    CusTensor<cusfloat> sum = a + b;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(sum[i], a[i] + b[i]));

    CusTensor<cusfloat> diff = a - b;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(diff[i], a[i] - b[i]));

    CusTensor<cusfloat> prod = a * b;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(prod[i], a[i] * b[i]));

    CusTensor<cusfloat> quot = b / a;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(quot[i], b[i] / a[i]));

    // Unary exp and chained expressions
    CusTensor<cusfloat> expb = nda_exp(b);
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(expb[i], std::exp(b[i])));

    CusTensor<cusfloat> exp_ab = nda_exp(a + b);
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(exp_ab[i], std::exp(a[i] + b[i])));

    CusTensor<cusfloat> a_plus_expb = a + nda_exp(b);
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(a_plus_expb[i], a[i] + std::exp(b[i])));

    CusTensor<cusfloat> expb_plus_a = nda_exp(b) + a;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(expb_plus_a[i], std::exp(b[i]) + a[i]));

    // Scalar interactions
    cusfloat s = static_cast<cusfloat>(1.375);
    CusTensor<cusfloat> a_minus_s = a - s;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(a_minus_s[i], a[i] - s));

    CusTensor<cusfloat> s_plus_a = s + a;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(s_plus_a[i], s + a[i]));

    CusTensor<cusfloat> expa_times_s = nda_exp(a) * s;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(expa_times_s[i], std::exp(a[i]) * s));

    CusTensor<cusfloat> s_times_expa = s * nda_exp(a);
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(s_times_expa[i], s * std::exp(a[i])));

    // Longer chain
    CusTensor<cusfloat> chain = a * b + a * b;
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(chain[i], 2.0 * (a[i] * b[i])));

    CusTensor<cusfloat> exp_chain = nda_exp(a) * nda_exp(s) * nda_exp(b);
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(exp_chain[i], std::exp(a[i]) * std::exp(s) * std::exp(b[i])));

    // 2D constructor and operator()
    CusTensor<cusfloat> m({2, 2});
    CusTensor<cusfloat> n2({2, 2});
    m(0, 0) = 0.0; m(0, 1) = 1.0; m(1, 0) = 2.0; m(1, 1) = 3.0;
    n2(0, 0) = 0.0; n2(0, 1) = 2.0; n2(1, 0) = 4.0; n2(1, 1) = 6.0;
    CusTensor<cusfloat> r = m + n2;
    assert(r.shape().size() == 2 && r.shape()[0] == 2 && r.shape()[1] == 2);
    assert(almost_equal(r[0], 0.0));
    assert(almost_equal(r[1], 3.0));
    assert(almost_equal(r[2], 6.0));
    assert(almost_equal(r[3], 9.0));

    // Copy constructor deep-copies
    CusTensor<cusfloat> a_copy(a);
    for (size_t i = 0; i < n; ++i)
        assert(almost_equal(a_copy[i], a[i]));
    a[0] += 10.0;
    assert(!almost_equal(a_copy[0], a[0])); // independent storage

    // Move constructor transfers ownership
    CusTensor<cusfloat> a_moved(std::move(a));
    assert(a.size() == 0);
    for (size_t i = 1; i < n; ++i)
        assert(almost_equal(a_moved[i], (i == 0 ? a_moved[i] : (i + 1)))) ; // sanity: moved data present

    // Move assignment transfers ownership
    CusTensor<cusfloat> tmp(n);
    tmp[0] = 7.0; tmp[1] = 8.0; tmp[2] = 9.0;
    CusTensor<cusfloat> target;
    target = std::move(tmp);
    assert(tmp.size() == 0);
    assert(almost_equal(target[0], 7.0));
    assert(almost_equal(target[1], 8.0));
    assert(almost_equal(target[2], 9.0));

    // std::vector emplace_back should move without aliasing
    std::vector<CusTensor<cusfloat>> vec;
    for (int k = 0; k < 20; ++k)
    {
        CusTensor<cusfloat> t(5);
        for (int i = 0; i < 5; ++i) t[i] = static_cast<cusfloat>(k * 10 + i);
        vec.emplace_back(std::move(t));
        assert(t.size() == 0);
    }
    for (int k = 0; k < 20; ++k)
    {
        for (int i = 0; i < 5; ++i)
        {
            assert(almost_equal(vec[k][i], static_cast<cusfloat>(k * 10 + i)));
        }
        assert(vec[k].shape().size() == 1 && vec[k].shape()[0] == 5);
    }

    // Holder object with CusTensor field, emplaced_back to trigger moves on reallocation
    struct TensorHolder
    {
        cut::CusTensor<cusfloat> data;
        int id = -1;

        TensorHolder(size_t n, int id_in) : data(n), id(id_in)
        {
            for (size_t i = 0; i < n; ++i)
            {
                data[i] = static_cast<cusfloat>(id * 100 + static_cast<int>(i));
            }
        }
    };

    std::vector<TensorHolder> holders;
    holders.reserve(2); // small reserve to encourage multiple reallocations
    for (int k = 0; k < 50; ++k)
    {
        holders.emplace_back(static_cast<size_t>(5), k);
    }

    for (int k = 0; k < 50; ++k)
    {
        assert(holders[k].data.shape().size() == 1 && holders[k].data.shape()[0] == 5);
        for (int i = 0; i < 5; ++i)
        {
            cusfloat expected = static_cast<cusfloat>(k * 100 + i);
            assert(almost_equal(holders[k].data[i], expected));
        }
    }

    // Stress: large 1D holders to trigger heavier moves
    struct Large1DHolder
    {
        cut::CusTensor<cusfloat> data;
        int id;
        Large1DHolder(size_t n, int id_in) : data(n), id(id_in)
        {
            for (size_t i = 0; i < n; ++i)
            {
                data[i] = static_cast<cusfloat>(id * 1000 + static_cast<int>(i));
            }
        }
    };

    std::vector<Large1DHolder> bigs;
    bigs.reserve(2);
    for (int k = 0; k < 40; ++k)
    {
        bigs.emplace_back(static_cast<size_t>(128), k);
    }
    for (int k = 0; k < 40; ++k)
    {
        assert(bigs[k].data.shape().size() == 1 && bigs[k].data.shape()[0] == 128);
        assert(almost_equal(bigs[k].data[0], static_cast<cusfloat>(k * 1000 + 0)));
        assert(almost_equal(bigs[k].data[127], static_cast<cusfloat>(k * 1000 + 127)));
    }

    // Stress: 2D holders with operator() access
    struct Tensor2DHolder
    {
        cut::CusTensor<cusfloat> data;
        size_t rows;
        size_t cols;
        int id;
        Tensor2DHolder(size_t r, size_t c, int id_in) : data({r, c}), rows(r), cols(c), id(id_in)
        {
            for (size_t i = 0; i < r; ++i)
            {
                for (size_t j = 0; j < c; ++j)
                {
                    data(i, j) = static_cast<cusfloat>(id * 100 + static_cast<int>(i * 10 + j));
                }
            }
        }
    };

    std::vector<Tensor2DHolder> mats;
    mats.reserve(2);
    for (int k = 0; k < 30; ++k)
    {
        mats.emplace_back(static_cast<size_t>(10), static_cast<size_t>(6), k);
    }
    for (int k = 0; k < 30; ++k)
    {
        auto shp = mats[k].data.shape();
        assert(shp.size() == 2 && shp[0] == 10 && shp[1] == 6);
        assert(almost_equal(mats[k].data(0, 0), static_cast<cusfloat>(k * 100 + 0)));
        assert(almost_equal(mats[k].data(9, 5), static_cast<cusfloat>(k * 100 + 95)));
    }

    // Stress: varying shape holders to check shape/values preserved across moves
    struct VaryHolder
    {
        cut::CusTensor<cusfloat> data;
        size_t rows;
        size_t cols;
        int id;
        VaryHolder(size_t r, size_t c, int id_in) : data({r, c}), rows(r), cols(c), id(id_in)
        {
            for (size_t i = 0; i < r; ++i)
            {
                for (size_t j = 0; j < c; ++j)
                {
                    data(i, j) = static_cast<cusfloat>(id * 1000 + static_cast<int>(i * 100 + j));
                }
            }
        }
    };

    std::vector<VaryHolder> varies;
    varies.reserve(1);
    for (int k = 0; k < 25; ++k)
    {
        size_t r = 2 + (k % 3);
        size_t c = 3 + (k % 4);
        varies.emplace_back(r, c, k);
    }
    for (int k = 0; k < 25; ++k)
    {
        auto shp = varies[k].data.shape();
        size_t r = 2 + (k % 3);
        size_t c = 3 + (k % 4);
        assert(shp.size() == 2 && shp[0] == r && shp[1] == c);
        assert(almost_equal(varies[k].data(0, 0), static_cast<cusfloat>(k * 1000 + 0)));
        assert(almost_equal(varies[k].data(r - 1, c - 1), static_cast<cusfloat>(k * 1000 + static_cast<int>((r - 1) * 100 + (c - 1)))));
    }

    // Equality operator checks for different data types
    {
        // Integral type
        cut::CusTensor<int> ti1(static_cast<size_t>(5));
        cut::CusTensor<int> ti2(static_cast<size_t>(5));
        for (size_t i = 0; i < 5; ++i) { ti1[i] = static_cast<int>(i); ti2[i] = static_cast<int>(i); }
        assert(ti1 == ti2);
        ti2[0] += 1;
        assert(!(ti1 == ti2));

        // Shape mismatch should be false
        cut::CusTensor<int> ti3({2, 3});
        for (size_t i = 0; i < 6; ++i) { ti3[i] = static_cast<int>(i); }
        assert(!(ti1 == ti3));
    }

    {
        // Floating type (cusfloat)
        cut::CusTensor<cusfloat> tf1(static_cast<size_t>(5));
        cut::CusTensor<cusfloat> tf2(static_cast<size_t>(5));
        for (size_t i = 0; i < 5; ++i) { tf1[i] = static_cast<cusfloat>(i * 0.5); tf2[i] = static_cast<cusfloat>(i * 0.5); }
        assert(tf1 == tf2);
        // Change larger than tolerance should make them unequal
        tf2[0] += static_cast<cusfloat>(EPS_PRECISION * 2.0);
        assert(!(tf1 == tf2));
    }

    {
        // Complex type (cuscomplex)
        cut::CusTensor<cuscomplex> tc1(static_cast<size_t>(4));
        cut::CusTensor<cuscomplex> tc2(static_cast<size_t>(4));
        for (size_t i = 0; i < 4; ++i)
        {
            cusfloat re = static_cast<cusfloat>(i * 0.1);
            cusfloat im = static_cast<cusfloat>(i * 0.2);
            tc1[i] = cuscomplex(re, im);
            tc2[i] = cuscomplex(re, im);
        }
        assert(tc1 == tc2);
        // Change larger than tolerance should make them unequal
        tc2[1] += cuscomplex(static_cast<cusfloat>(EPS_PRECISION * 2.0), static_cast<cusfloat>(0));
        assert(!(tc1 == tc2));
    }

    std::cout << "test_custensor: OK" << std::endl;
    return 0;
}
