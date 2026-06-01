"""
Test: compare dGdt_asymptotic against tabulated values in Gt.h5.

The asymptotic expression is valid for large beta; the comparison is therefore
restricted to beta >= BETA_MIN (default 17).

Tabulated file layout  (Gt.h5)
-------------------------------
  beta   : 1-D array, shape (N_beta,)   – non-dimensional time parameter
  mu     : 1-D array, shape (N_mu,)     – non-dimensional frequency parameter
  fcn    : 2-D array, shape (N_beta, N_mu) – reference dG/dt values

Usage
-----
    python test_dGdt_asymptotic.py                       # compare dGdt against Gt.h5
    python test_dGdt_asymptotic.py --plot                # also generate error plots
    python test_dGdt_asymptotic.py --compare-deriv-tables  # compare dGdtt/dGdttt/dGdtx/dGdttx against tabulated HDF5
    python test_dGdt_asymptotic.py --skip-table --test-deriv  # finite-difference derivative validation only
    python test_dGdt_asymptotic.py --beta-min 30 --compare-deriv-tables  # restrict beta range
"""

import os
import sys
import argparse
import numpy as np
import h5py

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT  = os.path.dirname(_TOOLS_DIR)

HDF5_PATH = os.path.join(
    _REPO_ROOT, 'aux_data', '0_integrals_database', '1_time_domain', 'Gt.h5'
)

_TD_DIR = os.path.join(_REPO_ROOT, 'aux_data', '0_integrals_database', '1_time_domain')

HDF5_DERIV_PATHS = {
    'dGdtt' : os.path.join(_TD_DIR, 'Gtt.h5'),
    'dGdttt': os.path.join(_TD_DIR, 'Gttt.h5'),
    'dGdtx' : os.path.join(_TD_DIR, 'Gtx.h5'),
    'dGdttx': os.path.join(_TD_DIR, 'Gttx.h5'),
}

# Add aux_tools to sys.path so the module is importable regardless of cwd
if _TOOLS_DIR not in sys.path:
    sys.path.insert(0, _TOOLS_DIR)

from time_domain_asymptotic import (          # noqa: E402
    dGdt_asymptotic,
    dGdtt_asymptotic,
    dGdttt_asymptotic,
    dGdtx_asymptotic,
    dGdttx_asymptotic,
    get_alpha,
    get_dmn,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BETA_MIN = 17.0   # asymptotic expression is valid for large beta


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _compute_asymptotic_grid(beta_arr, mu_arr):
    """
    Evaluate dGdt_asymptotic on the full (beta, mu) grid.

    Returns
    -------
    result : complex ndarray, shape (len(beta_arr), len(mu_arr))
    """
    n_beta = len(beta_arr)
    n_mu   = len(mu_arr)
    result = np.empty((n_beta, n_mu), dtype=float)

    for i, b in enumerate(beta_arr):
        if (i % 500) == 0:
            pct = 100.0 * i / n_beta
            print(f"  evaluating asymptotic: {pct:5.1f}%  (beta = {b:.2f})", flush=True)
        # dGdt_asymptotic is vectorisable over mu when beta is scalar
        result[i, :] = dGdt_asymptotic(b, mu_arr)

    print("  evaluating asymptotic: 100.0%  done", flush=True)
    return result


# ---------------------------------------------------------------------------
# Main comparison
# ---------------------------------------------------------------------------

def compare_dGdt_asymptotic(beta_min=BETA_MIN):
    """
    Load the reference table, evaluate the asymptotic formula for beta >= beta_min
    and return error statistics.

    Returns
    -------
    beta_sub  : ndarray (N_beta,)
    mu        : ndarray (N_mu,)
    fcn_sub   : ndarray (N_beta, N_mu) – reference values
    computed  : ndarray (N_beta, N_mu) – real part of asymptotic result
    abs_err   : ndarray (N_beta, N_mu) – |computed - reference|
    rel_err   : ndarray (N_beta, N_mu) – |computed - reference| / max(|reference|, tol)
    """
    # ------------------------------------------------------------------
    # Load reference data
    # ------------------------------------------------------------------
    print(f"Loading {HDF5_PATH} …")
    with h5py.File(HDF5_PATH, 'r') as fh:
        beta_all = fh['beta'][()]
        mu       = fh['mu'][()]
        fcn_all  = 2.0 * fh['fcn'][()]

    # ------------------------------------------------------------------
    # Restrict to beta >= beta_min
    # ------------------------------------------------------------------
    idx_start = np.searchsorted(beta_all, beta_min)
    beta_sub  = beta_all[idx_start:]
    fcn_sub   = fcn_all[idx_start:, :]

    print(f"beta range : {beta_sub[0]:.4f} … {beta_sub[-1]:.4f}  ({len(beta_sub)} points)")
    print(f"mu  range  : {mu[0]:.6f} … {mu[-1]:.6f}  ({len(mu)} points)")
    print(f"Grid size  : {len(beta_sub)} x {len(mu)} = {len(beta_sub)*len(mu):,} evaluations")

    # ------------------------------------------------------------------
    # Evaluate asymptotic expression
    # ------------------------------------------------------------------
    asym_complex = _compute_asymptotic_grid(beta_sub, mu)
    computed     = np.real(asym_complex)   # physical value must be real

    # import matplotlib.pyplot as plt
    # plt.plot( beta_sub, asym_complex[:, 0] )
    # plt.plot( beta_sub, fcn_sub[:, 0] )
    # plt.xlabel('beta')
    # plt.ylabel('dG/dt')
    # plt.title('Asymptotic Approximation')
    # plt.show()

    # ------------------------------------------------------------------
    # Errors
    # ------------------------------------------------------------------
    abs_err  = np.abs(computed - fcn_sub)
    tol_denom = np.where(np.abs(fcn_sub) > 1e-10, np.abs(fcn_sub), 1e-10)
    rel_err  = abs_err / tol_denom

    return beta_sub, mu, fcn_sub, computed, abs_err, rel_err


def print_statistics(beta_sub, mu, fcn_sub, computed, abs_err, rel_err,
                     label='dG/dt vs. Gt.h5'):
    """Print a summary of the comparison statistics."""
    print()
    print("=" * 60)
    print(f"  Comparison statistics  ({label})")
    print("=" * 60)

    def _row(label, val):
        print(f"  {label:<40s} {val:.4e}")

    _row("Max  |absolute error|",   abs_err.max())
    _row("Mean |absolute error|",   abs_err.mean())
    _row("Max  |relative error|",   rel_err.max())
    _row("Mean |relative error|",   rel_err.mean())

    print()

    # Worst-case points
    flat_abs  = abs_err.ravel()
    flat_rel  = rel_err.ravel()
    i_abs_max = np.argmax(flat_abs)
    i_rel_max = np.argmax(flat_rel)
    ib_abs, im_abs = np.unravel_index(i_abs_max, abs_err.shape)
    ib_rel, im_rel = np.unravel_index(i_rel_max, rel_err.shape)

    print(f"  Worst absolute error at  beta={beta_sub[ib_abs]:.4f}, mu={mu[im_abs]:.6f}")
    print(f"    reference  = {fcn_sub[ib_abs, im_abs]:.6e}")
    print(f"    asymptotic = {computed[ib_abs, im_abs]:.6e}")
    print(f"    |error|    = {abs_err[ib_abs, im_abs]:.6e}")

    print()
    print(f"  Worst relative error at  beta={beta_sub[ib_rel]:.4f}, mu={mu[im_rel]:.6f}")
    print(f"    reference  = {fcn_sub[ib_rel, im_rel]:.6e}")
    print(f"    asymptotic = {computed[ib_rel, im_rel]:.6e}")
    print(f"    |error|    = {abs_err[ib_rel, im_rel]:.6e}")
    print(f"    rel error  = {rel_err[ib_rel, im_rel]:.6e}")
    print("=" * 60)

    # Per-mu statistics (averaged over beta)
    print()
    print("  Per-mu mean absolute error (top 5 worst mu values):")
    mean_abs_per_mu = abs_err.mean(axis=0)
    worst_mu_idx = np.argsort(mean_abs_per_mu)[::-1][:5]
    for idx in worst_mu_idx:
        print(f"    mu = {mu[idx]:.6f}  mean|abs err| = {mean_abs_per_mu[idx]:.4e}")

    print()
    print("  Per-beta mean absolute error  (sampled every 100 beta steps):")
    mean_abs_per_beta = abs_err.mean(axis=1)
    for i in range(0, len(beta_sub), 100):
        print(f"    beta = {beta_sub[i]:.2f}  mean|abs err| = {mean_abs_per_beta[i]:.4e}")


# ---------------------------------------------------------------------------
# Generic table comparison for derivative functions
# ---------------------------------------------------------------------------

def _compare_derivative_table(h5_path, asym_fn, label, beta_min=BETA_MIN):
    """
    Load a tabulated HDF5 file, evaluate *asym_fn* on the same grid for
    beta >= beta_min, and return error statistics.

    The tabulated values are scaled by 2.0 before comparison (same convention
    as compare_dGdt_asymptotic / Gt.h5).

    Returns
    -------
    beta_sub, mu, fcn_sub, computed, abs_err, rel_err
        Same layout as the return value of compare_dGdt_asymptotic.
    """
    print(f"Loading {h5_path} …")
    with h5py.File(h5_path, 'r') as fh:
        beta_all = fh['beta'][()]
        mu       = fh['mu'][()]
        fcn_all  = 2.0 * fh['fcn'][()]

    idx_start = np.searchsorted(beta_all, beta_min)
    beta_sub  = beta_all[idx_start:]
    fcn_sub   = fcn_all[idx_start:, :]

    print(f"beta range : {beta_sub[0]:.4f} … {beta_sub[-1]:.4f}  ({len(beta_sub)} points)")
    print(f"mu  range  : {mu[0]:.6f} … {mu[-1]:.6f}  ({len(mu)} points)")
    print(f"Grid size  : {len(beta_sub)} x {len(mu)} = {len(beta_sub)*len(mu):,} evaluations")

    n_beta   = len(beta_sub)
    computed = np.empty((n_beta, len(mu)), dtype=float)
    for i, b in enumerate(beta_sub):
        if (i % 500) == 0:
            pct = 100.0 * i / n_beta
            print(f"  evaluating {label}: {pct:5.1f}%  (beta = {b:.2f})", flush=True)
        computed[i, :] = asym_fn(b, mu)
    print(f"  evaluating {label}: 100.0%  done", flush=True)

    abs_err   = np.abs(computed - fcn_sub)
    tol_denom = np.where(np.abs(fcn_sub) > 1e-10, np.abs(fcn_sub), 1e-10)
    rel_err   = abs_err / tol_denom

    return beta_sub, mu, fcn_sub, computed, abs_err, rel_err


def compare_dGdtt_asymptotic(beta_min=BETA_MIN):
    """Compare dGdtt_asymptotic against Gtt.h5."""
    return _compare_derivative_table(
        HDF5_DERIV_PATHS['dGdtt'], dGdtt_asymptotic, 'dGdtt', beta_min
    )


def compare_dGdttt_asymptotic(beta_min=BETA_MIN):
    """Compare dGdttt_asymptotic against Gttt.h5."""
    return _compare_derivative_table(
        HDF5_DERIV_PATHS['dGdttt'], dGdttt_asymptotic, 'dGdttt', beta_min
    )


def compare_dGdtx_asymptotic(beta_min=BETA_MIN):
    """Compare dGdtx_asymptotic against Gtx.h5."""
    return _compare_derivative_table(
        HDF5_DERIV_PATHS['dGdtx'], dGdtx_asymptotic, 'dGdtx', beta_min
    )


def compare_dGdttx_asymptotic(beta_min=BETA_MIN):
    """Compare dGdttx_asymptotic against Gttx.h5."""
    return _compare_derivative_table(
        HDF5_DERIV_PATHS['dGdttx'], dGdttx_asymptotic, 'dGdttx', beta_min
    )


def compare_all_derivative_tables(beta_min=BETA_MIN, plot=False):
    """
    Run table comparisons for all four derivative functions and print statistics.
    When plot=True, generate contourf error maps for all cases in one figure.
    """
    entries = [
        ('dGdtt',  compare_dGdtt_asymptotic),
        ('dGdttt', compare_dGdttt_asymptotic),
        ('dGdtx',  compare_dGdtx_asymptotic),
        ('dGdttx', compare_dGdttx_asymptotic),
    ]
    plot_data = []
    for label, fn in entries:
        result = fn(beta_min)
        print_statistics(*result, label=f'{label} vs. {os.path.basename(HDF5_DERIV_PATHS[label])}')
        if plot:
            # (label, beta_sub, mu, abs_err, rel_err)
            plot_data.append((label, result[0], result[1], result[4], result[5]))
    if plot:
        make_derivative_plots(plot_data)


# ---------------------------------------------------------------------------
# High-beta simplification tests: f0 negligibility and nmax convergence
# ---------------------------------------------------------------------------

def test_high_beta_simplifications(beta_min=55.0, beta_max=60.0,
                                    nmax_full=15, nmax_reduced=6,
                                    tol_f0=1e-3, tol_nmax=1e-8):
    """
    Check two simplifications that hold for beta in [beta_min, beta_max],
    applied to ALL five asymptotic derivative functions:
    dGdt, dGdtt, dGdttt, dGdtx, dGdttx.

    For each function the expansion is   g = f0 + f2   where:
      f0  — power series in 1/β (non-oscillatory correction)
      f2  — stationary-phase / oscillatory contribution

    Test 1 – f0 negligibility
    --------------------------
    Assert |f0| / (|f0| + |f2|) < tol_f0 at representative (β, μ) points
    where f2 is not exponentially suppressed (small μ, so exp(-β²μ/4) ~ 1).

    Test 2 – nmax convergence
    --------------------------
    Assert that reducing the inner-series truncation from nmax_full to
    nmax_reduced changes the f2 value by < tol_nmax (relative).  A
    convergence table across nmax = 2, 4, nmax_reduced, 8, nmax_full is
    printed for inspection.

    Parameters
    ----------
    beta_min, beta_max : float   Beta range to sample.  Default [55, 60].
    nmax_full : int              Reference nmax (default 15).
    nmax_reduced : int           Reduced nmax to assert sufficient (default 6).
    tol_f0 : float              Max |f0|/(|f0|+|f2|) ratio  (default 1e-3).
    tol_nmax : float            Max relative change when lowering nmax (default 1e-8).
    """
    import scipy.special as sp

    j = complex(0.0, 1.0)

    # ------------------------------------------------------------------
    # Per-derivative helpers: each returns (f0_part, f2_part) as scalars.
    # They mirror the internal structure of the production functions so
    # the split between f0 and f2 is exact.
    # ------------------------------------------------------------------

    def _f2a(beta, mu):
        """Shared oscillatory envelope (common to all five functions)."""
        eta = np.sqrt(1.0 - mu**2)
        return (
            complex(0.0, -4.0)
            * np.exp(-0.25 * beta**2 * mu)
            * np.exp(j * (0.25*beta**2*eta - 0.5*np.arccos(mu) + np.pi/4))
            / np.sqrt(2.0 * eta)
        )

    def _f0_f2_dGdt(beta, mu, nmax):
        eta = np.sqrt(1.0 - mu**2)
        z1  = j*mu + eta;  z2 = mu - j*eta
        f0  = sum(get_alpha(n) * beta**(-2*n-3) * sp.eval_legendre(n, mu)
                  for n in range(nmax)) * (-4.0)
        fa  = _f2a(beta, mu)
        f2b = 0.0
        for n in range(nmax):
            r0 = (j/eta)**n;  r1 = 0.0
            for m in range(n+1):
                p = 1.0 - 2.0*m - 2.0*n
                r1 += get_dmn(m, n) * (0.5*beta*z1)**p * z2**m
            f2b += r0 * r1
        return f0, (fa * f2b).real

    def _f0_f2_dGdtt(beta, mu, nmax):
        eta = np.sqrt(1.0 - mu**2)
        z1  = j*mu + eta;  z2 = mu - j*eta
        df0 = sum(get_alpha(n) * (2*n+3) * beta**(-2*n-4) * sp.eval_legendre(n, mu)
                  for n in range(nmax)) * 4.0
        fa    = _f2a(beta, mu)
        dfa   = fa * (0.5 * beta * (-mu + j*eta))
        f2b   = 0.0;  df2b = 0.0
        for n in range(nmax):
            r0 = (j/eta)**n;  r1 = 0.0;  dr1 = 0.0
            for m in range(n+1):
                p    = 1.0 - 2.0*m - 2.0*n
                base = get_dmn(m, n) * (0.5*beta*z1)**p * z2**m
                r1  += base;  dr1 += p * base
            f2b  += r0 * r1;  df2b += r0 * dr1
        df2b /= beta
        return df0, (dfa * f2b + fa * df2b).real

    def _f0_f2_dGdttt(beta, mu, nmax):
        eta = np.sqrt(1.0 - mu**2)
        z1  = j*mu + eta;  z2 = mu - j*eta
        d2f0 = sum(get_alpha(n) * (2*n+3) * (2*n+4) * beta**(-2*n-5) * sp.eval_legendre(n, mu)
                   for n in range(nmax)) * (-4.0)
        fa    = _f2a(beta, mu)
        hp    = 0.5 * beta * (-mu + j*eta)
        hpp   = 0.5 * (-mu + j*eta)
        dfa   = fa * hp
        d2fa  = fa * (hp**2 + hpp)
        S0 = S1 = S2 = 0.0
        for n in range(nmax):
            r0 = (j/eta)**n;  s0 = s1 = s2 = 0.0
            for m in range(n+1):
                p    = 1.0 - 2.0*m - 2.0*n
                base = get_dmn(m, n) * (0.5*beta*z1)**p * z2**m
                s0 += base;  s1 += p*base;  s2 += p*(p-1.0)*base
            S0 += r0*s0;  S1 += r0*s1;  S2 += r0*s2
        d2f2 = d2fa*S0 + 2.0*dfa*(S1/beta) + fa*(S2/beta**2)
        return d2f0, d2f2.real

    def _dPn(n, mu, eta):
        """P_n'(μ) via standard recurrence; requires n ≥ 1."""
        return n * (sp.eval_legendre(n-1, mu) - mu*sp.eval_legendre(n, mu)) / eta**2

    def _f0_f2_dGdtx(beta, mu, nmax):
        eta = np.sqrt(1.0 - mu**2)
        z1  = j*mu + eta;  z2 = mu - j*eta
        df0_dmu = sum(get_alpha(n) * beta**(-2*n-3) * _dPn(n, mu, eta)
                      for n in range(1, nmax)) * (-4.0)
        fa    = _f2a(beta, mu)
        c_mu  = mu/(2.0*eta**2) - 0.25*beta**2 + j*(2.0 - beta**2*mu)/(4.0*eta)
        dfa_mu = fa * c_mu
        f2b = 0.0;  df2b_mu = 0.0
        for n in range(nmax):
            r0 = (j/eta)**n;  r1 = 0.0;  dr1 = 0.0
            for m in range(n+1):
                p    = 1.0 - 2.0*m - 2.0*n
                base = get_dmn(m, n) * (0.5*beta*z1)**p * z2**m
                c    = -p * z2 / (eta*z1)
                if m > 0:
                    c += m * z1 / (eta*z2)
                r1  += base;  dr1 += base * c
            f2b     += r0 * r1
            df2b_mu += r0 * (n*mu/eta**2 * r1 + dr1)
        return df0_dmu, (dfa_mu*f2b + fa*df2b_mu).real

    def _f0_f2_dGdttx(beta, mu, nmax):
        eta = np.sqrt(1.0 - mu**2)
        z1  = j*mu + eta;  z2 = mu - j*eta
        d2f0 = sum(get_alpha(n) * (2*n+3) * beta**(-2*n-4) * _dPn(n, mu, eta)
                   for n in range(1, nmax)) * 4.0
        fa    = _f2a(beta, mu)
        hp    = 0.5 * beta * (-mu + j*eta)
        fa_b  = fa * hp
        c_mu  = mu/(2.0*eta**2) - 0.25*beta**2 + j*(2.0 - beta**2*mu)/(4.0*eta)
        fa_mu = fa * c_mu
        hp_mu = -0.5 * beta * z1 / eta
        fa_bmu = fa * (c_mu * hp + hp_mu)
        S0 = dS0 = S1 = dS1 = 0.0
        for n in range(nmax):
            r0 = (j/eta)**n;  r1 = dr1 = s1 = ds1 = 0.0
            for m in range(n+1):
                p    = 1.0 - 2.0*m - 2.0*n
                base = get_dmn(m, n) * (0.5*beta*z1)**p * z2**m
                c    = -p * z2 / (eta*z1)
                if m > 0:
                    c += m * z1 / (eta*z2)
                r1  += base;         dr1 += base*c
                s1  += p*base;       ds1 += p*base*c
            nm = n * mu / eta**2
            S0  += r0 * r1;          dS0 += r0*(nm*r1 + dr1)
            S1  += r0 * s1;          dS1 += r0*(nm*s1 + ds1)
        d2f2 = fa_bmu*S0 + fa_b*dS0 + fa_mu*(S1/beta) + fa*(dS1/beta)
        return d2f0, d2f2.real

    # dict preserves insertion order (Python 3.7+)
    _HELPERS = {
        'dGdt'  : _f0_f2_dGdt,
        'dGdtt' : _f0_f2_dGdtt,
        'dGdttt': _f0_f2_dGdttt,
        'dGdtx' : _f0_f2_dGdtx,
        'dGdttx': _f0_f2_dGdttx,
    }

    # Sample points: small μ so exp(-β²μ/4) ~ 1 and f2 is non-negligible.
    # At β=55: μ_crit ~ 4/55² ≈ 0.0013.  Use a range spanning below & above.
    betas      = np.linspace(beta_min, beta_max, 5)
    mus        = np.array([0.0001, 0.0005, 0.001, 0.003, 0.007])
    nmax_steps = sorted({2, 4, nmax_reduced, 8, nmax_full})
    b_ref      = betas[2]      # middle β for convergence table

    W = 82
    n_fail_global = 0

    for label, fn in _HELPERS.items():

        print()
        print('=' * W)
        print(f'  [{label}]  High-beta simplification  (β ∈ [{beta_min:.1f}, {beta_max:.1f}])')
        print('=' * W)

        # --------------------------------------------------------------
        # Test 1: f0 negligibility
        # --------------------------------------------------------------
        print()
        print(f'  Test 1: |f0| / (|f0|+|f2|) < {tol_f0:.0e}')
        hdr = (f"  {'beta':>7}  {'mu':>8}  {'|f0|':>12}  {'|f2|':>12}"
               f"  {'ratio':>10}  {'status':>6}")
        print(hdr)
        print('  ' + '-' * (len(hdr) - 2))

        n_fail_f0 = 0
        for b in betas:
            for m in mus:
                f0, f2 = fn(b, m, nmax_full)
                ratio  = abs(f0) / (abs(f0) + abs(f2) + 1e-300)
                ok     = ratio < tol_f0
                n_fail_f0 += 0 if ok else 1
                print(f"  {b:7.2f}  {m:8.5f}  {abs(f0):12.4e}  {abs(f2):12.4e}"
                      f"  {ratio:10.4e}  {'PASS' if ok else 'FAIL':>6}")

        # --------------------------------------------------------------
        # Test 2: nmax convergence
        # --------------------------------------------------------------
        print()
        print(f'  Test 2: nmax convergence — '
              f'|f2(nmax) − f2({nmax_full})| / |f2({nmax_full})|')
        print(f'  (sample: β = {b_ref:.2f})')

        col_hdr = ('  ' + f"{'mu':>8}"
                   + ''.join(f"  {'nmax='+str(k):>12}"
                              for k in nmax_steps if k != nmax_full))
        print(col_hdr)
        print('  ' + '-' * (len(col_hdr) - 2))

        for m in mus:
            _, f2_ref = fn(b_ref, m, nmax_full)
            row = f'  {m:8.5f}'
            for nk in nmax_steps:
                if nk == nmax_full:
                    continue
                _, f2_nk = fn(b_ref, m, nk)
                if abs(f2_ref) < 1e-30:
                    row += f'  {"(~0)":>12}'
                else:
                    row += f'  {abs(f2_ref - f2_nk)/abs(f2_ref):12.4e}'
            print(row)

        # assert nmax_reduced is sufficient across all (β, μ) pairs
        print()
        print(f'  Asserting nmax={nmax_reduced} sufficient '
              f'(rel. change < {tol_nmax:.0e}):')
        hdr2 = (f"  {'beta':>7}  {'mu':>8}"
                f"  {'f2(nmax='+str(nmax_full)+')':>18}"
                f"  {'f2(nmax='+str(nmax_reduced)+')':>18}"
                f"  {'rel_diff':>12}  {'status':>6}")
        print(hdr2)
        print('  ' + '-' * (len(hdr2) - 2))

        n_fail_nmax = 0
        for b in betas:
            for m in mus:
                _, f2_full = fn(b, m, nmax_full)
                _, f2_red  = fn(b, m, nmax_reduced)
                if abs(f2_full) < 1e-30:
                    print(f"  {b:7.2f}  {m:8.5f}"
                          f"  {'(~0)':>18}  {'(~0)':>18}"
                          f"  {'--':>12}  {'SKIP':>6}")
                    continue
                rel = abs(f2_full - f2_red) / abs(f2_full)
                ok  = rel < tol_nmax
                n_fail_nmax += 0 if ok else 1
                print(f"  {b:7.2f}  {m:8.5f}"
                      f"  {f2_full:18.10e}  {f2_red:18.10e}"
                      f"  {rel:12.4e}  {'PASS' if ok else 'FAIL':>6}")

        n_local = n_fail_f0 + n_fail_nmax
        n_fail_global += n_local
        print()
        if n_local == 0:
            print(f'  [{label}] PASSED  '
                  f'(f0 ratio < {tol_f0:.0e}, nmax={nmax_reduced} rel.err < {tol_nmax:.0e})')
        else:
            print(f'  [{label}] FAILED  ({n_local} check(s))')

    # ------------------------------------------------------------------
    # Global summary
    # ------------------------------------------------------------------
    print()
    print('=' * W)
    if n_fail_global == 0:
        print(f'  All simplification checks PASSED  '
              f'(all 5 functions, β ∈ [{beta_min:.1f}, {beta_max:.1f}])')
        print(f'  → f0 term can be dropped   '
              f'(|f0|/(|f0|+|f2|) < {tol_f0:.0e} for β ≥ {beta_min:.0f})')
        print(f'  → nmax can be reduced to {nmax_reduced}   '
              f'(rel. change vs nmax={nmax_full} < {tol_nmax:.0e})')
    else:
        print(f'  {n_fail_global} simplification check(s) FAILED')
    print('=' * W)
    return n_fail_global == 0


# ---------------------------------------------------------------------------
# Finite-difference consistency tests for derivative functions
# ---------------------------------------------------------------------------

def test_derivatives_fd(tol_rel=5e-3):
    """
    Validate dGdtt, dGdttt, dGdtx, dGdttx against centred finite differences.

    For each (beta, mu) sample point four checks are performed:

      dGdtt  ≈ (dGdt(β+h, μ) - dGdt(β-h, μ)) / (2h)
      dGdttt ≈ (dGdtt(β+h, μ) - dGdtt(β-h, μ)) / (2h)
      dGdtx  ≈ (dGdt(β, μ+h) - dGdt(β, μ-h)) / (2h)
      dGdttx ≈ (dGdtt(β, μ+h) - dGdtt(β, μ-h)) / (2h)   [primary]
             ≈ (dGdtx(β+h, μ) - dGdtx(β-h, μ)) / (2h)   [cross-check]
    Note: at very large beta and mu close to 1 the function values are
    astronomically small (~1e-130), causing FD in μ to suffer cancellation
    errors that can reach ~3e-3 relative.  The default tolerance is set to
    5e-3 to accommodate this; the cross-check in β always passes cleanly.    A failure is printed and counted; the function returns True if all pass.
    """
    betas  = [17.0, 25.0, 40.0]
    mus    = [0.2,  0.5,  0.8]
    h_beta = 1e-4
    h_mu   = 1e-4

    cases = [
        # (label, analytic_fn, fd_fn, vary, h)
        ('dGdtt',  dGdtt_asymptotic,  lambda b, m: (dGdt_asymptotic(b+h_beta, m) - dGdt_asymptotic(b-h_beta, m)) / (2*h_beta),  'beta', h_beta),
        ('dGdttt', dGdttt_asymptotic, lambda b, m: (dGdtt_asymptotic(b+h_beta, m) - dGdtt_asymptotic(b-h_beta, m)) / (2*h_beta), 'beta', h_beta),
        ('dGdtx',  dGdtx_asymptotic,  lambda b, m: (dGdt_asymptotic(b, m+h_mu) - dGdt_asymptotic(b, m-h_mu)) / (2*h_mu),        'mu',   h_mu),
        ('dGdttx', dGdttx_asymptotic, lambda b, m: (dGdtt_asymptotic(b, m+h_mu) - dGdtt_asymptotic(b, m-h_mu)) / (2*h_mu),      'mu',   h_mu),
        ('dGdttx(cross)', dGdttx_asymptotic, lambda b, m: (dGdtx_asymptotic(b+h_beta, m) - dGdtx_asymptotic(b-h_beta, m)) / (2*h_beta), 'beta', h_beta),
    ]

    header = f"{'function':<18} {'beta':>6} {'mu':>5} {'analytic':>15} {'FD':>15} {'rel_err':>10} {'status':>6}"
    print()
    print("=" * len(header))
    print("  Finite-difference derivative validation")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    n_fail = 0
    for label, fn_an, fn_fd, vary, h in cases:
        for b in betas:
            for m in mus:
                an  = fn_an(b, m)
                fd  = fn_fd(b, m)
                rel = abs(an - fd) / (abs(fd) + 1e-300)
                ok  = rel < tol_rel
                n_fail += 0 if ok else 1
                status = 'PASS' if ok else 'FAIL'
                print(f"  {label:<18} {b:6.1f} {m:5.2f} {an:15.6e} {fd:15.6e} {rel:10.2e} {status:>6}")
        print()

    print("=" * len(header))
    if n_fail == 0:
        print(f"  All checks PASSED  (tol_rel = {tol_rel:.0e})")
    else:
        print(f"  {n_fail} check(s) FAILED  (tol_rel = {tol_rel:.0e})")
    print("=" * len(header))

    return n_fail == 0


def make_derivative_plots(all_results):
    """
    Generate contourf error maps for all derivative table comparisons in one figure.

    Parameters
    ----------
    all_results : list of (label, beta_sub, mu, abs_err, rel_err)
        One entry per derivative function, as collected by compare_all_derivative_tables.

    Layout
    ------
    N rows × 2 columns, where N = len(all_results).
    Left column  : log10(|absolute error|)
    Right column : log10(|relative error|)
    Both axes use μ on the x-axis and β on the y-axis.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("matplotlib not available – skipping plots.")
        return

    n_cases = len(all_results)
    fig, axes = plt.subplots(
        n_cases, 2,
        figsize=(14, 4 * n_cases),
        squeeze=False,
    )
    fig.suptitle('Asymptotic derivative errors vs. tabulated reference', fontsize=13)

    for row, (label, beta_sub, mu, abs_err, rel_err) in enumerate(all_results):
        MU, BETA = np.meshgrid(mu, beta_sub)

        for col, (data, cmap, cbar_label, title_suffix) in enumerate([
            (abs_err, 'viridis', 'log$_{10}$(|abs error|)', 'absolute error'),
            (rel_err, 'plasma',  'log$_{10}$(|rel error|)', 'relative error'),
        ]):
            ax = axes[row, col]
            log_data = np.log10(np.maximum(data, 1e-16))

            # Build uniform levels spanning the data range
            vmin = np.nanpercentile(log_data, 2)    # avoid extreme outliers
            vmax = np.nanpercentile(log_data, 98)
            levels = np.linspace(vmin, vmax, 21)

            cf = ax.contourf(MU, BETA, log_data, levels=levels, cmap=cmap, extend='both')
            cb = fig.colorbar(cf, ax=ax)
            cb.set_label(cbar_label, fontsize=9)
            cb.locator = ticker.MaxNLocator(nbins=6)
            cb.update_ticks()

            ax.set_xlabel(r'$\mu$')
            ax.set_ylabel(r'$\beta$')
            ax.set_title(f'{label}  —  {title_suffix}')

    fig.tight_layout()
    plt.show()


def make_plots(beta_sub, mu, abs_err, rel_err):
    """Generate colour-map and line plots of the errors."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not available – skipping plots.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Absolute error map
    ax = axes[0]
    im = ax.pcolormesh(
        mu, beta_sub, np.log10(np.maximum(abs_err, 1e-16)),
        cmap='viridis', shading='auto'
    )
    fig.colorbar(im, ax=ax, label='log10(|abs error|)')
    ax.set_xlabel('mu')
    ax.set_ylabel('beta')
    ax.set_title('Absolute error  |asymptotic - reference|')

    # Relative error map
    ax = axes[1]
    im = ax.pcolormesh(
        mu, beta_sub, np.log10(np.maximum(rel_err, 1e-16)),
        cmap='plasma', shading='auto'
    )
    fig.colorbar(im, ax=ax, label='log10(|rel error|)')
    ax.set_xlabel('mu')
    ax.set_ylabel('beta')
    ax.set_title('Relative error  |asymptotic - reference| / |reference|')

    plt.tight_layout()
    plt.show()

    # Line plots: error vs beta for a few representative mu values
    mu_indices = [0, 24, 49, 74, 99]
    fig, ax = plt.subplots(figsize=(10, 5))
    for idx in mu_indices:
        ax.semilogy(beta_sub, abs_err[:, idx], label=f'mu={mu[idx]:.4f}')
    ax.set_xlabel('beta')
    ax.set_ylabel('|absolute error|')
    ax.set_title('Absolute error vs beta for selected mu values')
    ax.legend()
    ax.grid(True, which='both', ls='--', alpha=0.5)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--beta-min', type=float, default=BETA_MIN,
                        help=f'Minimum beta for comparison (default: {BETA_MIN})')
    parser.add_argument('--plot', action='store_true',
                        help='Generate error plots using matplotlib')
    parser.add_argument('--skip-table', action='store_true',
                        help='Skip the Gt.h5 table comparison (run only derivative FD tests)')
    parser.add_argument('--test-deriv', action='store_true',
                        help='Run finite-difference consistency tests for derivative functions')
    parser.add_argument('--compare-deriv-tables', action='store_true',
                        help='Compare dGdtt/dGdttt/dGdtx/dGdttx against their tabulated HDF5 files')
    parser.add_argument('--test-simplify', action='store_true',
                        help='Test high-beta simplifications: f0 negligibility and nmax convergence (beta 55-60)')
    args = parser.parse_args()

    if not args.skip_table:
        beta_sub, mu, fcn_sub, computed, abs_err, rel_err = compare_dGdt_asymptotic(
            beta_min=args.beta_min
        )
        print_statistics(beta_sub, mu, fcn_sub, computed, abs_err, rel_err)

        if args.plot:
            make_plots(beta_sub, mu, abs_err, rel_err)

    if args.compare_deriv_tables:
        compare_all_derivative_tables(beta_min=args.beta_min, plot=args.plot)

    all_ok = True
    if args.test_deriv:
        all_ok = test_derivatives_fd() and all_ok
    if args.test_simplify:
        all_ok = test_high_beta_simplifications() and all_ok
    if args.test_deriv or args.test_simplify:
        sys.exit(0 if all_ok else 1)
