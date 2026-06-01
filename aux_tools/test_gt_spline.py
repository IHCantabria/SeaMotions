#
# Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
#
# This file is part of SeaMotions Software.
#
# SeaMotions is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SeaMotions is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#

"""
Test 2D spline interpolation of the Gt tabulated function.

Dataset  (aux_data/0_integrals_database/1_time_domain/Gt.h5):
  - beta : 1-D array (6000,)   range [0, 59.99]    (uniform spacing)
  - mu   : 1-D array (100,)    range [0.0001, 0.9999]
  - fcn  : 2-D array (6000, 100)  = Gt(beta, mu)

Coordinate transform:
  mu spans ~4 decades so the spline is built in log10(mu) space to keep
  uniform resolution along the y axis.

Methods compared
  A. RectBivariateSpline  kx=3, ky=3  (cubic, interpolating  s=0)
  B. RectBivariateSpline  kx=5, ky=5  (quintic, interpolating s=0)
  C. RegularGridInterpolator  'linear'  (bilinear  – baseline)
  D. RegularGridInterpolator  'cubic'   (bicubic on regular grid)

Benchmarks reported
  1. Construction time  (tracemalloc peak + average over N_BUILD_REPEAT builds)
  2. Batch evaluation time vs batch size (N = 1 ... 10 000)
  3. Single-point evaluation time (timeit, N_SINGLE_CALLS repetitions)
  4. On-grid accuracy   - evaluated on the original 6000 x 100 nodes
  5. Off-grid accuracy  - 50 000 random interior points, ref = linear RGI
  6. Knot count (RectBivariateSpline only)
"""

import gc
import os
import time
import timeit
import tracemalloc

import math as _math

import h5py
import numpy as np
import scipy.interpolate as si
import scipy.special as _sps


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

ROOT    = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
H5_PATH = os.path.join(ROOT, "aux_data", "0_integrals_database", "1_time_domain", "Gt.h5")


# ---------------------------------------------------------------------------
# Data loader
# ---------------------------------------------------------------------------

def load_data():
    with h5py.File(H5_PATH, "r") as f:
        beta = f["beta"][:]
        mu   = f["mu"][:]
        fcn  = f["fcn"][:]
    return beta, mu, fcn


# ---------------------------------------------------------------------------
# Evaluator functions  (uniform signature: eval_fn(interp, pts) -> 1-D array)
# ---------------------------------------------------------------------------

def eval_rbs(interp, pts):
    """RectBivariateSpline:  pts shape (N, 2) = (beta, log10_mu)."""
    return interp.ev(pts[:, 0], pts[:, 1])


def eval_rgi(interp, pts):
    """RegularGridInterpolator:  pts shape (N, 2) = (beta, log10_mu)."""
    return interp(pts)


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def build_rbs(beta, log_mu, fcn, kx=3, ky=3):
    return si.RectBivariateSpline(beta, log_mu, fcn, kx=kx, ky=ky, s=0)


def build_rgi(beta, log_mu, fcn, method="linear"):
    return si.RegularGridInterpolator(
        (beta, log_mu), fcn, method=method, bounds_error=False)


# ---------------------------------------------------------------------------
# Benchmark helpers
# ---------------------------------------------------------------------------

BATCH_SIZES    = [1, 10, 100, 1_000, 10_000]
N_SINGLE_CALLS = 10_000
N_BUILD_REPEAT = 5


def time_build(builder, repeat=N_BUILD_REPEAT):
    """Return (mean_s, last_object)."""
    times = []
    obj   = None
    for _ in range(repeat):
        t0  = time.perf_counter()
        obj = builder()
        times.append(time.perf_counter() - t0)
    return float(np.mean(times)), obj


def peak_memory_build(builder):
    """Return (peak_bytes, object) using tracemalloc."""
    gc.collect()
    tracemalloc.start()
    obj = builder()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak, obj


def time_batch_eval(eval_fn, interp, pts_pool, repeat=5):
    """
    For each batch size, sample that many rows from pts_pool and time the
    evaluation averaged over `repeat` runs.
    Returns dict {batch_size: mean_seconds}.
    """
    rng = np.random.default_rng(42)
    results = {}
    for bs in BATCH_SIZES:
        idx  = rng.choice(len(pts_pool), size=bs, replace=(bs > len(pts_pool)))
        pts  = pts_pool[idx]
        runs = []
        for _ in range(repeat):
            t0 = time.perf_counter()
            eval_fn(interp, pts)
            runs.append(time.perf_counter() - t0)
        results[bs] = float(np.mean(runs))
    return results


# ---------------------------------------------------------------------------
# Accuracy helpers
# ---------------------------------------------------------------------------

def on_grid_accuracy(eval_fn, interp, beta, log_mu, fcn):
    """
    Evaluate at every original grid node; compare to tabulated fcn.
    Returns (max_abs, mean_abs, max_rel, mean_rel).
    """
    B, M = np.meshgrid(beta, log_mu, indexing="ij")
    pts  = np.column_stack([B.ravel(), M.ravel()])
    pred = eval_fn(interp, pts).reshape(fcn.shape)
    diff = np.abs(pred - fcn)
    ref  = np.abs(fcn)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(ref > 1e-14, diff / ref, diff)
    return diff.max(), diff.mean(), rel.max(), rel.mean()


def off_grid_accuracy(eval_fn, interp, ref_interp, beta, log_mu, n_pts=50_000, seed=7):
    """
    Sample n_pts random interior points, compare interpolant to linear-RGI
    reference. Returns (max_abs, mean_abs).
    """
    rng      = np.random.default_rng(seed)
    pts      = np.column_stack([
        rng.uniform(beta[0],   beta[-1],   n_pts),
        rng.uniform(log_mu[0], log_mu[-1], n_pts),
    ])
    ref_vals  = ref_interp(pts)
    pred_vals = eval_fn(interp, pts)
    diff      = np.abs(pred_vals - ref_vals)
    return diff.max(), diff.mean()


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def fmt_time(s):
    if s >= 1.0:    return f"{s:.3f} s"
    if s >= 1e-3:   return f"{s*1e3:.3f} ms"
    if s >= 1e-6:   return f"{s*1e6:.3f} us"
    return              f"{s*1e9:.1f} ns"


def fmt_bytes(b):
    for unit, thr in [("GB", 1e9), ("MB", 1e6), ("KB", 1e3)]:
        if b >= thr:
            return f"{b/thr:.2f} {unit}"
    return f"{b} B"


def fmt_e(v):
    return f"{v:.3e}"


def _section(title):
    print(f"\n{'='*72}\n  {title}\n{'='*72}")


def _row(label, *vals, lw=36):
    print(f"  {label:<{lw}}" + "".join(f"  {str(v):>15}" for v in vals))


# ---------------------------------------------------------------------------
# Section 7 helpers – asymptotic f0/f2 decomposition with configurable nmax
# ---------------------------------------------------------------------------

def _asym_alpha(n):
    """Expansion coefficient  α(n) = (2n+2)! / n! ."""
    return float(_math.factorial(2*n + 2)) / float(_math.factorial(n))


def _asym_dmn(m, n):
    """Double expansion coefficient d(m, n) from the asymptotic series."""
    if m == 0 and n == 0:
        return 1.0
    if n == 0:
        return 0.0
    return (float(_math.factorial(2*m + 2*n - 2))
            / (float(_math.factorial(2*n - 2)) * 4.0**m * float(_math.factorial(m))))


def _f0_term(beta, mu, nmax, order):
    """
    Algebraic f0 part of the dGdt asymptotic expansion or its beta-derivative.

    order=0: f0       = -4 · Σ α(n)·β^{-2n-3}·P_n(μ)                   (dGdt)
    order=1: df0/dβ   = +4 · Σ α(n)·(2n+3)·β^{-2n-4}·P_n(μ)           (dGdtt)
    order=2: d²f0/dβ² = -4 · Σ α(n)·(2n+3)·(2n+4)·β^{-2n-5}·P_n(μ)  (dGdttt)
    """
    val = 0.0
    for n in range(nmax):
        coef = _asym_alpha(n)
        if order >= 1:
            coef *= (2*n + 3)
        if order >= 2:
            coef *= (2*n + 4)
        val += coef * beta**(-2*n - 3 - order) * _sps.eval_legendre(n, mu)
    return (-4.0 if order % 2 == 0 else 4.0) * val


def _f2_dGdt(beta, mu, nmax):
    """
    Oscillatory part of the dGdt asymptotic expansion: Re[f2A · f2B].
    Uses configurable nmax series terms (reference uses nmax=15).
    """
    eta = np.sqrt(1.0 - mu**2)
    f2a = (-1j * 4.0
           * np.exp(-0.25 * beta**2 * mu)
           * np.exp(1j * (0.25*beta**2*eta - 0.5*np.arccos(mu) + np.pi/4.0))
           / np.sqrt(2.0 * eta))
    f2b = 0.0 + 0j
    for n in range(nmax):
        r0 = (1j / eta)**n
        r1 = 0.0 + 0j
        for m in range(n + 1):
            r1 += (_asym_dmn(m, n)
                   * (0.5 * beta * (1j*mu + eta))**(1.0 - 2.0*m - 2.0*n)
                   * (mu - 1j*eta)**m)
        f2b += r0 * r1
    return float((f2a * f2b).real)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    # ---- load ---------------------------------------------------------------
    _section("Loading data")
    t0 = time.perf_counter()
    beta, mu, fcn = load_data()
    load_time = time.perf_counter() - t0
    log_mu = np.log10(mu)

    n_beta, n_mu = beta.shape[0], mu.shape[0]
    print(f"  beta      : {beta.shape}  [{beta[0]:.4f}, {beta[-1]:.4f}]  "
          f"step={beta[1]-beta[0]:.4f}")
    print(f"  mu        : {mu.shape}    [{mu[0]:.4e}, {mu[-1]:.4f}]")
    print(f"  log10(mu) :             [{log_mu[0]:.3f}, {log_mu[-1]:.3f}]")
    print(f"  fcn       : {fcn.shape}  [{fcn.min():.4f}, {fcn.max():.4f}]")
    print(f"  Raw data size : {fmt_bytes(fcn.nbytes + beta.nbytes + mu.nbytes)}")
    print(f"  Load time     : {fmt_time(load_time)}")

    # ---- reference linear RGI (used in off-grid accuracy) -------------------
    ref_rgi = build_rgi(beta, log_mu, fcn, method="linear")

    # ---- pool of all grid points for batch timing ---------------------------
    B_g, M_g  = np.meshgrid(beta, log_mu, indexing="ij")
    pts_grid  = np.column_stack([B_g.ravel(), M_g.ravel()])

    # ---- method registry ----------------------------------------------------
    methods = [
        dict(label="RectBivSpline  cubic  (k=3,3)",
             builder=lambda: build_rbs(beta, log_mu, fcn, kx=3, ky=3),
             eval=eval_rbs, is_rbs=True),
        dict(label="RectBivSpline  quintic(k=5,5)",
             builder=lambda: build_rbs(beta, log_mu, fcn, kx=5, ky=5),
             eval=eval_rbs, is_rbs=True),
        dict(label="RegularGridInterp  linear",
             builder=lambda: build_rgi(beta, log_mu, fcn, method="linear"),
             eval=eval_rgi, is_rbs=False),
        dict(label="RegularGridInterp  cubic",
             builder=lambda: build_rgi(beta, log_mu, fcn, method="cubic"),
             eval=eval_rgi, is_rbs=False),
    ]

    # =========================================================================
    # 1.  Construction cost
    # =========================================================================
    _section(f"Construction cost  (average over {N_BUILD_REPEAT} builds)")
    _row("Method", "Build time (avg)", "Peak memory", "Knots (beta x mu)")
    print("  " + "-" * 72)
    for m in methods:
        peak_mem, interp = peak_memory_build(m["builder"])
        mean_t,  _       = time_build(m["builder"])
        m["interp"]      = interp
        m["peak_mem"]    = peak_mem

        if m["is_rbs"]:
            kx_k, ky_k = (len(k) for k in interp.get_knots())
            knot_str = f"{kx_k} x {ky_k}"
        else:
            knot_str = "N/A"

        _row(m["label"], fmt_time(mean_t), fmt_bytes(peak_mem), knot_str)

    # =========================================================================
    # 2.  Batch evaluation time
    # =========================================================================
    _section("Batch evaluation time  (mean of 5 runs per batch size)")
    _row("Method", *[f"N={bs}" for bs in BATCH_SIZES])
    print("  " + "-" * 72)
    for m in methods:
        times = time_batch_eval(m["eval"], m["interp"], pts_grid)
        _row(m["label"], *[fmt_time(times[bs]) for bs in BATCH_SIZES])

    # ---- throughput summary at largest batch --------------------------------
    _section("Throughput at N=10 000")
    _row("Method", "Total time", "Evals/sec")
    print("  " + "-" * 72)
    for m in methods:
        times = time_batch_eval(m["eval"], m["interp"], pts_grid)
        t     = times[10_000]
        _row(m["label"], fmt_time(t), f"{10_000/t:,.0f}")

    # =========================================================================
    # 3.  Single-point evaluation time
    # =========================================================================
    _section(f"Single-point evaluation time  (timeit, {N_SINGLE_CALLS:,} calls)")
    rng = np.random.default_rng(0)
    pt  = np.array([[
        float(rng.uniform(beta[0],   beta[-1])),
        float(rng.uniform(log_mu[0], log_mu[-1])),
    ]])
    _row("Method", "Mean per-call", "Calls/sec")
    print("  " + "-" * 72)
    for m in methods:
        ev    = m["eval"]
        itp   = m["interp"]
        total = timeit.timeit(lambda: ev(itp, pt), number=N_SINGLE_CALLS)
        tc    = total / N_SINGLE_CALLS
        _row(m["label"], fmt_time(tc), f"{1/tc:,.0f}")

    # =========================================================================
    # 4.  On-grid accuracy
    # =========================================================================
    _section(f"On-grid accuracy  ({n_beta} x {n_mu} = {n_beta*n_mu:,} nodes)")
    _row("Method", "max|abs|", "mean|abs|", "max|rel|", "mean|rel|")
    print("  " + "-" * 72)
    for m in methods:
        r = on_grid_accuracy(m["eval"], m["interp"], beta, log_mu, fcn)
        _row(m["label"], *[fmt_e(v) for v in r])

    # =========================================================================
    # 5.  Off-grid accuracy
    # =========================================================================
    _section("Off-grid accuracy  (50 000 random interior points, ref = linear RGI)")
    print("  Note: differences here reflect higher-order vs bilinear, not the true error.")
    _row("Method", "max|abs|", "mean|abs|")
    print("  " + "-" * 72)
    for m in methods:
        if m["label"].startswith("RegularGridInterp  linear"):
            _row(m["label"], "(reference)", "---")
            continue
        mae, mean_ae = off_grid_accuracy(
            m["eval"], m["interp"], ref_rgi, beta, log_mu)
        _row(m["label"], fmt_e(mae), fmt_e(mean_ae))

    # =========================================================================
    # 6.  Memory summary
    # =========================================================================
    _section("Memory summary")
    raw_size = fmt_bytes(fcn.nbytes + beta.nbytes + mu.nbytes)
    print("  Note: 'Peak memory' above (tracemalloc) is the most reliable measure.")
    print("  sys.getsizeof returns only the Python header size (48 B for all wrappers).")
    print("  The underlying C/Fortran arrays are shared with the raw table.")
    print()
    _row("Method", "Raw table (shared)", "Peak (tracemalloc)")
    print("  " + "-" * 72)
    # Recover peak from the already-captured peak_memory_build results
    for m in methods:
        _row(m["label"], raw_size, fmt_bytes(m["peak_mem"]))

    # =========================================================================
    # 7.  Asymptotic expansion in beta ∈ [55, 60]: f0 negligibility + nmax
    #
    # Gt(β,μ) ≈ f0 + Re[f2A·f2B]  where
    #   f0  is algebraic (decays as β^{-3})        ← current code omits this
    #   f2  is oscillatory (decays as exp(-β²μ/4)) ← current code returns this
    #
    # 7a verifies the absolute magnitude of f0 (and its β-derivatives).
    # 7b verifies that the f2 series converges with far fewer than 15 terms.
    # =========================================================================
    _section("Asymptotic expansion study for beta in [55, 60]")

    B55_60   = np.linspace(55.0, 60.0, 10)
    MU_GRID  = np.logspace(-4.0, np.log10(0.9999), 50)  # log-spaced [0.0001, 0.9999]
    NMAX_REF = 15

    # ---- 7a: absolute magnitude of f0 and its derivatives ------------------
    print("  7a.  Absolute magnitude of the f0 (algebraic) term  [nmax=15]")
    print("       Grid: 10 beta x 50 mu (log-spaced).\n")
    _row("Term", "max |value|", "beta at max", "mu at max")
    print("  " + "-"*62)
    for order, label in [(0, "f0         (dGdt  algebraic)"),
                         (1, "df0/dbeta  (dGdtt algebraic)"),
                         (2, "d2f0/dbeta2 (dGdttt algebraic)")]:
        grid = np.array([[abs(_f0_term(b, m, NMAX_REF, order))
                          for m in MU_GRID]
                         for b in B55_60])
        idx = np.unravel_index(grid.argmax(), grid.shape)
        _row(label, fmt_e(grid.max()),
             f"{B55_60[idx[0]]:.2f}",
             f"{MU_GRID[idx[1]]:.4e}")

    # ---- 7b: nmax convergence of f2 ----------------------------------------
    _section("  7b.  nmax convergence of f2 (oscillatory, dGdt)  –  ref nmax=15")
    print("  Active grid points = those where |f2(nmax=15)| > 1e-10.\n")

    t0_7b  = time.perf_counter()
    ref_f2 = np.array([[_f2_dGdt(b, m, NMAX_REF) for m in MU_GRID]
                        for b in B55_60])
    active = np.abs(ref_f2) > 1e-10
    print(f"  Active points : {active.sum()} / {active.size}  "
          f"(build ref took {fmt_time(time.perf_counter() - t0_7b)})\n")

    TOL_NMAX = 1e-5
    _row("nmax", "max |Δf2/f2_15|", "mean |Δf2/f2_15|", f"PASS (<{TOL_NMAX:.0e})")
    print("  " + "-"*62)
    all_nmax_pass = True
    for nm in [1, 2, 3, 4, 5, 8, 10]:
        f2_nm = np.array([[_f2_dGdt(b, m, nm) for m in MU_GRID]
                           for b in B55_60])
        with np.errstate(divide="ignore", invalid="ignore"):
            rel = np.where(active,
                           np.abs(f2_nm - ref_f2) / (np.abs(ref_f2) + 1e-300),
                           0.0)
        max_r  = float(rel.max())
        mean_r = float(rel[active].mean()) if active.any() else 0.0
        ok     = max_r < TOL_NMAX
        if nm <= 5:
            all_nmax_pass = all_nmax_pass and ok
        _row(f"nmax={nm:2d}", fmt_e(max_r), fmt_e(mean_r), "PASS" if ok else "FAIL")

    print()
    print(f"  >> All nmax <= 5 achieve < {TOL_NMAX:.0e} rel error ?  "
          f"{'PASS' if all_nmax_pass else 'FAIL'}")

    print("\nDone.\n")


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
