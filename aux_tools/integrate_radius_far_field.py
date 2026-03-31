"""Numerical integration for triple Hankel products over [1, inf)."""

from __future__ import annotations

from typing import Iterable, Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import hankel1, hankel2


def _romberg(
    func,
    a: float,
    b: float,
    *,
    max_iter: int = 10,
    rtol: float = 1e-8,
    atol: float = 1e-12,
) -> Tuple[float, float]:
    """Romberg integration for a real-valued function on [a, b]."""
    if a == b:
        return 0.0, 0.0

    r = np.zeros((max_iter, max_iter), dtype=float)
    h = b - a
    r[0, 0] = 0.5 * h * (func(a) + func(b))

    for i in range(1, max_iter):
        h *= 0.5
        step = 2 ** (i - 1)
        sum_f = 0.0
        for k in range(1, step + 1):
            sum_f += func(a + (2 * k - 1) * h)
        r[i, 0] = 0.5 * r[i - 1, 0] + h * sum_f

        for j in range(1, i + 1):
            factor = 4 ** j
            r[i, j] = r[i, j - 1] + (r[i, j - 1] - r[i - 1, j - 1]) / (factor - 1)

        if abs(r[i, i] - r[i - 1, i - 1]) <= atol + rtol * abs(r[i, i]):
            return r[i, i], abs(r[i, i] - r[i - 1, i - 1])

    return r[max_iter - 1, max_iter - 1], abs(r[max_iter - 1, max_iter - 1] - r[max_iter - 2, max_iter - 2])


def _quad_complex(
    func,
    a: float,
    b: float,
    *,
    method: str = "quad",
    **kwargs,
) -> Tuple[complex, float]:
    """Integrate a complex-valued function by splitting real/imag parts."""
    if method == "romberg":
        real_val, real_err = _romberg(lambda x: np.real(func(x)), a, b, **kwargs)
        imag_val, imag_err = _romberg(lambda x: np.imag(func(x)), a, b, **kwargs)
        return real_val + 1j * imag_val, real_err + imag_err

    real_val, real_err = quad(lambda x: np.real(func(x)), a, b, **kwargs)
    imag_val, imag_err = quad(lambda x: np.imag(func(x)), a, b, **kwargs)
    return real_val + 1j * imag_val, real_err + imag_err


def hankel_asymp(n: int, z: complex, kind: int) -> complex:
    """Asymptotic approximation for Hankel function H_n(z) for large |z|."""
    if z == 0:
        return 0.0
    if kind not in (1, 2):
        raise ValueError("kind must be 1 or 2.")
    phase = z - n * np.pi / 2.0 - np.pi / 4.0
    sign = 1.0 if kind == 1 else -1.0
    return np.sqrt(2.0 / (np.pi * z)) * np.exp(1j * sign * phase)


def _compute_rotated_contour_parameters(
    n1: int,
    n2: int,
    a: float,
    b: float,
    c: float,
    hkind0: int,
    hkind1: int,
    *,
    tail_cycles: float = 30.0,
) -> Tuple[float, float, float, float]:
    """Return sigma, sigma_max, segment length, and x_max for the rotated contour."""
    if hkind0 not in (1, 2) or hkind1 not in (1, 2):
        raise ValueError("hkind0 and hkind1 must be either 1 or 2.")

    sign0 = 1.0 if hkind0 == 1 else -1.0
    sign1 = 1.0 if hkind1 == 1 else -1.0

    sigma = sign0 * a + sign1 * b - c
    sigma_max = max(abs(a) + abs(b), 1e-15)
    segment_len = 2.0 * np.pi / sigma_max / 10.0

    abs_sigma = abs(sigma)
    tail_scale = (
        tail_cycles * 2.0 * np.pi / abs_sigma
        if abs_sigma > 1e-14
        else tail_cycles * 2.0 * np.pi / sigma_max
    )

    a_scale = max(abs(a), 1e-15)
    b_scale = max(abs(b), 1e-15)
    x_max = max(n1 / a_scale, n2 / b_scale, tail_scale)
    x_max = 1.0 if x_max < 1.1 else x_max

    return sigma, sigma_max, segment_len, x_max


def triple_hankel_integrand(
    r: float,
    n1: int,
    n2: int,
    n3: int,
    a: float,
    b: float,
    c: float,
    hkind0: int = 1,
    hkind1: int = 1,
) -> complex:
    """Return r * H_n1(a r) * H_n2(b r) * H_n3(c r) with optional damping."""
    if r <= 0.0:
        return 0.0
    if hkind0 not in (1, 2) or hkind1 not in (1, 2):
        raise ValueError("hkind0 and hkind1 must be either 1 or 2.")

    h0 = hankel1 if hkind0 == 1 else hankel2
    h1 = hankel1 if hkind1 == 1 else hankel2
    val = r * h0(n1, a * r) * h1(n2, b * r) * hankel2(n3, c * r)
    
    return val


def triple_hankel_integrand_asymp(
                                y: float,
                                x_max: float,
                                sigma: float
                            ) -> complex:
    """Return r * H_n1(a r) * H_n2(b r) * H_n3(c r) with optional damping."""
    if sigma >= 0.0:
        return 1j * np.exp( sigma * ( -y + 1j * x_max ) ) / np.sqrt( ( x_max + 1j * y ) )
    
    return -1j * np.conj( np.exp( sigma * ( y - 1j * x_max ) ) * np.conj( 1.0 / np.sqrt( ( x_max + 1j * y ) ) ) )


def triple_hankel_integrand_asymp_plot( 
                                            r: float,
                                            n1: int,
                                            n2: int,
                                            n3: int,
                                            a: float,
                                            b: float,
                                            c: float,
                                            hkind0: int = 1,
                                            hkind1: int = 1,
                                        ) -> complex:
    """Return r * H_n1(a r) * H_n2(b r) * H_n3(c r) with optional damping."""
    if r <= 0.0:
        return 0.0
    if hkind0 not in (1, 2) or hkind1 not in (1, 2):
        raise ValueError("hkind0 and hkind1 must be either 1 or 2.")

    val = r * hankel_asymp(n1, a * r, hkind0) * hankel_asymp(n2, b * r, hkind1) * hankel_asymp(n3, c * r, 2)
    
    return val
    

def integrate_triple_hankel(
    n1: int,
    n2: int,
    n3: int,
    a: float,
    b: float,
    c: float,
    hkind0: int = 1,
    hkind1: int = 1,
    *,
    r_min: float = 1.0,
    rtol: float = 1e-6,
    atol: float = 1e-10,
    max_segments: int = 10000,
    min_segments: int = 4,
    segment_cycles: float = 2.0,
    quad_limit: int = 200,
) -> complex:
    """Integrate r * H_n1(a r) * H_n2(b r) * H_n3(c r) over [1, inf).

    The integral is approximated by summing fixed-length segments. The segment
    length is based on the shortest oscillation period from a, b, c.

    Parameters
    ----------
    n1, n2, n3
        Hankel orders (0..32 recommended).
    a, b, c
        Hankel arguments scaling factors.
    r_min
        Lower cutoff to avoid the r=0 singularity; use >= 1 for [1, inf).
    rtol, atol
        Convergence tolerances for segment contributions.
    max_segments, min_segments
        Segment iteration controls.
    segment_cycles
        Number of oscillation cycles per segment.
    damping
        Optional exponential damping factor to improve convergence.
    quad_limit
        Subinterval limit passed to quad.
    """
    k_max = max(abs(a), abs(b), abs(c))
    if k_max == 0.0:
        return 0.0

    period = 2.0 * np.pi / k_max
    segment_len = max(period * segment_cycles, r_min)

    total = 0.0 + 0.0j
    r_start = max(r_min, 0.0)

    for seg_idx in range(max_segments):
        r_end = r_start + segment_len
        seg_val, _ = _quad_complex(
            lambda r: triple_hankel_integrand(
                r, n1, n2, n3, a, b, c, hkind0, hkind1
            ),
            r_start,
            r_end,
            limit=quad_limit,
        )
        total += seg_val

        if seg_idx + 1 >= min_segments:
            if abs(seg_val) <= atol + rtol * abs(total):
                break

        r_start = r_end

    return total


def integrate_triple_hankel_mod(
                                            l: int,
                                            m: int,
                                            n: int,
                                            alpha: float,
                                            beta: float,
                                            gamma: float,
                                            hkind0: int,
                                            hkind1: int,
                                            *,
                                            r_min: float = 1.0,
                                            rtol: float = 1e-6,
                                            atol: float = 1e-10,
                                            max_segments: int = 10000,
                                            min_segments: int = 4,
                                            segment_cycles: float = 2.0,
                                            quad_limit: int = 200,
                                            verbose: bool = False,
                                        ) -> complex:
    """Integrate r * H_n1(a r) * H_n2(b r) * H_n3(c r) over [1, inf).

    The integral is approximated by summing fixed-length segments. The segment
    length is based on the shortest oscillation period from a, b, c.

    Parameters
    ----------
    n1, n2, n3
        Hankel orders (0..32 recommended).
    a, b, c
        Hankel arguments scaling factors.
    r_min
        Lower cutoff to avoid the r=0 singularity; use >= 1 for [1, inf).
    rtol, atol
        Convergence tolerances for segment contributions.
    max_segments, min_segments
        Segment iteration controls.
    segment_cycles
        Number of oscillation cycles per segment.
    damping
        Optional exponential damping factor to improve convergence.
    quad_limit
        Subinterval limit passed to quad.
    """

    sigma, sigma_max, segment_len, x_max = _compute_rotated_contour_parameters(
        l,
        m,
        alpha,
        beta,
        gamma,
        hkind0,
        hkind1,
        tail_cycles=30.0,
    )

    sign0 = 1.0 if hkind0 == 1 else -1.0
    sign1 = 1.0 if hkind1 == 1 else -1.0

    # Integrate first part from 1 to x_max
    finint = 0.0 + 0.0j
    head_start = max(r_min, 1.0)
    if x_max > head_start:
        dx              = x_max - head_start
        segments_np     = int(np.ceil(dx / segment_len))
        segment_len_fin = dx / segments_np
        for i in range(segments_np):
            start       = head_start + i * segment_len_fin
            end         = head_start + (i + 1) * segment_len_fin
            finint_, _  = _quad_complex(
                                            lambda x: triple_hankel_integrand(
                                                                                    x, l, m, n, alpha, beta, gamma, hkind0, hkind1
                                                                                ),
                                            start,
                                            end,
                                            limit=quad_limit,
                                        )
            finint += finint_

            if verbose:
                print( f"Head segment {i}: {finint_} - Total: {finint} - start: {start}, end: {end}" )

    # Integrate semi-infinite part from x_max to infinity
    finsint = 0.0 + 0.0j
    r_start = 0.0
    for seg_idx in range(max_segments):
        r_end = r_start + segment_len
        seg_val, _ = _quad_complex(
            lambda r: triple_hankel_integrand_asymp(
                                                       r, x_max, sigma
                                                    ),
            r_start,
            r_end,
            limit=quad_limit,
        )
        finsint += seg_val

        if verbose:
            print( f"Tail segment {seg_idx}: {seg_val}" )

        if seg_idx + 1 >= min_segments:
            if abs(seg_val) <= atol + rtol * abs(finsint):
                break

        r_start = r_end

    # Calculate scaling for semi-infinite part
    order_phase = ( -1 * sign0 * l ) + ( -1 * sign1 * m ) + n
    one_phase   = ( -1 * sign0 ) + ( -1 * sign1 ) + 1.0
    semi_scale  = (
                    (2.0 / np.pi ) 
                    * 
                    np.sqrt( 2.0 / np.pi )
                    *
                    np.sqrt( 1.0 / ( alpha * beta * gamma ) )
                    *
                    np.exp( 
                                1j * order_phase * np.pi / 2.0
                                +
                                1j * one_phase * np.pi / 4.0
                            )
                )

    if verbose:
        print( finint, semi_scale * finsint, finint + semi_scale * finsint )

    return finint + semi_scale * finsint


def plot_triple_hankel_integrand(
    n1: int,
    n2: int,
    n3: int,
    a: float,
    b: float,
    c: float,
    hkind0: int = 1,
    hkind1: int = 1,
    *,
    r_max: float,
    num_points: int = 2000,
    parts: Tuple[str, ...] = ("real", "imag", "abs"),
    show: bool = True,
):
    """Plot the integrand r * H_n1(a r) * H_n2(b r) * H_n3(c r).

    Parameters
    ----------
    r_max
        Upper plot limit for r (lower bound is fixed at 1).
    parts
        Any of: "real", "imag", "abs", "phase".
    """

    r = np.linspace(1.0, r_max, num_points)
    vals = np.array(
        [triple_hankel_integrand(x, n1, n2, n3, a, b, c, hkind0, hkind1) for x in r]
    )
    valsv = np.array(
        [triple_hankel_integrand_asymp_plot(x, n1, n2, n3, a, b, c, hkind0, hkind1) for x in r]
    )

    plt.figure(figsize=(10, 6))
    if "real" in parts:
        plt.plot(r, np.real(vals), label="real")
        plt.plot(r, np.real(valsv), label="real (asymp)")
    if "imag" in parts:
        plt.plot(r, np.imag(vals), label="imag")
        plt.plot(r, np.imag(valsv), label="imag (asymp)")
    if "abs" in parts:
        plt.plot(r, np.abs(vals), label="abs")
        plt.plot(r, np.abs(valsv), label="abs (asymp)")
    if "phase" in parts:
        plt.plot(r, np.angle(vals), label="phase")
        plt.plot(r, np.angle(valsv), label="phase (asymp)")

    plt.xlabel("r")
    plt.ylabel("integrand")
    plt.title("Triple Hankel integrand")
    plt.grid(True, alpha=0.3)
    plt.legend()
    if show:
        plt.show()
    return plt.gca()


def integrate_triple_hankel_orders(
    orders: Iterable[int],
    a: float,
    b: float,
    c: float,
    **kwargs,
) -> Dict[Tuple[int, int, int], complex]:
    """Compute the integral for all combinations of orders in the iterable."""
    orders_list = list(orders)
    results: Dict[Tuple[int, int, int], complex] = {}
    for n1 in orders_list:
        for n2 in orders_list:
            for n3 in orders_list:
                results[(n1, n2, n3)] = integrate_triple_hankel(
                    n1, n2, n3, a, b, c, **kwargs
                )
    return results


def compare_hankel_with_asymptotic(
    n: int,
    k: float,
    *,
    kind: int = 1,
    r_min: float = 1.0,
    r_max: float = 50.0,
    num_points: int = 2000,
    show: bool = True,
) -> Dict[str, np.ndarray]:
    """Compare H_n^{(kind)}(k r) with its asymptotic approximation over r in [r_min, r_max].

    Returns arrays for exact/asymptotic values and absolute/relative errors.
    """
    if kind not in (1, 2):
        raise ValueError("kind must be 1 or 2.")
    if r_min <= 0.0 or r_max <= r_min:
        raise ValueError("Require 0 < r_min < r_max.")
    if k == 0.0:
        raise ValueError("k must be non-zero.")

    r = np.linspace(r_min, r_max, num_points)
    z = k * r

    h_exact = hankel1(n, z) if kind == 1 else hankel2(n, z)

    # Correct asymptotic phase sign:
    # H_n^(1) ~ exp(+i*phase), H_n^(2) ~ exp(-i*phase)
    sign = 1.0 if kind == 1 else -1.0
    phase = z - n * np.pi / 2.0 - np.pi / 4.0
    h_asymp = np.sqrt(2.0 / (np.pi * z)) * np.exp(1j * sign * phase)

    abs_err = np.abs(h_exact - h_asymp)
    rel_err = abs_err / np.maximum(np.abs(h_exact), 1e-15)

    if show:
        fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax[0].plot(r, np.real(h_exact), label="Re exact")
        ax[0].plot(r, np.real(h_asymp), "--", label="Re asymp")
        ax[0].plot(r, np.imag(h_exact), label="Im exact")
        ax[0].plot(r, np.imag(h_asymp), "--", label="Im asymp")
        ax[0].set_ylabel("value")
        ax[0].set_title(f"Hankel comparison: n={n}, kind={kind}, k={k}")
        ax[0].grid(True, alpha=0.3)
        ax[0].legend()

        ax[1].semilogy(r, abs_err, label="abs error")
        ax[1].semilogy(r, rel_err, label="rel error")
        ax[1].set_xlabel("r")
        ax[1].set_ylabel("error")
        ax[1].grid(True, alpha=0.3)
        ax[1].legend()

        plt.tight_layout()
        plt.show()

    return {
        "r": r,
        "z": z,
        "exact": h_exact,
        "asymp": h_asymp,
        "abs_err": abs_err,
        "rel_err": rel_err,
    }


def integrate_triple_hankel_asymp_plot(
    n1: int,
    n2: int,
    n3: int,
    a: float,
    b: float,
    c: float,
    hkind0: int = 1,
    hkind1: int = 1,
    *,
    r_min: float = 1.0,
    rtol: float = 1e-6,
    atol: float = 1e-10,
    max_segments: int = 10000,
    min_segments: int = 4,
    segment_cycles: float = 2.0,
    quad_limit: int = 200,
) -> complex:
    """Integrate asymptotic triple-Hankel integrand on [r_min, inf)."""
    k_max = max(abs(a), abs(b), abs(c))
    if k_max == 0.0:
        return 0.0

    period = 2.0 * np.pi / k_max
    segment_len = max(period * segment_cycles, r_min)

    total = 0.0 + 0.0j
    r_start = max(r_min, 0.0)

    for seg_idx in range(max_segments):
        r_end = r_start + segment_len
        seg_val, _ = _quad_complex(
            lambda r: triple_hankel_integrand_asymp_plot(
                r, n1, n2, n3, a, b, c, hkind0, hkind1
            ),
            r_start,
            r_end,
            limit=quad_limit,
        )
        total += seg_val

        if seg_idx + 1 >= min_segments:
            if abs(seg_val) <= atol + rtol * abs(total):
                break

        r_start = r_end

    return total


def compare_triple_hankel_mod_with_asymp_plot(
    n1: int,
    n2: int,
    n3: int,
    a: float,
    b: float,
    c: float,
    hkind0: int = 1,
    hkind1: int = 1,
    *,
    r_min: float = 1.0,
    rtol: float = 1e-6,
    atol: float = 1e-10,
    max_segments: int = 20000,
    min_segments: int = 4,
    segment_cycles: float = 2.0,
    quad_limit: int = 200,
    verbose: bool = True,
) -> Dict[str, complex | float]:
    """Compare integrate_triple_hankel_mod with asymptotic-plot integrand integration."""
    val_mod = integrate_triple_hankel_mod(
        n1,
        n2,
        n3,
        a,
        b,
        c,
        hkind0,
        hkind1,
        r_min=r_min,
        rtol=rtol,
        atol=atol,
        max_segments=max_segments,
        min_segments=min_segments,
        segment_cycles=segment_cycles,
        quad_limit=quad_limit,
    )

    val_asymp_plot = integrate_triple_hankel(
        n1,
        n2,
        n3,
        a,
        b,
        c,
        hkind0,
        hkind1,
        r_min=r_min,
        rtol=rtol,
        atol=atol,
        max_segments=max_segments,
        min_segments=min_segments,
        segment_cycles=segment_cycles,
        quad_limit=quad_limit,
    )

    diff = val_mod - val_asymp_plot
    abs_diff = abs(diff)
    rel_diff = abs_diff / max(abs(val_mod), 1e-15)

    if verbose:
        print(f"integrate_triple_hankel_mod       = {val_mod}")
        print(f"integrate_triple_hankel_asymp_plot = {val_asymp_plot}")
        print(f"absolute difference                = {abs_diff:.6e}")
        print(f"relative difference                = {rel_diff:.6e}")

    return {
        "mod": val_mod,
        "asymp_plot": val_asymp_plot,
        "difference": diff,
        "abs_difference": abs_diff,
        "rel_difference": rel_diff,
    }


if __name__ == "__main__":
    import time
    # Example usage
    nu0         = 1
    nu1         = 1
    nu2         = 1
    alpha       = 2.0 * np.pi / 100.0 * 100
    beta        = 2.0 * np.pi / 100.0 * 100
    gamma       = 2.0 * np.pi / 100.0 * 100
    
    # t0          = time.perf_counter( )
    # val_ref     = integrate_triple_hankel(nu0, nu1, nu2, alpha, beta, gamma, 1, 1)
    # t1          = time.perf_counter( )
    # val_mod     = integrate_triple_hankel_mod(nu0, nu1, nu2, alpha, beta, gamma, 1, 1)
    # t2          = time.perf_counter( )

    # print(f"Reference integral: {val_ref:.6e} computed in {t1 - t0:.2f} seconds")
    # print(f"Modified integral:  {val_mod:.6e} computed in {t2 - t1:.2f} seconds")

    val_mod     = compare_triple_hankel_mod_with_asymp_plot(nu0, nu1, nu2, alpha, beta, gamma, 1, 1)
    plot_triple_hankel_integrand(nu0, nu1, nu2, alpha, beta, gamma, 1, 1, r_max=20.0)

    # compare_hankel_with_asymptotic( 1, 1 )


