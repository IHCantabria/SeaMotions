"""Generate reference databases for triple Hankel integrals."""

from __future__ import annotations

import os
import time
from itertools import product
from typing import Any, Iterable

import h5py
import numpy as np

from integrate_radius_far_field import integrate_triple_hankel, integrate_triple_hankel_mod


def get_root_fopath() -> str:
    """Return the repository root folder."""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def get_triple_hankel_reference_database_fopath() -> str:
    """Return the default folder used to store triple-Hankel reference datasets."""
    return os.path.join(
        get_root_fopath(),
        "aux_data",
        "2_check_database_files",
        "0_frequency_domain",
    )


def _as_array(values: Iterable[float] | Iterable[int], dtype: Any) -> np.ndarray:
    array = np.asarray(list(values), dtype=dtype)
    if array.ndim != 1 or array.size == 0:
        raise ValueError("Expected a non-empty one-dimensional iterable.")
    return array


def _format_seconds(seconds: float) -> str:
    total = max(0, int(seconds))
    h = total // 3600
    m = (total % 3600) // 60
    s = total % 60
    if h > 0:
        return f"{h:d}:{m:02d}:{s:02d}"
    return f"{m:02d}:{s:02d}"


def _print_progress(completed: int, total: int, t0: float, bar_width: int = 30) -> None:
    frac = completed / max(total, 1)
    filled = int(round(bar_width * frac))
    bar = "#" * filled + "-" * (bar_width - filled)
    elapsed = time.time() - t0
    rate = completed / elapsed if elapsed > 0 else 0.0
    remaining = (total - completed) / rate if rate > 0 else 0.0
    msg = (
        f"\r[{bar}] {100.0 * frac:6.2f}% "
        f"({completed}/{total}) "
        f"elapsed {_format_seconds(elapsed)} "
        f"eta {_format_seconds(remaining)}"
    )
    print(msg, end="", flush=True)


def create_triple_hankel_reference_database(
    file_path: str,
    orders_n1: Iterable[int],
    orders_n2: Iterable[int],
    orders_n3: Iterable[int],
    alpha_values: Iterable[float],
    beta_values: Iterable[float],
    gamma_values: Iterable[float],
    *,
    hkind0: int = 2,
    hkind1: int = 2,
    include_full_integral: bool = True,
    overwrite: bool = False,
    r_min: float = 1.0,
    rtol: float = 1e-6,
    atol: float = 1e-10,
    max_segments: int = 10000,
    min_segments: int = 4,
    segment_cycles: float = 2.0,
    quad_limit: int = 200,
    show_progress: bool = True,
) -> str:
    """Create an HDF5 database of rotated-contour triple Hankel integrals.

    The database stores one regular grid of values computed with
    `integrate_triple_hankel_mod` and, optionally, a second grid from
    `integrate_triple_hankel` for cross-checking.
    """
    orders1 = _as_array(orders_n1, np.int32)
    orders2 = _as_array(orders_n2, np.int32)
    orders3 = _as_array(orders_n3, np.int32)
    alpha = _as_array(alpha_values, np.float64)
    beta = _as_array(beta_values, np.float64)
    gamma = _as_array(gamma_values, np.float64)

    if not overwrite and os.path.exists(file_path):
        raise FileExistsError(f"Database already exists: {file_path}")

    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    grid_shape = (
        orders1.size,
        orders2.size,
        orders3.size,
        alpha.size,
        beta.size,
        gamma.size,
    )
    values_mod = np.zeros(grid_shape, dtype=np.complex128)
    values_full = np.zeros(grid_shape, dtype=np.complex128) if include_full_integral else None

    total_points = int(np.prod(grid_shape))
    t0 = time.time()
    last_update = t0
    if show_progress:
        print(f"Generating triple Hankel database: {total_points} points")

    for step, idx in enumerate(product(
        range(orders1.size),
        range(orders2.size),
        range(orders3.size),
        range(alpha.size),
        range(beta.size),
        range(gamma.size),
    ), start=1):
        i1, i2, i3, ia, ib, ic = idx
        n1 = int(orders1[i1])
        n2 = int(orders2[i2])
        n3 = int(orders3[i3])
        a = float(alpha[ia])
        b = float(beta[ib])
        c = float(gamma[ic])

        values_mod[idx] = integrate_triple_hankel_mod(
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
            verbose=False,
        )

        if include_full_integral:
            assert values_full is not None
            values_full[idx] = integrate_triple_hankel(
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

        if show_progress:
            now = time.time()
            if step == total_points or (now - last_update) >= 0.25:
                _print_progress(step, total_points, t0)
                last_update = now

    if show_progress:
        print()

    with h5py.File(file_path, "w") as fid:
        fid.attrs["generator"] = "generate_triple_hankel_reference_database"
        fid.attrs["hkind0"] = hkind0
        fid.attrs["hkind1"] = hkind1
        fid.attrs["r_min"] = r_min
        fid.attrs["rtol"] = rtol
        fid.attrs["atol"] = atol
        fid.attrs["max_segments"] = max_segments
        fid.attrs["min_segments"] = min_segments
        fid.attrs["segment_cycles"] = segment_cycles
        fid.attrs["quad_limit"] = quad_limit

        fid.create_dataset("orders_n1", data=orders1)
        fid.create_dataset("orders_n2", data=orders2)
        fid.create_dataset("orders_n3", data=orders3)
        fid.create_dataset("alpha", data=alpha)
        fid.create_dataset("beta", data=beta)
        fid.create_dataset("gamma", data=gamma)
        fid.create_dataset("values_mod_real", data=np.real(values_mod))
        fid.create_dataset("values_mod_imag", data=np.imag(values_mod))

        if include_full_integral and values_full is not None:
            abs_diff = np.abs(values_mod - values_full)
            rel_diff = abs_diff / np.maximum(np.abs(values_full), 1e-15)
            fid.create_dataset("values_full_real", data=np.real(values_full))
            fid.create_dataset("values_full_imag", data=np.imag(values_full))
            fid.create_dataset("abs_difference", data=abs_diff)
            fid.create_dataset("rel_difference", data=rel_diff)

    return file_path


if __name__ == "__main__":
    output_path = os.path.join(
        get_triple_hankel_reference_database_fopath(),
        "triple_hankel_reference.h5",
    )

    k_min = 2.0 * np.pi / 1500.0
    k_max = 2.0 * np.pi
    k_values = np.linspace(k_min, k_max, 10)

    create_triple_hankel_reference_database(
        output_path,
        orders_n1=[0, 1, 5, 10, 15],
        orders_n2=[0, 1, 5, 10, 15],
        orders_n3=[0, 1, 5, 10, 15],
        alpha_values=k_values,
        beta_values=k_values,
        gamma_values=k_values,
        overwrite=True,
        include_full_integral=False
    )