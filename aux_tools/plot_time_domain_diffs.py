"""
Visualise the approximation errors of the time-domain Green's function evaluators.

Usage
-----
1. Build the test executable (CMake, ninja).

2. Run the executable with --diff-output to generate CSV files:

   cd <build_dir>
   tests\\green\\test_time_domain_evaluator_cmd.exe ^
       <dat_dir>\\dGdt.dat <dat_dir>\\dGdtx.dat <dat_dir>\\dGdtxx.dat ^
       <dat_dir>\\dGdtt.dat <dat_dir>\\dGdttx.dat <dat_dir>\\dGdttxx.dat ^
       <dat_dir>\\dGdt_G0.dat  <dat_dir>\\dGdtx_G0.dat  <dat_dir>\\dGdtxx_G0.dat ^
       <dat_dir>\\dGdtt_G0.dat <dat_dir>\\dGdttx_G0.dat <dat_dir>\\dGdttxx_G0.dat ^
       <dat_dir>\\dGdt_residual.dat  <dat_dir>\\dGdtx_residual.dat ^
       <dat_dir>\\dGdtxx_residual.dat <dat_dir>\\dGdtt_residual.dat ^
       <dat_dir>\\dGdttx_residual.dat <dat_dir>\\dGdttxx_residual.dat ^
       --diff-output <output_csv_dir>

   Replace <dat_dir> with the absolute path to tests/tests_data/time_domain/
   and <output_csv_dir> with any directory (e.g. aux_data/check_integrals).

3. Run this script:

   python aux_tools/plot_time_domain_diffs.py --diff-dir <output_csv_dir>

   Optional flags:
     --metric  abs_err | rel_err        (default: abs_err)
     --part    total | G0 | residual | all  (default: all)
     --log     log10 colour scale for the error column
     --save    save figures as PNG instead of showing interactively

Layout
------
One figure per integral function (dGdt, dGdtx, …).
Rows  = decomposition parts selected by --part  (total / G0 / residual).
Cols  = [reference value,  computed value,  error metric].

The reference and computed columns share the same colour scale and use a
diverging colormap (RdBu_r) when the data contains both positive and negative
values, or viridis when all values share the same sign.  The error column uses
a "hot_r" sequential map; failing test points are marked with a cyan cross.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FUNCTIONS   = ["dGdt", "dGdtx", "dGdtxx", "dGdtt", "dGdttx", "dGdttxx"]
PARTS       = ["total", "G0", "residual"]
PART_SUFFIX = {"total": "", "G0": "_G0", "residual": "_residual"}
PART_LABEL  = {"total": "Total", "G0": "G₀", "residual": "Residual"}

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def csv_name(fcn: str, part: str) -> str:
    return f"{fcn}{PART_SUFFIX[part]}_diff.csv"


def load_csv(path: str) -> dict:
    """Return a dict of column arrays loaded from a diff CSV."""
    data = np.genfromtxt(path, delimiter=",", names=True)
    return {name: data[name] for name in data.dtype.names}

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _value_norm(values: np.ndarray):
    """Build a Normalize for the reference/computed value columns.
    Uses TwoSlopeNorm (diverging) when both signs are present, otherwise
    a plain linear norm.
    """
    finite = values[np.isfinite(values)]
    vmin, vmax = float(finite.min()), float(finite.max())
    if vmin < 0 < vmax:
        # diverging: centre at zero
        half = max(abs(vmin), abs(vmax))
        return matplotlib.colors.TwoSlopeNorm(vmin=-half, vcenter=0.0, vmax=half)
    return matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)


def _error_norm(values: np.ndarray, use_log: bool):
    """Build a Normalize for the |error| column."""
    abs_vals = np.abs(values)
    finite   = abs_vals[np.isfinite(abs_vals)]
    pos      = finite[finite > 0]
    vmax     = float(finite.max()) if len(finite) else 1.0
    if use_log:
        lo = float(pos.min()) if len(pos) else 1e-20
        return LogNorm(vmin=max(lo, 1e-20), vmax=max(vmax, 1e-20))
    return matplotlib.colors.Normalize(vmin=0.0, vmax=vmax)


def _scatter(ax, beta, mu, values, title, cmap, norm, label=""):
    """Generic scatter in (β, log10 μ) space."""
    log_mu = np.log10(mu)
    sc = ax.scatter(beta, log_mu, c=values, cmap=cmap, norm=norm,
                    s=18, edgecolors="none")
    cb = plt.colorbar(sc, ax=ax)
    if label:
        cb.set_label(label, fontsize=7)
    ax.set_xlabel("β", fontsize=8)
    ax.set_ylabel("log₁₀(μ)", fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    return sc


def _plot_row(axes_row, d: dict, fcn: str, part: str, metric: str, use_log: bool):
    """Fill one row: [reference, computed, error]."""
    beta   = d["beta"]
    mu     = d["mu"]
    ref    = d["expected"]
    comp   = d["computed"]
    err    = d[metric]
    fail   = ~d["pass"].astype(bool)

    part_label = PART_LABEL[part]

    # Shared value normalisation / colormap for ref and computed
    all_vals  = np.concatenate([ref, comp])
    val_norm  = _value_norm(all_vals)
    val_cmap  = "RdBu_r" if isinstance(val_norm, TwoSlopeNorm) else "viridis"

    # Column 0 – reference values
    _scatter(axes_row[0], beta, mu, ref,
             f"{part_label} — reference", val_cmap, val_norm)

    # Column 1 – computed values (same scale as reference for direct comparison)
    _scatter(axes_row[1], beta, mu, comp,
             f"{part_label} — computed", val_cmap, val_norm)

    # Column 2 – error metric with failing points marked
    err_norm = _error_norm(err, use_log)
    sc = _scatter(axes_row[2], beta, mu, np.abs(err),
                  f"{part_label} — {metric}",
                  "hot_r", err_norm, label=metric)
    if fail.any():
        log_mu = np.log10(mu)
        axes_row[2].scatter(beta[fail], log_mu[fail], marker="x", color="cyan",
                            s=40, linewidths=0.8, zorder=3,
                            label=f"{fail.sum()} fail")
        axes_row[2].legend(fontsize=7, loc="upper right")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Plot time-domain evaluator diff CSVs")
    parser.add_argument("--diff-dir", required=True,
                        help="Directory containing the *_diff.csv files")
    parser.add_argument("--metric", default="abs_err",
                        choices=["abs_err", "rel_err"],
                        help="Error metric to plot (default: abs_err)")
    parser.add_argument("--part", default="all",
                        choices=["total", "G0", "residual", "all"],
                        help="Which decomposition parts to show as rows (default: all)")
    parser.add_argument("--log", action="store_true",
                        help="Log10 colour scale for the error column")
    parser.add_argument("--save", action="store_true",
                        help="Save figures as PNG instead of showing interactively")
    args = parser.parse_args()

    if not os.path.isdir(args.diff_dir):
        sys.exit(f"ERROR: directory not found: {args.diff_dir}")

    parts_to_plot = PARTS if args.part == "all" else [args.part]
    n_rows = len(parts_to_plot)
    n_cols = 3  # reference | computed | error

    scale_label = "log" if args.log else "linear"

    for fcn in FUNCTIONS:
        # Collect data for every requested part first
        rows_data = {}
        for part in parts_to_plot:
            csv_path = os.path.join(args.diff_dir, csv_name(fcn, part))
            if os.path.isfile(csv_path):
                rows_data[part] = load_csv(csv_path)
            else:
                print(f"  WARNING: {csv_path} not found — skipping")

        if not rows_data:
            print(f"No CSV files found for {fcn} in {args.diff_dir}")
            continue

        fig, axes = plt.subplots(n_rows, n_cols,
                                 figsize=(n_cols * 5, n_rows * 3.8 + 0.6),
                                 squeeze=False)

        fig.suptitle(
            f"{fcn}   [{args.metric}, {scale_label} error scale]",
            fontsize=12, fontweight="bold"
        )

        for ri, part in enumerate(parts_to_plot):
            if part not in rows_data:
                for ci in range(n_cols):
                    axes[ri, ci].set_visible(False)
                continue
            _plot_row(axes[ri], rows_data[part], fcn, part, args.metric, args.log)

        fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

        if args.save:
            out = os.path.join(args.diff_dir, f"{fcn}_{args.metric}.png")
            fig.savefig(out, dpi=150, bbox_inches="tight")
            print(f"Saved: {out}")
        else:
            plt.show()


if __name__ == "__main__":
    main()

