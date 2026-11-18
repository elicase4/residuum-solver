"""Plot the Q1 h-convergence study written by Benchmark_Equation_Heateq_ConvergenceStudy.

Usage:
    python plot_convergence.py <path/to/HeatEq_ConvergenceStudy_Q1.csv> [output_prefix]

Writes <prefix>_errors.png (L2 and H1-seminorm error vs h, with reference slopes and the
observed rate fitted over the finest three meshes) and <prefix>_cg.png (CG iterations vs h).
"""

import csv
import math
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAVY = "#1B3A6B"
BLUE = "#1F5FAD"
GRAY = "#8FA3BA"
BODY = "#2E4A70"
GRID = "#E3E8EE"


def read(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {
        "h": [float(r["h"]) for r in rows],
        "L2": [float(r["error_L2"]) for r in rows],
        "H1": [float(r["error_H1_semi"]) for r in rows],
        "it": [int(r["cg_iterations"]) for r in rows],
        "dofs": [int(r["free_dofs"]) for r in rows],
    }


def fitted_rate(h, e, last=3):
    """Least-squares slope of log(e) vs log(h) over the finest `last` meshes."""
    xs = [math.log(v) for v in h[-last:]]
    ys = [math.log(v) for v in e[-last:]]
    xm, ym = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sum((x - xm) ** 2 for x in xs)


def style(ax):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", color=GRID, lw=1)
    ax.grid(False, which="minor")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#C9D1DB")
    ax.tick_params(colors=BODY, labelsize=13)
    ax.invert_xaxis()  # refinement reads left to right


def slope_triangle(ax, h, e, order, color, below=True):
    """Small reference triangle of the given order anchored near the fine end of a series."""
    h1, h0 = h[-2], h[-3]
    y0 = e[-3] * (0.22 if below else 3.2)
    y1 = y0 * (h1 / h0) ** order
    ax.plot([h0, h1, h0, h0], [y0, y1, y1, y0], color=color, lw=1.4)
    ax.text(h0 * 1.06, math.sqrt(y0 * y1), f"{order}", color=color, fontsize=13, ha="right", va="center")


def main():
    path = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "convergence_q1"
    d = read(path)
    rL2, rH1 = fitted_rate(d["h"], d["L2"]), fitted_rate(d["h"], d["H1"])

    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": BODY, "axes.labelcolor": BODY})

    # --- error vs h ---
    fig, ax = plt.subplots(figsize=(8, 5.6), dpi=200)
    style(ax)
    ax.plot(d["h"], d["L2"], "o-", color=BLUE, lw=2.5, ms=7, label=f"L2 error  (observed rate {rL2:.2f})")
    ax.plot(d["h"], d["H1"], "s-", color=NAVY, lw=2.5, ms=7, label=f"H1 seminorm error  (observed rate {rH1:.2f})")
    slope_triangle(ax, d["h"], d["L2"], 2, BLUE, below=True)
    slope_triangle(ax, d["h"], d["H1"], 1, NAVY, below=False)
    ax.set_xlabel("element size h", fontsize=14)
    ax.set_ylabel("error", fontsize=14)
    ax.set_title("Q1 elements match theory: O(h²) in L2, O(h) in H1", fontsize=15, color=NAVY, loc="left")
    ax.legend(frameon=False, fontsize=12, loc="lower left")
    fig.tight_layout()
    fig.savefig(f"{prefix}_errors.png", facecolor="white")

    # --- CG iterations vs h ---
    fig, ax = plt.subplots(figsize=(8, 5.6), dpi=200)
    style(ax)
    rIt = -fitted_rate(d["h"], [float(i) for i in d["it"]])
    ax.plot(d["h"], d["it"], "o-", color=BLUE, lw=2.5, ms=7, label=f"Jacobi-preconditioned CG  (fitted: iterations \u221d h^-{rIt:.2f})")
    ax.set_xlabel("element size h", fontsize=14)
    ax.set_ylabel("CG iterations to 1e-13 relative residual", fontsize=14)
    ax.set_title("Iterations double with each refinement", fontsize=15, color=NAVY, loc="left")
    ax.legend(frameon=False, fontsize=12, loc="upper left")
    fig.tight_layout()
    fig.savefig(f"{prefix}_cg.png", facecolor="white")

    print(f"L2 rate {rL2:.3f}, H1 rate {rH1:.3f}, CG iterations ~ h^-{rIt:.2f}")


if __name__ == "__main__":
    main()
