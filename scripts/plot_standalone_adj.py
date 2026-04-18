#!/usr/bin/env python3
"""Plot standalone adjacency heatmap SVGs from CSV matrices.

Reads adjacency CSVs (exported by gen_standalone_adj.m) and produces
publication-quality SVGs via matplotlib. Text and axes are vector;
the heatmap and colorbar are rasterized at 600 DPI. The output SVGs
import cleanly into Inkscape.

Usage:
    python3 scripts/plot_standalone_adj.py            # all
    python3 scripts/plot_standalone_adj.py --study doi # DOI only
"""
import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Figure mapping: study -> figure directory
_FIG_MAP = {"doi": "fig3", "ket": "fig5"}

DPI = 600

# Nature/NPP styling constants
FONT_FAMILY = "Arial"
TICK_SIZE = 6
LABEL_SIZE = 7
TITLE_SIZE = 7
AXES_LW = 0.6

# Diverging colormap (blue-white-red) matching MATLAB's diverging_cmap
def _diverging_cmap(n=256):
    half = n // 2
    top = n - half
    blue = np.column_stack([
        np.linspace(0.11, 1, half),
        np.linspace(0.30, 1, half),
        np.linspace(0.60, 1, half),
    ])
    red = np.column_stack([
        np.linspace(1, 0.75, top),
        np.linspace(1, 0.10, top),
        np.linspace(1, 0.15, top),
    ])
    colors = np.vstack([blue, red])
    return mcolors.ListedColormap(colors)

DIVERGING = _diverging_cmap()


def plot_adj(A: np.ndarray, title: str, symmetric: bool, out_path: Path):
    """Plot one adjacency heatmap and save as SVG + PNG."""
    fig, ax = plt.subplots(figsize=(3.2, 2.8))

    if symmetric:
        finite = A[np.isfinite(A)]
        vmax = np.max(np.abs(finite)) if len(finite) else 1.0
        if vmax == 0:
            vmax = 1.0
        im = ax.imshow(A, aspect="equal", origin="upper",
                       cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                       rasterized=True)
    else:
        finite = A[np.isfinite(A)]
        vmin = np.min(finite) if len(finite) else 0
        vmax = np.max(finite) if len(finite) else 1
        if vmin == vmax:
            vmin -= 1; vmax += 1
        im = ax.imshow(A, aspect="equal", origin="upper",
                       cmap="parula" if "parula" in plt.colormaps() else "viridis",
                       vmin=vmin, vmax=vmax,
                       rasterized=True)

    ax.set_title(title, fontsize=TITLE_SIZE, fontweight="bold", fontfamily=FONT_FAMILY)
    ax.set_xlabel("Channel index", fontsize=LABEL_SIZE, fontfamily=FONT_FAMILY)
    ax.set_ylabel("Channel index", fontsize=LABEL_SIZE, fontfamily=FONT_FAMILY)
    ax.tick_params(labelsize=TICK_SIZE, direction="out", width=AXES_LW)
    for spine in ax.spines.values():
        spine.set_linewidth(AXES_LW)

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Peak cross-correlation (z)",
                 fontsize=LABEL_SIZE, fontweight="bold", fontfamily=FONT_FAMILY)
    cb.ax.tick_params(labelsize=TICK_SIZE, width=AXES_LW)
    cb.outline.set_linewidth(AXES_LW)

    fig.tight_layout()

    # SVG with rasterized heatmap, vector text/axes
    svg_path = out_path.with_suffix(".svg")
    fig.savefig(svg_path, format="svg", dpi=DPI, bbox_inches="tight")

    # PNG for quick preview
    png_path = out_path.with_suffix(".png")
    fig.savefig(png_path, format="png", dpi=DPI, bbox_inches="tight")

    plt.close(fig)
    svg_kb = svg_path.stat().st_size / 1024
    print(f"  {svg_path.name} ({svg_kb:.0f} KB) + {png_path.name}")


PANELS = [
    ("baseline",  "Baseline",                     False),
    ("treatment", None,                            False),  # title filled per-study
    ("delta",     "\u0394 (treatment \u2212 baseline)", True),
]

STUDY_LABELS = {
    "doi": {"treatment": "DOI"},
    "ket": {"treatment": "Ketanserin"},
}


def process_study(study: str):
    print(f"\n=== {study.upper()} ===")
    labels = STUDY_LABELS[study]
    panel_dir = ROOT / "output" / _FIG_MAP[study] / "panels"

    for suffix, title, symmetric in PANELS:
        csv = panel_dir / f"{study}_adj_{suffix}.csv"
        if not csv.exists():
            print(f"  SKIP {csv.name} — not found")
            continue

        A = np.loadtxt(csv, delimiter=",")
        if title is None:
            title = labels.get(suffix, suffix.capitalize())

        out = panel_dir / f"{study}_adj_{suffix}"
        plot_adj(A, title, symmetric, out)


def main():
    parser = argparse.ArgumentParser(description="Plot adjacency heatmap SVGs")
    parser.add_argument("--study", choices=["doi", "ket"],
                        help="Process one study (default: both)")
    args = parser.parse_args()

    studies = [args.study] if args.study else ["doi", "ket"]
    for s in studies:
        process_study(s)

    print("\nDone.")


if __name__ == "__main__":
    main()
