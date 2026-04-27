#!/usr/bin/env python3
"""Plot adjacency heatmaps from CSV data exported by MATLAB.

Reads *_data.csv + *_meta.json pairs under output/fig{3,5}/heatmap/
and produces publication-quality PDFs + PNGs via matplotlib. Text and
axes are vector; the heatmap is rasterized at 600 DPI.

Usage:
    python3 scripts/plot_heatmaps.py            # all
    python3 scripts/plot_heatmaps.py --study doi # DOI only
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# Embed text as TrueType so Inkscape can edit tick labels, titles, etc.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

_FIG_MAP = {"doi": "fig3", "ket": "fig5"}

DPI = 600

# Nature/NPP styling constants
FONT_FAMILY = "Arial"
TICK_SIZE = 6
LABEL_SIZE = 7
TITLE_SIZE = 7
AXES_LW = 0.6


def _diverging_cmap(n=256):
    """Blue-to-red diverging colormap without a white midpoint."""
    colors = np.column_stack([
        np.linspace(0.11, 0.75, n),   # R: blue-low  -> red-high
        np.linspace(0.30, 0.10, n),   # G: falls off both ends
        np.linspace(0.60, 0.15, n),   # B: blue-high -> red-low
    ])
    return mcolors.ListedColormap(colors)


DIVERGING = _diverging_cmap()


def plot_heatmap(csv_path: Path, meta_path: Path):
    """Plot one adjacency heatmap and save as PDF + PNG."""
    A = np.loadtxt(csv_path, delimiter=",")
    with open(meta_path) as f:
        meta = json.load(f)

    symmetric = meta["symmetric"]
    title = meta["title"]
    tick_pos = [int(t) - 1 for t in meta["tick_pos"]]  # 0-indexed
    tick_labels = meta["tick_labels"]
    n_ch = A.shape[0]

    # Thin out tick labels if there are too many to fit legibly.
    # Show every 5th label when >20 channels; otherwise show all.
    if len(tick_pos) > 20:
        step = 5
        tick_pos = tick_pos[::step]
        tick_labels = tick_labels[::step]

    fig, ax = plt.subplots(figsize=(3.6, 3.0))

    if symmetric:
        finite = A[np.isfinite(A)]
        vmax = float(np.max(np.abs(finite))) if len(finite) else 1.0
        if vmax == 0:
            vmax = 1.0
        im = ax.imshow(A, aspect="equal", origin="upper",
                       cmap=DIVERGING, vmin=-vmax, vmax=vmax,
                       rasterized=True)
    else:
        finite = A[np.isfinite(A)]
        vmin = float(np.min(finite)) if len(finite) else 0
        vmax = float(np.max(finite)) if len(finite) else 1
        if vmin == vmax:
            vmin -= 1
            vmax += 1
        im = ax.imshow(A, aspect="equal", origin="upper",
                       cmap="viridis", vmin=vmin, vmax=vmax,
                       rasterized=True)

    ax.set_title(title, fontsize=TITLE_SIZE, fontweight="bold",
                 fontfamily=FONT_FAMILY)
    ax.set_xlabel("Electrode (MCS)", fontsize=LABEL_SIZE,
                  fontfamily=FONT_FAMILY)
    ax.set_ylabel("Electrode (MCS)", fontsize=LABEL_SIZE,
                  fontfamily=FONT_FAMILY)

    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, rotation=45, fontfamily=FONT_FAMILY)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels(tick_labels, fontfamily=FONT_FAMILY)

    ax.tick_params(labelsize=TICK_SIZE, direction="out", width=AXES_LW)
    for spine in ax.spines.values():
        spine.set_linewidth(AXES_LW)

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Peak cross-correlation (z)",
                 fontsize=LABEL_SIZE, fontweight="bold",
                 fontfamily=FONT_FAMILY)
    cb.ax.tick_params(labelsize=TICK_SIZE, width=AXES_LW)
    for tl in cb.ax.get_yticklabels():
        tl.set_fontfamily(FONT_FAMILY)
    cb.outline.set_linewidth(AXES_LW)

    fig.tight_layout()

    stem = csv_path.stem.replace("_data", "")
    out_dir = csv_path.parent

    pdf_path = out_dir / f"{stem}.pdf"
    fig.savefig(pdf_path, format="pdf", dpi=DPI, bbox_inches="tight")

    png_path = out_dir / f"{stem}.png"
    fig.savefig(png_path, format="png", dpi=DPI, bbox_inches="tight")

    plt.close(fig)
    print(f"  {pdf_path.name} + {png_path.name}")


def process_study(study: str):
    print(f"\n=== {study.upper()} heatmaps ===")
    heat_dir = ROOT / "output" / _FIG_MAP[study] / "heatmap"
    if not heat_dir.exists():
        print(f"  SKIP - {heat_dir} not found")
        return

    csv_files = sorted(heat_dir.glob("*_data.csv"))
    if not csv_files:
        print(f"  SKIP - no CSV files in {heat_dir}")
        return

    for csv_path in csv_files:
        meta_path = csv_path.with_name(
            csv_path.name.replace("_data.csv", "_meta.json"))
        if not meta_path.exists():
            print(f"  SKIP {csv_path.name} - no matching meta JSON")
            continue
        plot_heatmap(csv_path, meta_path)


def main():
    parser = argparse.ArgumentParser(
        description="Plot adjacency heatmaps from CSV data")
    parser.add_argument("--study", choices=["doi", "ket"],
                        help="Process one study (default: both)")
    args = parser.parse_args()

    studies = [args.study] if args.study else ["doi", "ket"]
    for s in studies:
        process_study(s)

    print("\nDone.")


if __name__ == "__main__":
    main()
