#!/usr/bin/env python3
"""Plot cross-correlograms from CSV data exported by MATLAB.

Reads *_data.csv + *_meta.json pairs under output/fig{3,5}/correlogram/
and produces publication-quality PDFs + PNGs via matplotlib. Text and
axes are vector; the line plot is rendered at 600 DPI.

Usage:
    python3 scripts/plot_correlograms.py            # all
    python3 scripts/plot_correlograms.py --study doi # DOI only
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
# Embed text as TrueType so Inkscape can edit all text elements.
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# Correlograms now live in SI (supplementary), not in the main figure dirs.
_SI_CORR_DIR = ROOT / "output" / "SI" / "correlogram"

DPI = 600

# Nature/NPP styling
FONT_FAMILY = "Arial"
TICK_SIZE = 6
LABEL_SIZE = 7
TITLE_SIZE = 7
AXES_LW = 0.6


def plot_correlogram(csv_path: Path, meta_path: Path):
    """Plot one cross-correlogram and save as PDF + PNG."""
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    with open(meta_path) as f:
        meta = json.load(f)

    lags = data["lag_ms"]
    baseline = data["baseline"]
    treatment = data["treatment"]

    bc = meta["baseline_color"]
    tc = meta["treatment_color"]

    fig, ax = plt.subplots(figsize=(2.4, 1.8))

    ax.plot(lags, baseline, color=bc, linewidth=1.2,
            label=meta["baseline_label"])
    ax.plot(lags, treatment, color=tc, linewidth=1.2,
            label=meta["treatment_label"])
    ax.axvline(0, color=(0.6, 0.6, 0.6), linestyle=":", linewidth=0.6)

    ax.set_xlabel("Lag (ms)", fontsize=LABEL_SIZE, fontfamily=FONT_FAMILY)
    ax.set_ylabel("Correlation (z)", fontsize=LABEL_SIZE,
                  fontfamily=FONT_FAMILY)
    ax.set_title(meta["title"], fontsize=TITLE_SIZE, fontweight="bold",
                 fontfamily=FONT_FAMILY)
    ax.set_xlim(-meta["max_lag_ms"], meta["max_lag_ms"])
    ax.legend(fontsize=TICK_SIZE, frameon=False, loc="best")

    ax.tick_params(labelsize=TICK_SIZE, direction="out", width=AXES_LW)
    for spine in ax.spines.values():
        spine.set_linewidth(AXES_LW)

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
    print(f"\n=== {study.upper()} correlograms ===")
    corr_dir = _SI_CORR_DIR
    if not corr_dir.exists():
        print(f"  SKIP — {corr_dir} not found")
        return

    # Filter CSV files by study prefix (e.g. doi_*, ket_*)
    csv_files = sorted(corr_dir.glob(f"{study}_*_data.csv"))
    if not csv_files:
        print(f"  SKIP — no {study}_* CSV files in {corr_dir}")
        return

    for csv_path in csv_files:
        meta_path = csv_path.with_name(
            csv_path.name.replace("_data.csv", "_meta.json"))
        if not meta_path.exists():
            print(f"  SKIP {csv_path.name} — no matching meta JSON")
            continue
        plot_correlogram(csv_path, meta_path)


def main():
    parser = argparse.ArgumentParser(
        description="Plot cross-correlograms from CSV data")
    parser.add_argument("--study", choices=["doi", "ket"],
                        help="Process one study (default: both)")
    args = parser.parse_args()

    studies = [args.study] if args.study else ["doi", "ket"]
    for s in studies:
        process_study(s)

    print("\nDone.")


if __name__ == "__main__":
    main()
