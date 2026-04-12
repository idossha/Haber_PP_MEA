#!/usr/bin/env python3
"""aggregate_panels.py - assemble individual figure PNGs into composite SVGs.

Each composite is specified by an in-source LAYOUTS dict with:
  - overall page size (mm)
  - one Panel entry per rectangle, in mm, with:
      .file   : PNG filename under figures_out/
      .rect   : (x, y, w, h) in mm
      .label  : (text, x, y, font_size) for the panel letter (or None)

Run:

    python3 scripts/aggregate_panels.py            # build all
    python3 scripts/aggregate_panels.py --list     # list available composites
    python3 scripts/aggregate_panels.py Figure2_DOI_rates   # build one

Output goes to paper/figures/<name>.svg. The PNG content is base64-
embedded so the resulting SVG is a single self-contained file.
"""
from __future__ import annotations
import argparse
import base64
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

REPO_ROOT   = Path(__file__).resolve().parent.parent
FIGURES_IN  = REPO_ROOT / "figures_out"
FIGURES_OUT = REPO_ROOT / "paper" / "figures"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class Panel:
    file:  str
    rect:  Tuple[float, float, float, float]   # (x, y, w, h) in mm
    label: Optional[Tuple[str, float, float, float]] = None
    # label = (text, x_mm, y_mm, font_size_mm)


@dataclass
class Layout:
    name:    str
    title:   str
    size_mm: Tuple[float, float]     # (width, height)
    panels:  List[Panel] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Composite layouts — edit these to change page geometry or add figures
# ---------------------------------------------------------------------------

def _layout_rates_composite(study: str, fig_num: int, title: str) -> Layout:
    """3-column rates composite: A traces / B raster wide-left, C/D spike+burst
    paired rates narrow-middle, E/F channel-category + percent-change violins
    wide-right."""
    W, H = 340, 220
    m   = 12.0   # top margin
    ls  = 8.0    # label inset
    lsz = 6.0    # label font size mm

    # wide-left block (traces + raster) ~42 % width
    leftX, leftW = 10, 140
    aY, aH = m,            90
    bY, bH = m + aH + 6,   95

    # middle narrow column (~18 % width)
    midX, midW = leftX + leftW + 10, 56
    cY, cH = m,            90
    dY, dH = m + aH + 6,   95

    # right wide column (~34 % width) split into two sub-panels (cats + violin)
    rightX      = midX + midW + 10
    rightW      = W - rightX - 10
    subGap      = 4
    subW        = (rightW - subGap) / 2
    eLeftX      = rightX
    eRightX     = rightX + subW + subGap
    eY, eH      = m,            90
    fY, fH      = m + aH + 6,   95

    panels = [
        # Panel A - traces
        Panel(file=f"{study}_rate_panel_A_traces.png",
              rect=(leftX, aY, leftW, aH),
              label=("A.", leftX - ls, aY - 2, lsz)),
        # Panel B - raster
        Panel(file=f"{study}_rate_panel_B_raster.png",
              rect=(leftX, bY, leftW, bH),
              label=("B.", leftX - ls, bY - 2, lsz)),
        # Panel C - spike rate paired (single dataset view)
        Panel(file=f"{study}_spike_rate_paired.png",
              rect=(midX, cY, midW, cH),
              label=("C.", midX - ls, cY - 2, lsz)),
        # Panel D - burst rate paired
        Panel(file=f"{study}_burst_rate_paired.png",
              rect=(midX, dY, midW, dH),
              label=("D.", midX - ls, dY - 2, lsz)),
        # Panel E - spike rates across all datasets (cats + violin)
        Panel(file=f"{study}_spike_rate_change_categories.png",
              rect=(eLeftX, eY, subW, eH),
              label=("E.", eLeftX - ls, eY - 2, lsz)),
        Panel(file=f"{study}_spike_pct_change_violin.png",
              rect=(eRightX, eY, subW, eH),
              label=None),
        # Panel F - burst rates across all datasets
        Panel(file=f"{study}_burst_rate_change_categories.png",
              rect=(eLeftX, fY, subW, fH),
              label=("F.", eLeftX - ls, fY - 2, lsz)),
        Panel(file=f"{study}_burst_pct_change_violin.png",
              rect=(eRightX, fY, subW, fH),
              label=None),
    ]
    return Layout(name=f"Figure{fig_num}_{study.upper()}_rates",
                  title=title, size_mm=(W, H), panels=panels)


def _layout_new_angles_composite() -> Layout:
    """Figure 6 - Information flow and criticality composite.

    Two-row (DOI top, ket bottom) x five-column (A..E) grid of paired-slope
    and distribution panels for the four exploratory new-angle metrics plus
    an avalanche-classification table. All trend-level panels; see
    paper/figures/Figure6_NewAngles_caption.md for the honest "N of M wells"
    language and the pointer to the Limitations paragraph.

    JNE double-column tall format. Panel widths compressed to 31 mm (vs the
    32 mm target) so that 5 columns + 4*3 mm gutters fit inside the 168 mm
    content strip on a 176 mm canvas. Row height 92 mm (vs 85 mm target) so
    that both rows + 12 mm inter-row gap + 8 mm top margin + 15 mm bottom
    margin sum to 219 mm, just inside the 220 mm canvas ceiling.
    """
    # -- canvas --------------------------------------------------------
    W, H = 176, 220               # mm; <=178 width, <=230 height (JNE)

    # -- grid step values ---------------------------------------------
    panel_w   = 31.0              # 5*31 + 4*3 = 167 mm content strip
    panel_h   = 92.0              # 8 + 92 + 12 + 92 + 15 = 219 mm tall stack
    col_gap   = 3.0               # inter-column gutter
    row_gap   = 12.0              # inter-row gutter (holds row label)
    margin_l  = 5.0               # left margin
    margin_t  = 8.0               # top margin (holds "Figure 6" title)

    # -- label style ---------------------------------------------------
    ls        = 6.0               # label inset (mm) left of panel rect
    lsz       = 5.0               # label font size (mm)

    # -- column x origins ---------------------------------------------
    col_x = [margin_l + i * (panel_w + col_gap) for i in range(5)]
    # => [5.0, 39.0, 73.0, 107.0, 141.0]

    # -- row y origins ------------------------------------------------
    row1_y = margin_t                                      # 8.0
    row2_y = margin_t + panel_h + row_gap                  # 8 + 92 + 12 = 112.0

    letters = ["A", "B", "C", "D", "E"]
    panels: List[Panel] = []

    # --- top row: DOI arm, n = 6 paired wells -----------------------
    for i, L in enumerate(letters):
        panels.append(Panel(
            file=f"doi_new_angles_{L}_{_new_angles_slug(L)}.png",
            rect=(col_x[i], row1_y, panel_w, panel_h),
            label=(f"{L}.", col_x[i] - ls, row1_y - 1.0, lsz),
        ))

    # --- bottom row: ketanserin arm, n = 3 paired wells --------------
    # Uses A' .. E' labels so reviewers don't mistake them for the DOI row.
    for i, L in enumerate(letters):
        panels.append(Panel(
            file=f"ket_new_angles_{L}_{_new_angles_slug(L)}.png",
            rect=(col_x[i], row2_y, panel_w, panel_h),
            label=(f"{L}\u2032.", col_x[i] - ls, row2_y - 1.0, lsz),
        ))

    return Layout(
        name="Figure6_NewAngles",
        title="Figure 6 - Information flow, burst synchrony, entropy, "
              "and avalanches (exploratory)",
        size_mm=(W, H),
        panels=panels,
    )


def _new_angles_slug(letter: str) -> str:
    return {
        "A": "avalanches_size",
        "B": "burstsync",
        "C": "entropy",
        "D": "te",
        "E": "avalanche_class",
    }[letter]


def _layout_connectivity_composite(study: str, figNum: int, title: str) -> Layout:
    W, H = 220, 260
    ls, lsz = 8.0, 6.0
    return Layout(
        name=f"Figure{figNum}_{study.upper()}_connectivity",
        title=title,
        size_mm=(W, H),
        panels=[
            Panel(file=f"{study}_connectivity_exemplar.png",
                  rect=(10, 12, W - 20, 120),
                  label=("A.", 10 - ls, 10, lsz)),
            Panel(file=f"{study}_connectivity_summary.png",
                  rect=(10, 140, W - 20, 110),
                  label=("B.", 10 - ls, 138, lsz)),
        ],
    )


LAYOUTS: Dict[str, Layout] = {
    "Figure2_DOI_rates": _layout_rates_composite(
        "doi", 2,
        "Figure 2 - DOI effects on cortical network activity"),
    "Figure3_KET_rates": _layout_rates_composite(
        "ket", 3,
        "Figure 3 - DOI + Ketanserin effects on cortical network activity"),
    "Figure4_DOI_connectivity": _layout_connectivity_composite(
        "doi", 4, "Figure 4 - Functional connectivity under DOI"),
    "Figure5_KET_connectivity": _layout_connectivity_composite(
        "ket", 5, "Figure 5 - Functional connectivity under DOI + Ketanserin"),
    "Figure6_NewAngles": _layout_new_angles_composite(),
}

# Rename helper: "Figure2_DOI_rates" but we also register the layout under
# "Figure2" pointing at the DOI composite (alias) for convenience.


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def png_dimensions(path: Path) -> Tuple[int, int]:
    with path.open("rb") as fh:
        signature = fh.read(8)
        if signature != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"{path} is not a PNG file")
        fh.read(4)
        chunk_type = fh.read(4)
        if chunk_type != b"IHDR":
            raise ValueError(f"{path} missing IHDR chunk")
        ihdr = fh.read(8)
        width, height = struct.unpack(">II", ihdr)
        return width, height


def png_to_data_uri(path: Path) -> str:
    blob = path.read_bytes()
    b64 = base64.b64encode(blob).decode("ascii")
    return f"data:image/png;base64,{b64}"


def panel_fit_rect(panel_w_px: int, panel_h_px: int,
                   box_x: float, box_y: float,
                   box_w: float, box_h: float
                   ) -> Tuple[float, float, float, float]:
    if panel_w_px == 0 or panel_h_px == 0:
        return box_x, box_y, box_w, box_h
    aspect_panel = panel_w_px / panel_h_px
    aspect_box   = box_w / box_h
    if aspect_panel > aspect_box:
        w = box_w
        h = box_w / aspect_panel
    else:
        h = box_h
        w = box_h * aspect_panel
    x = box_x + (box_w - w) / 2
    y = box_y + (box_h - h) / 2
    return x, y, w, h


def build_svg(layout: Layout) -> str:
    W, H = layout.size_mm
    parts = []
    parts.append(
        f'<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{W}mm" height="{H}mm" '
        f'viewBox="0 0 {W} {H}" version="1.1">\n'
    )
    parts.append(
        f'  <title>{layout.title}</title>\n'
        f'  <rect x="0" y="0" width="{W}" height="{H}" fill="white"/>\n'
        f'  <text x="{W/2}" y="6" font-family="Arial, sans-serif" '
        f'font-size="4" font-weight="bold" text-anchor="middle">'
        f'{layout.title}</text>\n'
    )
    for panel in layout.panels:
        src = FIGURES_IN / panel.file
        if not src.exists():
            print(f"  [warn] missing panel: {panel.file}", file=sys.stderr)
            continue
        try:
            pw_px, ph_px = png_dimensions(src)
        except Exception as exc:
            print(f"  [warn] {panel.file}: {exc}", file=sys.stderr)
            pw_px, ph_px = 1000, 1000
        bx, by, bw, bh = panel.rect
        px, py, pw, ph = panel_fit_rect(pw_px, ph_px, bx, by, bw, bh)
        data_uri = png_to_data_uri(src)
        parts.append(f'  <g id="panel_{panel.file}">\n')
        if panel.label:
            text, lx, ly, lsz = panel.label
            parts.append(
                f'    <text x="{lx}" y="{ly}" '
                f'font-family="Arial, sans-serif" font-size="{lsz}" '
                f'font-weight="bold">{text}</text>\n'
            )
        parts.append(
            f'    <image x="{px}" y="{py}" width="{pw}" height="{ph}" '
            f'xlink:href="{data_uri}" '
            f'preserveAspectRatio="xMidYMid meet"/>\n'
        )
        parts.append('  </g>\n')
    parts.append("</svg>\n")
    return "".join(parts)


def build_composite(name: str) -> Path:
    if name not in LAYOUTS:
        raise KeyError(f"unknown layout: {name}")
    layout = LAYOUTS[name]
    FIGURES_OUT.mkdir(parents=True, exist_ok=True)
    svg_text = build_svg(layout)
    out_path = FIGURES_OUT / f"{layout.name}.svg"
    out_path.write_text(svg_text, encoding="utf-8")
    return out_path


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("names", nargs="*",
                    help="composite layout names (default: all)")
    ap.add_argument("--list", action="store_true",
                    help="list available composites and exit")
    args = ap.parse_args(argv)

    if args.list:
        for name, layout in LAYOUTS.items():
            print(f"  {name}  -  {layout.title}")
        return 0

    names = args.names if args.names else list(LAYOUTS.keys())
    rc = 0
    for name in names:
        try:
            out = build_composite(name)
            size_kb = out.stat().st_size / 1024
            print(f"  [ok] {out.relative_to(REPO_ROOT)}  ({size_kb:.1f} KB)")
        except KeyError as exc:
            print(f"  [err] {exc}", file=sys.stderr)
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
