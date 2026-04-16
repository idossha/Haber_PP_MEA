#!/usr/bin/env python3
"""Compose publication-ready composite figures with inline SVG panel content.

Each panel's vector content is embedded as a <g> group (Inkscape layer),
NOT as base64 <image> data URIs. Every path, text element, and shape from
each panel remains individually selectable and editable in Inkscape.

Pipeline steps (all run by this script):
  1. Convert panel PDFs → SVGs via Inkscape  (figures/panels/ → figures/_panel_svgs/)
  2. Inline panel SVGs into composite SVGs    (→ figures/composite_svgs/)
  3. Export composites → PDF + PNG via Inkscape (→ figures/composite_pdfs/, composite_pngs/)

Usage:
    python3 scripts/compose_svg_inline.py                  # all figures
    python3 scripts/compose_svg_inline.py --fig 2          # Figure 2 only
    python3 scripts/compose_svg_inline.py --skip-convert    # skip PDF→SVG step
    python3 scripts/compose_svg_inline.py --skip-export     # skip SVG→PDF/PNG step
"""
import argparse
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
PANELS_PDF = ROOT / "figures" / "panels"
PANELS_SVG = ROOT / "figures" / "_panel_svgs"
PANELS_PNG = ROOT / "figures" / "_panel_pngs"
OUT_SVG = ROOT / "figures" / "composite_svgs"
OUT_PDF = ROOT / "figures" / "composite_pdfs"
OUT_PNG = ROOT / "figures" / "composite_pngs"

INKSCAPE = "/Applications/Inkscape.app/Contents/MacOS/inkscape"

# Composite canvas in mm (JNE double-column = 180 mm)
PAGE_W = 183.0  # NPP double-column width

# Layout constants
LEFT_FRAC = 0.42
GAP_MM = 2.0
PANEL_GAP = 1.5
TOP_MARGIN = 5.0
BOTTOM_MARGIN = 4.0

# Panel label style (bold, 8 pt, lowercase per NPP specs)
LABEL_FONT = "Arial"
LABEL_SIZE_PT = 8
LABEL_SIZE_MM = LABEL_SIZE_PT * 0.3528  # pt → mm

DPI = 600

# ── Namespace handling ───────────────────────────────────────────────────
NS = {
    "svg": "http://www.w3.org/2000/svg",
    "xlink": "http://www.w3.org/1999/xlink",
    "inkscape": "http://www.inkscape.org/namespaces/inkscape",
    "sodipodi": "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd",
}
for prefix, uri in NS.items():
    ET.register_namespace(prefix, uri)
ET.register_namespace("", NS["svg"])


# ── Figure definitions ───────────────────────────────────────────────────

def rate_panels(prefix):
    """Panel file map for a rate figure (Figs 2/3)."""
    return {
        "traces":      f"{prefix}_rate_panel_A_traces",
        "raster":      f"{prefix}_rate_panel_B_raster",
        "scatter_spike": f"{prefix}_spike_rate_scatter",
        "scatter_burst": f"{prefix}_burst_rate_scatter",
        "cats_spike":  f"{prefix}_spike_rate_change_categories",
        "cats_burst":  f"{prefix}_burst_rate_change_categories",
        "violin_spike": f"{prefix}_spike_pct_change_violin",
        "violin_burst": f"{prefix}_burst_pct_change_violin",
    }

def conn_panels(prefix):
    """Panel file map for a connectivity figure (Figs 4/5)."""
    return {
        "exemplar": f"{prefix}_connectivity_exemplar",
        "summary":  f"{prefix}_connectivity_summary",
    }

FIGURES = {
    2: {"type": "rate", "prefix": "doi", "title": "Figure 2 — DOI firing and burst rates"},
    3: {"type": "rate", "prefix": "ket", "title": "Figure 3 — Ketanserin firing and burst rates"},
    4: {"type": "conn", "prefix": "doi", "title": "Figure 4 — DOI functional connectivity"},
    5: {"type": "conn", "prefix": "ket", "title": "Figure 5 — Ketanserin functional connectivity"},
}


# ── Step 1: PDF → SVG conversion ────────────────────────────────────────

def convert_pdfs_to_svgs(fig_nums):
    """Convert panel PDFs to SVGs using Inkscape."""
    # Collect all needed panel basenames
    needed = set()
    for n in fig_nums:
        info = FIGURES[n]
        panels = rate_panels(info["prefix"]) if info["type"] == "rate" else conn_panels(info["prefix"])
        needed.update(panels.values())

    PANELS_SVG.mkdir(parents=True, exist_ok=True)
    PANELS_PNG.mkdir(parents=True, exist_ok=True)

    converted = 0
    for basename in sorted(needed):
        pdf = PANELS_PDF / f"{basename}.pdf"
        svg = PANELS_SVG / f"{basename}.svg"
        png = PANELS_PNG / f"{basename}.png"

        if not pdf.exists():
            print(f"  WARNING: {pdf.name} not found — skipping")
            continue

        # Convert PDF → SVG
        if not svg.exists() or pdf.stat().st_mtime > svg.stat().st_mtime:
            print(f"  Converting {pdf.name} → SVG")
            subprocess.run([
                INKSCAPE, str(pdf),
                "--export-type=svg",
                f"--export-filename={svg}",
            ], capture_output=True, check=True)
            converted += 1
        else:
            print(f"  {svg.name} up to date")

        # Also copy PNG if it exists in panels/
        src_png = PANELS_PDF / f"{basename}.png"
        if src_png.exists() and (not png.exists() or src_png.stat().st_mtime > png.stat().st_mtime):
            import shutil
            shutil.copy2(src_png, png)

    print(f"  PDF→SVG: {converted} converted")


# ── SVG parsing and composition helpers ──────────────────────────────────

def parse_panel(svg_path: Path):
    """Parse a panel SVG, return (viewBox_w, viewBox_h, children, defs)."""
    tree = ET.parse(svg_path)
    root = tree.getroot()

    vb = root.get("viewBox")
    if vb:
        parts = vb.split()
        vb_w, vb_h = float(parts[2]), float(parts[3])
    else:
        vb_w = float(root.get("width", "100"))
        vb_h = float(root.get("height", "100"))

    children = []
    defs = []
    for child in root:
        tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
        if tag == "defs":
            defs.append(child)
        elif tag in ("namedview", "metadata"):
            continue
        else:
            children.append(child)

    return vb_w, vb_h, children, defs


def make_panel_group(svg_path: Path, panel_id: str, x: float, y: float,
                     target_w: float, target_h: float = None):
    """Create a <g> element containing the panel content, scaled/translated to fit."""
    vb_w, vb_h, children, defs = parse_panel(svg_path)

    sx = target_w / vb_w
    if target_h is not None:
        sy = target_h / vb_h
        s = min(sx, sy)
    else:
        s = sx

    actual_h = vb_h * s

    g = ET.Element("g", attrib={
        "id": panel_id,
        "transform": f"translate({x:.2f},{y:.2f}) scale({s:.6f})",
        "{http://www.inkscape.org/namespaces/inkscape}label": panel_id,
        "{http://www.inkscape.org/namespaces/inkscape}groupmode": "layer",
    })
    for child in children:
        g.append(child)

    return g, defs, actual_h


def make_label(text: str, x: float, y: float):
    """Create a bold panel label text element."""
    t = ET.Element("text", attrib={
        "x": f"{x:.2f}",
        "y": f"{y:.2f}",
        "font-family": LABEL_FONT,
        "font-size": f"{LABEL_SIZE_MM:.2f}",
        "font-weight": "bold",
        "fill": "black",
    })
    t.text = text
    return t


def assemble_svg(groups, labels, all_defs, title_text, page_h):
    """Build the root <svg> element from composed parts."""
    svg_attrs = {
        "xmlns": NS["svg"],
        "xmlns:xlink": NS["xlink"],
        "xmlns:inkscape": NS["inkscape"],
        "xmlns:sodipodi": NS["sodipodi"],
        "width": f"{PAGE_W}mm",
        "height": f"{page_h:.1f}mm",
        "viewBox": f"0 0 {PAGE_W} {page_h:.1f}",
        "version": "1.1",
    }
    root = ET.Element("svg", attrib=svg_attrs)

    title = ET.SubElement(root, "title")
    title.text = title_text

    ET.SubElement(root, "rect", attrib={
        "x": "0", "y": "0",
        "width": str(PAGE_W), "height": f"{page_h:.1f}",
        "fill": "white",
    })

    # Merge defs
    if all_defs:
        merged = ET.SubElement(root, "defs")
        seen_ids = set()
        for d in all_defs:
            for child in d:
                cid = child.get("id", "")
                if cid and cid in seen_ids:
                    continue
                seen_ids.add(cid)
                merged.append(child)

    for g in groups:
        root.append(g)

    label_layer = ET.SubElement(root, "g", attrib={
        "id": "panel_labels",
        "{http://www.inkscape.org/namespaces/inkscape}label": "Panel Labels",
        "{http://www.inkscape.org/namespaces/inkscape}groupmode": "layer",
    })
    for lbl in labels:
        label_layer.append(lbl)

    return root


def write_svg(root, out_path):
    """Write an ElementTree SVG to disk."""
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        tree.write(f, encoding="utf-8", xml_declaration=True)
    size_mb = out_path.stat().st_size / 1e6
    print(f"  Written: {out_path.name} ({size_mb:.1f} MB)")


# ── Step 2: Compose rate figures (Figs 2, 3) ────────────────────────────

def compose_rate(fig_num, prefix, title_text):
    """Compose a rate figure: left col A/B, right col C/D/E over F/G/H."""
    panels = rate_panels(prefix)
    svg_files = {k: PANELS_SVG / f"{v}.svg" for k, v in panels.items()}

    missing = [k for k, v in svg_files.items() if not v.exists()]
    if missing:
        print(f"  SKIP Figure {fig_num}: missing panels {missing}")
        return None

    left_w = PAGE_W * LEFT_FRAC
    right_w = PAGE_W - left_w - GAP_MM
    right_x = left_w + GAP_MM
    rpanel_w = (right_w - 2 * PANEL_GAP) / 3

    lbl_ox = 1.0
    lbl_oy = -1.0

    groups, labels, all_defs = [], [], []

    # Left column: A (traces) over B (raster)
    y = TOP_MARGIN
    gA, dA, hA = make_panel_group(svg_files["traces"], "panel_a_traces", 0, y, left_w)
    groups.append(gA); all_defs.extend(dA)
    labels.append(make_label("a", lbl_ox, y + lbl_oy + LABEL_SIZE_MM))

    y_B = y + hA + GAP_MM
    gB, dB, hB = make_panel_group(svg_files["raster"], "panel_b_raster", 0, y_B, left_w)
    groups.append(gB); all_defs.extend(dB)
    labels.append(make_label("b", lbl_ox, y_B + lbl_oy + LABEL_SIZE_MM))

    left_h = y + hA + GAP_MM + hB

    # Right column: 2 rows × 3 panels, matched to left height
    row_h = (left_h - TOP_MARGIN - GAP_MM) / 2

    for i, (key, pid, lbl) in enumerate([
        ("scatter_spike", "panel_c_scatter_spike", "c"),
        ("cats_spike",    "panel_d_cats_spike",    "d"),
        ("violin_spike",  "panel_e_violin_spike",  "e"),
    ]):
        px = right_x + i * (rpanel_w + PANEL_GAP)
        g, d, h = make_panel_group(svg_files[key], pid, px, TOP_MARGIN, rpanel_w, row_h)
        groups.append(g); all_defs.extend(d)
        labels.append(make_label(lbl, px + lbl_ox, TOP_MARGIN + lbl_oy + LABEL_SIZE_MM))

    y_bot = TOP_MARGIN + row_h + GAP_MM
    for i, (key, pid, lbl) in enumerate([
        ("scatter_burst", "panel_f_scatter_burst", "f"),
        ("cats_burst",    "panel_g_cats_burst",    "g"),
        ("violin_burst",  "panel_h_violin_burst",  "h"),
    ]):
        px = right_x + i * (rpanel_w + PANEL_GAP)
        g, d, h = make_panel_group(svg_files[key], pid, px, y_bot, rpanel_w, row_h)
        groups.append(g); all_defs.extend(d)
        labels.append(make_label(lbl, px + lbl_ox, y_bot + lbl_oy + LABEL_SIZE_MM))

    page_h = max(left_h, y_bot + row_h) + BOTTOM_MARGIN
    return assemble_svg(groups, labels, all_defs, title_text, page_h)


# ── Step 2b: Compose connectivity figures (Figs 4, 5) ───────────────────

def compose_conn(fig_num, prefix, title_text):
    """Compose a connectivity figure: A (exemplar) over B (summary), full width."""
    panels = conn_panels(prefix)
    svg_files = {k: PANELS_SVG / f"{v}.svg" for k, v in panels.items()}

    missing = [k for k, v in svg_files.items() if not v.exists()]
    if missing:
        print(f"  SKIP Figure {fig_num}: missing panels {missing}")
        return None

    lbl_ox = 1.0
    lbl_oy = -1.0

    groups, labels, all_defs = [], [], []

    y = TOP_MARGIN
    gA, dA, hA = make_panel_group(svg_files["exemplar"], "panel_a_exemplar", 0, y, PAGE_W)
    groups.append(gA); all_defs.extend(dA)
    labels.append(make_label("a", lbl_ox, y + lbl_oy + LABEL_SIZE_MM))

    y_B = y + hA + GAP_MM
    gB, dB, hB = make_panel_group(svg_files["summary"], "panel_b_summary", 0, y_B, PAGE_W)
    groups.append(gB); all_defs.extend(dB)
    labels.append(make_label("b", lbl_ox, y_B + lbl_oy + LABEL_SIZE_MM))

    page_h = y_B + hB + BOTTOM_MARGIN
    return assemble_svg(groups, labels, all_defs, title_text, page_h)


# ── Step 3: Export SVG → PDF + PNG ───────────────────────────────────────

def export_composites(fig_nums):
    """Export composite SVGs to PDF and PNG via Inkscape."""
    OUT_PDF.mkdir(parents=True, exist_ok=True)
    OUT_PNG.mkdir(parents=True, exist_ok=True)

    for n in fig_nums:
        svg = OUT_SVG / f"Figure{n}_composite.svg"
        if not svg.exists():
            print(f"  SKIP export Figure {n}: composite SVG not found")
            continue

        pdf = OUT_PDF / f"Figure{n}_composite.pdf"
        png = OUT_PNG / f"Figure{n}_composite.png"

        print(f"  Exporting Figure {n} → PDF")
        subprocess.run([
            INKSCAPE, str(svg),
            "--export-type=pdf",
            f"--export-filename={pdf}",
        ], capture_output=True, check=True)

        print(f"  Exporting Figure {n} → PNG ({DPI} DPI)")
        subprocess.run([
            INKSCAPE, str(svg),
            "--export-type=png",
            f"--export-filename={png}",
            f"--export-dpi={DPI}",
        ], capture_output=True, check=True)


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compose inline SVG figures")
    parser.add_argument("--fig", type=int, choices=[2, 3, 4, 5],
                        help="Compose a single figure (default: all)")
    parser.add_argument("--skip-convert", action="store_true",
                        help="Skip PDF→SVG conversion step")
    parser.add_argument("--skip-export", action="store_true",
                        help="Skip SVG→PDF/PNG export step")
    args = parser.parse_args()

    fig_nums = [args.fig] if args.fig else [2, 3, 4, 5]

    # Step 1: PDF → SVG
    if not args.skip_convert:
        print("\n── Step 1: Convert panel PDFs → SVGs ──")
        convert_pdfs_to_svgs(fig_nums)

    # Step 2: Compose
    print("\n── Step 2: Compose inline SVG composites ──")
    OUT_SVG.mkdir(parents=True, exist_ok=True)
    for n in fig_nums:
        info = FIGURES[n]
        print(f"  Composing Figure {n} ({info['prefix']}, {info['type']})...")
        if info["type"] == "rate":
            svg_root = compose_rate(n, info["prefix"], info["title"])
        else:
            svg_root = compose_conn(n, info["prefix"], info["title"])

        if svg_root is not None:
            out_path = OUT_SVG / f"Figure{n}_composite.svg"
            write_svg(svg_root, out_path)

    # Step 3: Export
    if not args.skip_export:
        print("\n── Step 3: Export composites → PDF + PNG ──")
        export_composites(fig_nums)

    print("\nDone.")


if __name__ == "__main__":
    main()
