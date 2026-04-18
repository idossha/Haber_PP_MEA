#!/usr/bin/env python3
"""Compose publication-ready composite figures with inline SVG panel content.

Each panel's vector content is embedded as a <g> group (Inkscape layer),
NOT as base64 <image> data URIs. Every path, text element, and shape from
each panel remains individually selectable and editable in Inkscape.

Pipeline steps (all run by this script):
  1. Convert panel PDFs → SVGs via Inkscape  (output/fig{N}/panels/ → output/_intermediate/panel_svgs/)
  2. Inline panel SVGs into composite SVGs    (→ output/fig{N}/composite/)
  3. Export composites → PDF + PNG via Inkscape (→ output/fig{N}/composite/)

Usage:
    python3 scripts/compose_svg_inline.py                  # all figures
    python3 scripts/compose_svg_inline.py --fig 2          # Figure 2 only
    python3 scripts/compose_svg_inline.py --skip-convert    # skip PDF→SVG step
    python3 scripts/compose_svg_inline.py --skip-export     # skip SVG→PDF/PNG step
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET
from PIL import Image

Image.MAX_IMAGE_PIXELS = None  # 600-DPI panels are large but legitimate

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
INTERMEDIATE = OUTPUT / "_intermediate"
PANEL_SVGS = INTERMEDIATE / "panel_svgs"
PANEL_PNGS = INTERMEDIATE / "panel_pngs"

# Figure mapping: figure number -> output subdirectory
_FIG_DIR = {2: "fig2", 3: "fig3", 4: "fig4", 5: "fig5"}

def _panels_dir(fig_num: int) -> Path:
    return OUTPUT / _FIG_DIR[fig_num] / "panels"

def _composite_dir(fig_num: int) -> Path:
    return OUTPUT / _FIG_DIR[fig_num] / "composite"

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
    3: {"type": "conn", "prefix": "doi", "title": "Figure 3 — DOI functional connectivity"},
    4: {"type": "rate", "prefix": "ket", "title": "Figure 4 — Ketanserin firing and burst rates"},
    5: {"type": "conn", "prefix": "ket", "title": "Figure 5 — Ketanserin functional connectivity"},
}


# ── Step 1: PDF → SVG conversion ────────────────────────────────────────

def convert_pdfs_to_svgs(fig_nums):
    """Convert panel PDFs to SVGs using Inkscape."""
    import shutil

    # Collect all needed (fig_num, basename) pairs
    needed = []
    for n in fig_nums:
        info = FIGURES[n]
        panels = rate_panels(info["prefix"]) if info["type"] == "rate" else conn_panels(info["prefix"])
        for basename in panels.values():
            needed.append((n, basename))

    PANEL_SVGS.mkdir(parents=True, exist_ok=True)
    PANEL_PNGS.mkdir(parents=True, exist_ok=True)

    converted = 0
    for fig_num, basename in sorted(set(needed)):
        panels_dir = _panels_dir(fig_num)
        pdf = panels_dir / f"{basename}.pdf"
        svg = PANEL_SVGS / f"{basename}.svg"
        png = PANEL_PNGS / f"{basename}.png"

        if not pdf.exists():
            print(f"  WARNING: {pdf.name} not found in {panels_dir} — skipping")
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
        src_png = panels_dir / f"{basename}.png"
        if src_png.exists() and (not png.exists() or src_png.stat().st_mtime > png.stat().st_mtime):
            shutil.copy2(src_png, png)

    print(f"  PDF→SVG: {converted} converted")


# ── SVG parsing and composition helpers ──────────────────────────────────

# Panels that should be embedded as raster PNGs (too heavy as vector)
RASTER_PANELS = {"traces", "raster"}


def get_png_dimensions(png_path: Path):
    """Get pixel dimensions of a PNG without loading full image."""
    with Image.open(png_path) as img:
        return img.width, img.height


def make_raster_panel(png_path: Path, panel_id: str, x: float, y: float,
                      target_w: float, target_h: float = None,
                      composite_dir: Path = None):
    """Link a PNG as an external <image> inside an Inkscape layer group.

    Uses a relative path so the SVG stays portable within the repo.
    The PNG must remain at its path relative to the composite SVG.
    """
    px_w, px_h = get_png_dimensions(png_path)

    sx = target_w / px_w
    if target_h is not None:
        sy = target_h / px_h
        s = min(sx, sy)
    else:
        s = sx
    actual_h = px_h * s

    # Relative path from the composite directory to the PNG
    ref_dir = composite_dir if composite_dir else png_path.parent
    rel = os.path.relpath(png_path, ref_dir)

    g = ET.Element("g", attrib={
        "id": panel_id,
        "{http://www.inkscape.org/namespaces/inkscape}label": panel_id,
        "{http://www.inkscape.org/namespaces/inkscape}groupmode": "layer",
    })
    ET.SubElement(g, "image", attrib={
        "x": f"{x:.2f}",
        "y": f"{y:.2f}",
        "width": f"{target_w:.2f}",
        "height": f"{actual_h:.2f}",
        "href": rel,
        "preserveAspectRatio": "xMidYMid meet",
    })

    return g, [], actual_h


def _prefix_ids(elem, prefix):
    """Recursively prefix all 'id' attributes and update references (url(#...), href=#...)."""
    import re
    url_re = re.compile(r'url\(#([^)]+)\)')
    href_re = re.compile(r'^#(.+)$')

    # Collect all original IDs first
    old_ids = set()
    for e in elem.iter():
        eid = e.get("id")
        if eid:
            old_ids.add(eid)

    if not old_ids:
        return

    # Build mapping
    id_map = {oid: f"{prefix}_{oid}" for oid in old_ids}

    # Apply
    for e in elem.iter():
        # Rename id
        eid = e.get("id")
        if eid and eid in id_map:
            e.set("id", id_map[eid])

        # Update all attributes that may reference IDs
        for attr, val in list(e.attrib.items()):
            # url(#id) references (clip-path, fill, stroke, filter, mask, marker)
            new_val = url_re.sub(lambda m: f"url(#{id_map.get(m.group(1), m.group(1))})", val)
            # href="#id" references (xlink:href, href)
            href_m = href_re.match(new_val)
            if href_m and href_m.group(1) in id_map:
                new_val = f"#{id_map[href_m.group(1)]}"
            if new_val != val:
                e.set(attr, new_val)

        # Also check style attribute for url(#...) refs
        style = e.get("style")
        if style:
            new_style = url_re.sub(lambda m: f"url(#{id_map.get(m.group(1), m.group(1))})", style)
            if new_style != style:
                e.set("style", new_style)


def parse_panel(svg_path: Path, id_prefix: str = ""):
    """Parse a panel SVG, return (viewBox_w, viewBox_h, children, defs).

    If id_prefix is given, all IDs and their references are prefixed to
    avoid collisions when multiple panels are inlined in one document.
    """
    tree = ET.parse(svg_path)
    root = tree.getroot()

    vb = root.get("viewBox")
    if vb:
        parts = vb.split()
        vb_w, vb_h = float(parts[2]), float(parts[3])
    else:
        vb_w = float(root.get("width", "100"))
        vb_h = float(root.get("height", "100"))

    # Prefix all IDs in the entire tree before splitting
    if id_prefix:
        _prefix_ids(root, id_prefix)

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
    vb_w, vb_h, children, defs = parse_panel(svg_path, id_prefix=panel_id)

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
    # Namespace prefixes are handled by ET.register_namespace() calls above.
    # Only set non-namespace attributes here to avoid duplicate xmlns decls.
    root = ET.Element(f"{{{NS['svg']}}}svg", attrib={
        "width": f"{PAGE_W}mm",
        "height": f"{page_h:.1f}mm",
        "viewBox": f"0 0 {PAGE_W} {page_h:.1f}",
        "version": "1.1",
    })

    title = ET.SubElement(root, "title")
    title.text = title_text

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

def embed_panel(panel_key, basename, panel_id, x, y, target_w,
                target_h=None, panels_dir=None, composite_dir=None):
    """Route to raster or vector embedding based on panel type."""
    if panel_key in RASTER_PANELS:
        png = PANEL_PNGS / f"{basename}.png"
        if not png.exists() and panels_dir is not None:
            png = panels_dir / f"{basename}.png"
        if not png.exists():
            print(f"    WARNING: {basename}.png not found for raster embed")
            return None, [], 0
        print(f"    {panel_id}: raster (PNG)")
        return make_raster_panel(png, panel_id, x, y, target_w, target_h,
                                 composite_dir=composite_dir)
    else:
        svg = PANEL_SVGS / f"{basename}.svg"
        if not svg.exists():
            print(f"    WARNING: {basename}.svg not found for vector embed")
            return None, [], 0
        return make_panel_group(svg, panel_id, x, y, target_w, target_h)


def compose_rate(fig_num, prefix, title_text):
    """Compose a rate figure: left col A/B, right col C/D/E over F/G/H."""
    panels = rate_panels(prefix)
    panels_dir = _panels_dir(fig_num)

    # Check availability: raster panels need PNGs, others need SVGs
    missing = []
    for key, basename in panels.items():
        if key in RASTER_PANELS:
            if not (PANEL_PNGS / f"{basename}.png").exists() and \
               not (panels_dir / f"{basename}.png").exists():
                missing.append(key)
        else:
            if not (PANEL_SVGS / f"{basename}.svg").exists():
                missing.append(key)
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

    composite_dir = _composite_dir(fig_num)

    # Left column: A (traces, raster PNG) over B (raster, raster PNG)
    y = TOP_MARGIN
    gA, dA, hA = embed_panel("traces", panels["traces"], "panel_a_traces", 0, y, left_w,
                             panels_dir=panels_dir, composite_dir=composite_dir)
    groups.append(gA); all_defs.extend(dA)
    labels.append(make_label("a", lbl_ox, y + lbl_oy + LABEL_SIZE_MM))

    y_B = y + hA + GAP_MM
    gB, dB, hB = embed_panel("raster", panels["raster"], "panel_b_raster", 0, y_B, left_w,
                             panels_dir=panels_dir, composite_dir=composite_dir)
    groups.append(gB); all_defs.extend(dB)
    labels.append(make_label("b", lbl_ox, y_B + lbl_oy + LABEL_SIZE_MM))

    left_h = y + hA + GAP_MM + hB

    # Right column: 2 rows × 3 panels (vector SVG), matched to left height
    row_h = (left_h - TOP_MARGIN - GAP_MM) / 2

    for i, (key, pid, lbl) in enumerate([
        ("scatter_spike", "panel_c_scatter_spike", "c"),
        ("cats_spike",    "panel_d_cats_spike",    "d"),
        ("violin_spike",  "panel_e_violin_spike",  "e"),
    ]):
        px = right_x + i * (rpanel_w + PANEL_GAP)
        g, d, h = embed_panel(key, panels[key], pid, px, TOP_MARGIN, rpanel_w, row_h)
        groups.append(g); all_defs.extend(d)
        labels.append(make_label(lbl, px + lbl_ox, TOP_MARGIN + lbl_oy + LABEL_SIZE_MM))

    y_bot = TOP_MARGIN + row_h + GAP_MM
    for i, (key, pid, lbl) in enumerate([
        ("scatter_burst", "panel_f_scatter_burst", "f"),
        ("cats_burst",    "panel_g_cats_burst",    "g"),
        ("violin_burst",  "panel_h_violin_burst",  "h"),
    ]):
        px = right_x + i * (rpanel_w + PANEL_GAP)
        g, d, h = embed_panel(key, panels[key], pid, px, y_bot, rpanel_w, row_h)
        groups.append(g); all_defs.extend(d)
        labels.append(make_label(lbl, px + lbl_ox, y_bot + lbl_oy + LABEL_SIZE_MM))

    page_h = max(left_h, y_bot + row_h) + BOTTOM_MARGIN
    return assemble_svg(groups, labels, all_defs, title_text, page_h)


# ── Step 2b: Compose connectivity figures (Figs 4, 5) ───────────────────

def compose_conn(fig_num, prefix, title_text):
    """Compose a connectivity figure: A (exemplar) over B (summary), full width."""
    panels = conn_panels(prefix)
    svg_files = {k: PANEL_SVGS / f"{v}.svg" for k, v in panels.items()}

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
    for n in fig_nums:
        comp_dir = _composite_dir(n)
        comp_dir.mkdir(parents=True, exist_ok=True)

        svg = comp_dir / f"Figure{n}_composite.svg"
        if not svg.exists():
            print(f"  SKIP export Figure {n}: composite SVG not found")
            continue

        pdf = comp_dir / f"Figure{n}_composite.pdf"
        png = comp_dir / f"Figure{n}_composite.png"

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
    for n in fig_nums:
        info = FIGURES[n]
        comp_dir = _composite_dir(n)
        comp_dir.mkdir(parents=True, exist_ok=True)

        print(f"  Composing Figure {n} ({info['prefix']}, {info['type']})...")
        if info["type"] == "rate":
            svg_root = compose_rate(n, info["prefix"], info["title"])
        else:
            svg_root = compose_conn(n, info["prefix"], info["title"])

        if svg_root is not None:
            out_path = comp_dir / f"Figure{n}_composite.svg"
            write_svg(svg_root, out_path)

    # Step 3: Export
    if not args.skip_export:
        print("\n── Step 3: Export composites → PDF + PNG ──")
        export_composites(fig_nums)

    print("\nDone.")


if __name__ == "__main__":
    main()
