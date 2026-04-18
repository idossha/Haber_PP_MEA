#!/usr/bin/env python3
"""Adjust font sizes and stroke widths in a composite SVG to Nature/NPP specs.

Walks the SVG XML tree, computes cumulative transforms per element, and
adjusts font-size / stroke-width so rendered sizes match Nature guidelines.

Nature/NPP specs:
  - Panel labels: 8pt bold, lowercase, Arial
  - Axis/tick labels: 5–7pt, Arial
  - Axis lines (spines, ticks): 0.5pt (0.18mm)
  - Data lines: >= 0.5pt

Usage:
    python3 scripts/fix_svg_nature_style.py fig2.svg -o fig2_nature.svg
    python3 scripts/fix_svg_nature_style.py fig3.svg --skip-stroke-groups g1202,g1203,g1204
"""
import argparse
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

NS = {"svg": "http://www.w3.org/2000/svg"}
PT_TO_MM = 25.4 / 72  # 0.3528 mm/pt

# ── Nature target sizes in pt ───────────────────────────────────────────
TARGETS = {
    "tick":       5,   # tick numbers, subscripts, superscripts
    "axis":       6,   # axis labels
    "title":      7,   # panel titles, group titles (Baseline/DOI)
    "panel":      8,   # panel labels (a, b, c ...)
    "annotation": 5,   # chi-square text, legend, p-values
    "axis_stroke": 0.5,
    "data_stroke": 0.5,
}


def parse_scale(transform_str):
    """Extract uniform scale from a transform attribute. Returns 1.0 if none."""
    if not transform_str:
        return 1.0
    m = re.match(r'matrix\(([0-9.e+-]+),\s*0,\s*0,\s*([0-9.e+-]+)', transform_str)
    if m:
        sx, sy = float(m.group(1)), float(m.group(2))
        return (abs(sx) + abs(sy)) / 2
    m = re.match(r'scale\(([0-9.e+-]+)', transform_str)
    if m:
        return float(m.group(1))
    return 1.0


def parse_font_size(style_str):
    """Extract font-size in px from a style string."""
    m = re.search(r'font-size:([0-9.]+)px', style_str)
    return float(m.group(1)) if m else None


def set_font_size(style_str, new_size):
    """Replace font-size value in a style string."""
    return re.sub(r'font-size:[0-9.]+px', f'font-size:{new_size:.4g}px', style_str)


def parse_stroke_width(style_str):
    """Extract stroke-width from a style string."""
    m = re.search(r'stroke-width:([0-9.]+)', style_str)
    return float(m.group(1)) if m else None


def set_stroke_width(style_str, new_width):
    """Replace stroke-width value in a style string."""
    return re.sub(r'stroke-width:([0-9.]+)', f'stroke-width:{new_width:.4g}', style_str)


def classify_font(current_mm):
    """Map a rendered font size in mm to a Nature target category."""
    current_pt = current_mm / PT_TO_MM
    if current_pt < 2.5:
        return "tick"       # very small → tick labels
    elif current_pt < 4.0:
        return "tick"       # small → tick labels
    elif current_pt < 6.0:
        return "axis"       # medium → axis labels
    elif current_pt < 9.0:
        return "title"      # larger → titles
    elif current_pt < 14.0:
        return "panel"      # panel labels
    else:
        return "title"      # very large → title (capped)


def classify_font_by_raw(raw_size, eff_scale):
    """Classify font based on raw SVG size and effective scale."""
    rendered_mm = raw_size * eff_scale
    rendered_pt = rendered_mm / PT_TO_MM

    # Very small rendered sizes (< 3pt) — these are in highly-scaled panels
    if rendered_pt < 3:
        # Use the raw size to infer role within the panel
        if raw_size < 5.0:
            return "tick"
        elif raw_size < 5.8:
            return "tick"
        elif raw_size < 6.5:
            return "tick"
        elif raw_size < 7.5:
            return "axis"
        elif raw_size < 10:
            return "title"
        else:
            return "title"

    # Larger rendered sizes — direct classification
    return classify_font(rendered_mm)


def target_size_for_role(role, eff_scale):
    """Compute the raw SVG font-size needed to render at the target pt size."""
    target_pt = TARGETS[role]
    target_mm = target_pt * PT_TO_MM
    return target_mm / eff_scale


def classify_stroke(raw_width, eff_scale):
    """Classify stroke: axis/tick lines vs data lines."""
    rendered_mm = raw_width * eff_scale
    rendered_pt = rendered_mm / PT_TO_MM
    if rendered_pt < 0.3:
        return "axis_stroke"  # thin → axis lines
    elif rendered_pt < 0.8:
        return "axis_stroke"
    else:
        return "data_stroke"  # thick → data


def target_stroke_for_role(role, eff_scale):
    """Compute raw SVG stroke-width for a Nature target."""
    target_pt = TARGETS[role]
    target_mm = target_pt * PT_TO_MM
    return target_mm / eff_scale


def fix_font_family(style_str):
    """Replace DejaVu Sans with Arial."""
    style_str = style_str.replace("font-family:'DejaVu Sans'", "font-family:Arial")
    style_str = style_str.replace(
        "-inkscape-font-specification:'DejaVu Sans'",
        "-inkscape-font-specification:'Arial'")
    style_str = style_str.replace(
        "-inkscape-font-specification:'DejaVu Sans Bold'",
        "-inkscape-font-specification:'Arial Bold'")
    return style_str


def walk_and_fix(elem, cumulative_scale, skip_stroke_ids, stats):
    """Recursively walk the SVG tree, fixing font sizes and stroke widths."""
    # Compute this element's contribution to scale
    transform = elem.get("transform", "")
    local_scale = parse_scale(transform)
    eff_scale = cumulative_scale * local_scale

    elem_id = elem.get("id", "")
    in_skip_group = elem_id in skip_stroke_ids

    # Fix style on this element
    style = elem.get("style", "")
    if style:
        modified = False

        # Font size
        fs = parse_font_size(style)
        if fs is not None and eff_scale > 0.001:
            role = classify_font_by_raw(fs, eff_scale)
            new_fs = target_size_for_role(role, eff_scale)
            if abs(new_fs - fs) > 0.01:
                style = set_font_size(style, new_fs)
                modified = True
                stats["fonts"] += 1

        # Stroke width (skip if in a protected group)
        if not in_skip_group:
            sw = parse_stroke_width(style)
            if sw is not None and sw > 0 and eff_scale > 0.001:
                role = classify_stroke(sw, eff_scale)
                new_sw = target_stroke_for_role(role, eff_scale)
                if abs(new_sw - sw) / max(sw, 0.001) > 0.05:
                    style = set_stroke_width(style, new_sw)
                    modified = True
                    stats["strokes"] += 1

        # Font family
        if "DejaVu" in style:
            style = fix_font_family(style)
            modified = True

        if modified:
            elem.set("style", style)

    # Propagate skip_stroke into children if this element is in the skip set
    child_skip = skip_stroke_ids
    if in_skip_group:
        child_skip = skip_stroke_ids | {"__all_children__"}
    if "__all_children__" in skip_stroke_ids:
        # We're inside a skipped group — mark this element as skipped too
        in_skip_group = True

    # Recurse into children
    for child in elem:
        child_eff_skip = child_skip
        if in_skip_group:
            child_eff_skip = child_skip | {"__all_children__"}
        walk_and_fix(child, eff_scale, child_eff_skip, stats)


def process_svg(input_path, output_path, skip_stroke_groups):
    """Process the SVG file using XML tree walking."""
    # Register namespaces to preserve them in output
    namespaces = {
        "": "http://www.w3.org/2000/svg",
        "svg": "http://www.w3.org/2000/svg",
        "xlink": "http://www.w3.org/1999/xlink",
        "inkscape": "http://www.inkscape.org/namespaces/inkscape",
        "sodipodi": "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd",
    }
    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)

    tree = ET.parse(input_path)
    root = tree.getroot()

    # Determine viewBox scale (viewBox units to mm)
    width_str = root.get("width", "190mm")
    vb_str = root.get("viewBox", "0 0 190 180")
    width_mm = float(re.match(r'([0-9.]+)', width_str).group(1))
    vb_w = float(vb_str.split()[2])
    base_scale = width_mm / vb_w  # viewBox units per mm

    print(f"Document: {width_str} wide, viewBox {vb_str}")
    print(f"Base scale: {base_scale:.4f} (viewBox unit = {base_scale:.4f} mm)")
    print(f"Skip stroke adjustment in groups: {skip_stroke_groups or 'none'}")
    print()

    skip_ids = set(skip_stroke_groups) if skip_stroke_groups else set()
    stats = {"fonts": 0, "strokes": 0}

    walk_and_fix(root, base_scale, skip_ids, stats)

    print(f"Modified {stats['fonts']} font-size attributes")
    print(f"Modified {stats['strokes']} stroke-width attributes")

    # Write output
    ET.indent(tree, space="  ")
    with open(output_path, "wb") as f:
        tree.write(f, encoding="utf-8", xml_declaration=True)

    print(f"Written: {output_path}")
    print(f"  Size: {output_path.stat().st_size / 1024:.0f} KB")


def main():
    parser = argparse.ArgumentParser(description="Fix SVG to Nature/NPP style")
    parser.add_argument("input", type=Path, help="Input SVG file")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output SVG (default: <input>_nature.svg)")
    parser.add_argument("--skip-stroke-groups", type=str, default=None,
                        help="Comma-separated group IDs to skip stroke adjustment "
                             "(e.g. g1202,g1203,g1204 for network panels)")
    args = parser.parse_args()

    if args.output is None:
        args.output = args.input.with_stem(args.input.stem + "_nature")

    skip_groups = args.skip_stroke_groups.split(",") if args.skip_stroke_groups else []
    process_svg(args.input, args.output, skip_groups)


if __name__ == "__main__":
    main()
