#!/usr/bin/env python3
"""Replace raster colorbar images in MATLAB-exported SVGs with native SVG gradients.

MATLAB's Batik SVG generator renders colorbars as base64-encoded PNG <image>
elements. This script detects those raster colorbars, samples their colors,
and replaces them with native SVG <linearGradient> + <rect> elements — making
the colorbar resolution-independent and editable.

Usage:
    python3 scripts/fix_svg_colorbar.py figures/panels/standalone_adj/*.svg
    python3 scripts/fix_svg_colorbar.py --inplace figures/panels/standalone_adj/doi_adj_baseline.svg
    python3 scripts/fix_svg_colorbar.py --test   # run self-test with synthetic SVG

The default writes to <input>_fixed.svg; use --inplace to overwrite the original.
"""
import argparse
import base64
import io
import re
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

try:
    from PIL import Image
except ImportError:
    Image = None

# ── Namespace handling ───────────────────────────────────────────────────
NS = {
    "svg":    "http://www.w3.org/2000/svg",
    "xlink":  "http://www.w3.org/1999/xlink",
}
for prefix, uri in NS.items():
    ET.register_namespace(prefix, uri)
ET.register_namespace("", NS["svg"])


def _tag(local: str) -> str:
    return f"{{{NS['svg']}}}{local}"


# ── Clip-path geometry parsing ───────────────────────────────────────────

def _parse_clip_rect(clip_elem):
    """Extract bounding box (x, y, w, h) from a clipPath's child <path> or <rect>.

    Returns None if the clip path is not a simple rectangle.
    """
    for child in clip_elem:
        tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag

        if tag == "rect":
            x = float(child.get("x", 0))
            y = float(child.get("y", 0))
            w = float(child.get("width", 0))
            h = float(child.get("height", 0))
            return (x, y, w, h)

        if tag == "path":
            d = child.get("d", "")
            # Parse simple rectangular paths like "M0 370 L16 370 L16 0 L0 0 Z"
            coords = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)', d)
            if len(coords) >= 8:
                xs = [float(coords[i]) for i in range(0, len(coords), 2)]
                ys = [float(coords[i]) for i in range(1, len(coords), 2)]
                x_min, x_max = min(xs), max(xs)
                y_min, y_max = min(ys), max(ys)
                return (x_min, y_min, x_max - x_min, y_max - y_min)
    return None


def _is_tall_narrow(rect, min_aspect=5.0):
    """Check if a rectangle is tall and narrow (colorbar-shaped)."""
    if rect is None:
        return False
    _, _, w, h = rect
    if w == 0:
        return False
    return h / w >= min_aspect


# ── Color sampling from embedded PNG ────────────────────────────────────

def _decode_base64_png(data_uri: str):
    """Decode a base64 data-URI PNG and return a PIL Image."""
    if Image is None:
        raise ImportError("Pillow is required: pip install Pillow")
    # Strip the data URI prefix
    if "base64," in data_uri:
        b64 = data_uri.split("base64,", 1)[1]
    else:
        b64 = data_uri
    raw = base64.b64decode(b64)
    return Image.open(io.BytesIO(raw)).convert("RGBA")


def _sample_gradient_colors(img, n_stops=16):
    """Sample colors along the vertical center of the image.

    Returns list of (offset, hex_color) from bottom to top (offset 0 = bottom).
    """
    w, h = img.size
    cx = w // 2
    stops = []
    for i in range(n_stops):
        # Sample from bottom to top (SVG gradient y1=1 → y2=0)
        frac = i / (n_stops - 1)
        y = int((1.0 - frac) * (h - 1))
        y = max(0, min(h - 1, y))
        r, g, b, a = img.getpixel((cx, y))
        hex_col = f"#{r:02x}{g:02x}{b:02x}"
        stops.append((frac, hex_col))
    return stops


# ── SVG gradient construction ───────────────────────────────────────────

def _make_gradient_defs(gradient_id: str, stops):
    """Create a <linearGradient> element from color stops."""
    grad = ET.Element(_tag("linearGradient"), attrib={
        "id":             gradient_id,
        "x1":             "0",
        "y1":             "1",
        "x2":             "0",
        "y2":             "0",
        "gradientUnits":  "objectBoundingBox",
    })
    for offset, color in stops:
        ET.SubElement(grad, _tag("stop"), attrib={
            "offset":     f"{offset:.4f}",
            "stop-color": color,
        })
    return grad


# ── Main processing ────────────────────────────────────────────────────

def fix_svg_colorbar(svg_path: Path, inplace: bool = False, n_stops: int = 16) -> Path:
    """Process one SVG file: find raster colorbar, replace with native gradient."""
    tree = ET.parse(svg_path)
    root = tree.getroot()

    # 1. Index all clipPaths by ID
    clip_rects = {}
    for elem in root.iter():
        tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag
        if tag == "clipPath":
            cid = elem.get("id", "")
            rect = _parse_clip_rect(elem)
            if cid and rect:
                clip_rects[cid] = rect

    # 2. Find tall-narrow clipPaths (colorbar candidates)
    colorbar_clips = {cid: r for cid, r in clip_rects.items() if _is_tall_narrow(r)}
    if not colorbar_clips:
        print(f"  {svg_path.name}: no tall-narrow clip-path found — skipping")
        return svg_path

    # 3. Find <image> elements referencing colorbar clip-paths
    url_re = re.compile(r'url\(#([^)]+)\)')
    replaced = 0

    for elem in list(root.iter()):
        tag = elem.tag.split("}")[-1] if "}" in elem.tag else elem.tag

        # Check both <image> directly and <g> groups that contain <image>
        if tag == "image":
            parent_g = _find_parent(root, elem)
            clip_ref = _get_clip_id(elem) or (
                _get_clip_id(parent_g) if parent_g is not None else None
            )
            if clip_ref and clip_ref in colorbar_clips:
                replaced += _replace_image_with_gradient(
                    root, elem, parent_g, clip_ref, colorbar_clips[clip_ref], n_stops
                )
        elif tag == "g":
            clip_ref = _get_clip_id(elem)
            if clip_ref and clip_ref in colorbar_clips:
                # Find child <image> in this group
                for child in list(elem):
                    child_tag = child.tag.split("}")[-1] if "}" in child.tag else child.tag
                    if child_tag == "image":
                        replaced += _replace_image_with_gradient(
                            root, child, elem, clip_ref, colorbar_clips[clip_ref], n_stops
                        )

    if replaced == 0:
        print(f"  {svg_path.name}: no raster colorbar image matched clip-paths — skipping")
        return svg_path

    # 4. Write output
    if inplace:
        out_path = svg_path
    else:
        out_path = svg_path.with_stem(svg_path.stem + "_fixed")

    ET.indent(tree, space="  ")
    with open(out_path, "wb") as f:
        tree.write(f, encoding="utf-8", xml_declaration=True)

    print(f"  {svg_path.name}: replaced {replaced} colorbar(s) → {out_path.name}")
    return out_path


def _find_parent(root, target):
    """Find parent element of target in the tree."""
    for parent in root.iter():
        for child in parent:
            if child is target:
                return parent
    return None


def _get_clip_id(elem):
    """Extract clip-path ID from an element's style or attribute."""
    if elem is None:
        return None
    url_re = re.compile(r'url\(#([^)]+)\)')

    # Check clip-path attribute
    cp = elem.get("clip-path", "")
    m = url_re.search(cp)
    if m:
        return m.group(1)

    # Check style attribute
    style = elem.get("style", "")
    m = url_re.search(style)
    if m and "clip-path" in style:
        # Extract just the clip-path value
        for part in style.split(";"):
            if "clip-path" in part:
                m2 = url_re.search(part)
                if m2:
                    return m2.group(1)
    return None


def _replace_image_with_gradient(root, image_elem, parent_g, clip_id, clip_rect, n_stops):
    """Replace one <image> with a gradient <rect>. Returns 1 on success, 0 on failure."""
    # Get the image's data URI
    href = (image_elem.get("href") or
            image_elem.get(f"{{{NS['xlink']}}}href") or
            image_elem.get("{http://www.w3.org/1999/xlink}href", ""))
    if "base64," not in href:
        return 0

    try:
        img = _decode_base64_png(href)
    except Exception as e:
        print(f"    Warning: could not decode colorbar PNG: {e}")
        return 0

    stops = _sample_gradient_colors(img, n_stops)
    gradient_id = f"colorbar_gradient_{clip_id}"

    # Add gradient to <defs>
    defs = root.find(_tag("defs"))
    if defs is None:
        defs = ET.SubElement(root, _tag("defs"))
        # Move defs to be the first child
        root.remove(defs)
        root.insert(0, defs)
    # Also check nested defs
    for d in root.iter():
        dtag = d.tag.split("}")[-1] if "}" in d.tag else d.tag
        if dtag == "defs":
            defs = d
            break

    grad_elem = _make_gradient_defs(gradient_id, stops)
    defs.append(grad_elem)

    # Create replacement <rect> with same position/size
    cx, cy, cw, ch = clip_rect

    # Get the image's own x, y, width, height
    img_x = image_elem.get("x", "0")
    img_y = image_elem.get("y", "0")
    img_w = image_elem.get("width", str(cw))
    img_h = image_elem.get("height", str(ch))

    rect = ET.Element(_tag("rect"), attrib={
        "x":      img_x,
        "y":      img_y,
        "width":  img_w,
        "height": img_h,
        "fill":   f"url(#{gradient_id})",
    })

    # Copy clip-path from the image if it had one
    style = image_elem.get("style", "")
    if "clip-path" in style:
        rect.set("style", style)

    # Replace image with rect in parent
    if parent_g is not None:
        idx = list(parent_g).index(image_elem)
        parent_g.remove(image_elem)
        parent_g.insert(idx, rect)
    else:
        # Try root level
        idx = list(root).index(image_elem)
        root.remove(image_elem)
        root.insert(idx, rect)

    return 1


# ── Self-test ───────────────────────────────────────────────────────────

def _create_test_svg() -> str:
    """Create a minimal SVG mimicking MATLAB Batik output with a raster colorbar."""
    # Create a tiny 4x40 PNG gradient (blue→cyan→yellow) as base64
    if Image is None:
        raise ImportError("Pillow is required for testing")

    img = Image.new("RGBA", (4, 40))
    # Simple gradient: blue at bottom, yellow at top (parula-like)
    for y in range(40):
        t = y / 39.0  # 0 at top (yellow), 1 at bottom (blue)
        r = int(255 * (1 - t))
        g = int(255 * (1 - 0.3 * t))
        b = int(255 * t)
        for x in range(4):
            img.putpixel((x, y), (r, g, b, 255))

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    b64 = base64.b64encode(buf.getvalue()).decode("ascii")
    data_uri = f"data:image/png;base64,{b64}"

    # Also create a small heatmap PNG
    heatmap = Image.new("RGBA", (100, 100), (128, 200, 128, 255))
    buf2 = io.BytesIO()
    heatmap.save(buf2, format="PNG")
    b64_hm = base64.b64encode(buf2.getvalue()).decode("ascii")
    data_uri_hm = f"data:image/png;base64,{b64_hm}"

    svg = f"""<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink"
     width="252" height="227">
<!--Generated by the Batik Graphics2D SVG Generator-->
<defs id="genericDefs"/>
<g>
  <defs id="defs1">
    <clipPath clipPathUnits="userSpaceOnUse" id="clipPath1">
      <path d="M0 0 L252 0 L252 227 L0 227 L0 0 Z"/>
    </clipPath>
    <clipPath clipPathUnits="userSpaceOnUse" id="clipPath2">
      <path d="M0 0 L200 0 L200 200 L0 200 Z"/>
    </clipPath>
    <clipPath clipPathUnits="userSpaceOnUse" id="clipPath3">
      <path d="M0 370 L16 370 L16 0 L0 0 Z"/>
    </clipPath>
  </defs>
  <!-- Heatmap image -->
  <g transform="matrix(0.5,0,0,0.5,33,23.5)" style="clip-path:url(#clipPath2);">
    <image x="0" y="0" width="200" height="200"
           xlink:href="{data_uri_hm}" preserveAspectRatio="none"/>
  </g>
  <!-- Colorbar image (raster - this is what we want to replace) -->
  <g transform="translate(220,17)">
    <image x="0" y="0" width="16" height="370"
           xlink:href="{data_uri}"
           style="clip-path:url(#clipPath3);"
           preserveAspectRatio="none"/>
  </g>
  <!-- Colorbar tick marks -->
  <g style="stroke:rgb(38,38,38); stroke-width:0.6;">
    <line x1="236" y1="17" x2="240" y2="17"/>
    <line x1="236" y1="202" x2="240" y2="202"/>
  </g>
  <!-- Tick labels -->
  <text x="242" y="20" font-family="Arial" font-size="6">1.0</text>
  <text x="242" y="205" font-family="Arial" font-size="6">0.0</text>
</g>
</svg>"""
    return svg


def run_test():
    """Run self-test with a synthetic MATLAB-style SVG."""
    import tempfile

    print("=== fix_svg_colorbar self-test ===\n")

    svg_content = _create_test_svg()
    with tempfile.NamedTemporaryFile(suffix=".svg", delete=False, mode="w") as f:
        f.write(svg_content)
        test_path = Path(f.name)

    print(f"  Test input: {test_path}")
    print(f"  Input size: {test_path.stat().st_size} bytes")

    # Verify the input has an <image> with base64 data
    tree = ET.parse(test_path)
    images = list(tree.iter(_tag("image")))
    print(f"  Input <image> count: {len(images)}")

    # Run the fix
    out_path = fix_svg_colorbar(test_path, inplace=False, n_stops=10)

    # Verify output
    tree2 = ET.parse(out_path)
    images2 = list(tree2.iter(_tag("image")))
    grads = list(tree2.iter(_tag("linearGradient")))
    rects = [r for r in tree2.iter(_tag("rect")) if "url(#colorbar_gradient" in r.get("fill", "")]

    print(f"\n  Output: {out_path}")
    print(f"  Output <image> count: {len(images2)} (should be 1 — heatmap only)")
    print(f"  Output <linearGradient> count: {len(grads)} (should be 1)")
    print(f"  Output gradient <rect> count: {len(rects)} (should be 1)")

    # Check gradient stops
    if grads:
        stops = list(grads[0])
        print(f"  Gradient stops: {len(stops)}")
        for s in stops[:3]:
            print(f"    offset={s.get('offset')}  color={s.get('stop-color')}")
        if len(stops) > 3:
            print(f"    ... ({len(stops) - 3} more)")

    ok = len(images2) == 1 and len(grads) == 1 and len(rects) == 1
    print(f"\n  Result: {'PASS' if ok else 'FAIL'}")

    # Cleanup
    test_path.unlink()
    if out_path != test_path:
        # Keep the output for inspection
        print(f"  Output kept at: {out_path}")

    return ok


# ── CLI ─────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Replace raster colorbar images in MATLAB SVGs with native SVG gradients"
    )
    parser.add_argument("files", nargs="*", help="SVG files to process")
    parser.add_argument("--inplace", action="store_true",
                        help="Overwrite input files (default: write *_fixed.svg)")
    parser.add_argument("--stops", type=int, default=16,
                        help="Number of gradient color stops (default: 16)")
    parser.add_argument("--test", action="store_true",
                        help="Run self-test with synthetic SVG")
    args = parser.parse_args()

    if args.test:
        ok = run_test()
        sys.exit(0 if ok else 1)

    if not args.files:
        parser.error("No input files specified. Use --test for self-test.")

    for fpath in args.files:
        p = Path(fpath)
        if not p.exists():
            print(f"  WARNING: {p} not found — skipping")
            continue
        if p.suffix.lower() != ".svg":
            print(f"  WARNING: {p} is not an SVG — skipping")
            continue
        fix_svg_colorbar(p, inplace=args.inplace, n_stops=args.stops)


if __name__ == "__main__":
    main()
