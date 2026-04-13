#!/usr/bin/env python3
"""Compose publication-ready composite figures from individual panels.

Layout for Figs 2-3 (horizontal):
  Left column:  A (traces) over B (raster)
  Right column: C (paired firing) | D (categories firing) | E (violin firing)
                F (paired burst)  | G (categories burst)  | H (violin burst)

Layout for Figs 4-5:
  A (exemplar) over B (summary)
"""
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures_out"

DPI = 600
COL_WIDTH_PX = 4020
LABEL_SIZE = 80
PAD = 30
BG = (255, 255, 255)

try:
    FONT = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", LABEL_SIZE)
except Exception:
    FONT = ImageFont.load_default()


def load(name):
    p = FIG / name
    if not p.exists():
        print(f"  WARNING: {name} not found")
        return Image.new("RGB", (800, 600), BG)
    return Image.open(p).convert("RGB")


def add_label(img, label, x=15, y=5):
    draw = ImageDraw.Draw(img)
    draw.text((x, y), label, fill=(0, 0, 0), font=FONT)
    return img


def hstack(images, pad=PAD):
    h = max(im.height for im in images)
    w = sum(im.width for im in images) + pad * (len(images) - 1)
    out = Image.new("RGB", (w, h), BG)
    x = 0
    for im in images:
        out.paste(im, (x, (h - im.height) // 2))
        x += im.width + pad
    return out


def vstack(images, pad=PAD):
    w = max(im.width for im in images)
    h = sum(im.height for im in images) + pad * (len(images) - 1)
    out = Image.new("RGB", (w, h), BG)
    y = 0
    for im in images:
        out.paste(im, ((w - im.width) // 2, y))
        y += im.height + pad
    return out


def scale_to_width(img, target_w):
    ratio = target_w / img.width
    return img.resize((target_w, int(img.height * ratio)), Image.LANCZOS)


def scale_to_height(img, target_h):
    ratio = target_h / img.height
    return img.resize((int(img.width * ratio), target_h), Image.LANCZOS)


def compose_rate_figure(prefix, fig_num):
    """Compose a horizontal rate figure (Figs 2 or 3).

    Left column:  A (traces) over B (raster)
    Right column: C (paired firing) | D (categories firing) | E (violin firing)
                  F (paired burst)  | G (categories burst)  | H (violin burst)
    """
    print(f"Composing Figure {fig_num} ({prefix})...")

    traces = load(f"{prefix}_rate_panel_A_traces.png")
    raster = load(f"{prefix}_rate_panel_B_raster.png")
    scatter_spike = load(f"{prefix}_spike_rate_scatter.png")
    scatter_burst = load(f"{prefix}_burst_rate_scatter.png")
    cats_spike = load(f"{prefix}_spike_rate_change_categories.png")
    cats_burst = load(f"{prefix}_burst_rate_change_categories.png")
    violin_spike = load(f"{prefix}_spike_pct_change_violin.png")
    violin_burst = load(f"{prefix}_burst_pct_change_violin.png")

    # Target: right column takes 55% of total width, left takes 45%
    left_w = int(COL_WIDTH_PX * 0.42)
    right_w = COL_WIDTH_PX - left_w - PAD

    # Left column: traces (A) over raster (B), equal height
    traces_s = scale_to_width(traces, left_w)
    raster_s = scale_to_width(raster, left_w)
    left_col = vstack([add_label(traces_s, "A"), add_label(raster_s, "B")])

    # Right column: 2 rows × 3 panels
    # Each panel gets 1/3 of right_w
    panel_w = (right_w - 2 * PAD) // 3

    # Top row: firing rate panels
    c = scale_to_width(scatter_spike, panel_w)
    d = scale_to_width(cats_spike, panel_w)
    e = scale_to_width(violin_spike, panel_w)
    top_row = hstack([add_label(c, "C"), add_label(d, "D"), add_label(e, "E")])

    # Bottom row: burst rate panels
    f = scale_to_width(scatter_burst, panel_w)
    g = scale_to_width(cats_burst, panel_w)
    h = scale_to_width(violin_burst, panel_w)
    bot_row = hstack([add_label(f, "F"), add_label(g, "G"), add_label(h, "H")])

    # Stack right column rows
    right_col = vstack([top_row, bot_row])

    # Match heights: scale right to match left height
    target_h = left_col.height
    right_col = scale_to_height(right_col, target_h)

    # Combine left + right
    fig = hstack([left_col, right_col])
    fig = scale_to_width(fig, COL_WIDTH_PX)

    out = FIG / f"Figure{fig_num}_composite.png"
    fig.save(str(out), dpi=(DPI, DPI))
    print(f"  Saved {out} ({fig.width}x{fig.height})")


def compose_conn_figure(prefix, fig_num):
    """Compose a connectivity figure (Figs 4 or 5): A (exemplar) over B (summary)."""
    print(f"Composing Figure {fig_num} ({prefix})...")
    a = add_label(load(f"{prefix}_connectivity_exemplar.png"), "A")
    b = add_label(load(f"{prefix}_connectivity_summary.png"), "B")

    a_s = scale_to_width(a, COL_WIDTH_PX)
    b_s = scale_to_width(b, COL_WIDTH_PX)
    fig = vstack([a_s, b_s])
    fig = scale_to_width(fig, COL_WIDTH_PX)

    out = FIG / f"Figure{fig_num}_composite.png"
    fig.save(str(out), dpi=(DPI, DPI))
    print(f"  Saved {out} ({fig.width}x{fig.height})")


if __name__ == "__main__":
    compose_rate_figure("doi", 2)
    compose_rate_figure("ket", 3)
    compose_conn_figure("doi", 4)
    compose_conn_figure("ket", 5)
    print("All composites done.")
