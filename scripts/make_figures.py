#!/usr/bin/env python3
"""Regenerate the README figures in docs/figures/ (light and dark SVG pairs).

The numbers are the measured ones from the README's performance table (real CTCF
ChIP-seq, 5.0 M treatment + 5.0 M control reads, median of 3, release build) and
from `oracle/run_golden.sh`. Edit them here when either is re-measured.

Usage:  python3 scripts/make_figures.py
"""
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "docs" / "figures"

# (workload, macs3 seconds, macs3-rs seconds, macs3 MB, macs3-rs MB)
BENCH = [
    ("Single-end, bedGraph out (-B)", 42.7, 11.9, 304, 186),
    ("Single-end, --SPMR", 32.3, 8.1, 299, 173),
    ("Single-end, model mode", 44.4, 11.9, 296, 183),
    ("Single-end, no control", 22.9, 6.7, 169, 88),
    ("Paired-end BAM", 1.41, 0.41, 80, 31),
]

TILES = [
    ("7,204 / 7,204", "recorded runs byte-identical"),
    ("22,456 / 22,456", "output files byte-identical"),
    ("0", "exit-status mismatches"),
    ("3.4–4.0×", "faster on real ChIP-seq"),
]

THEMES = {
    "light": dict(surface="#fcfcfb", border="#e6e5e0", text="#0b0b0b", muted="#52514e",
                  grid="#e6e5e0", ref="#8a8983", bar="#2a78d6"),
    "dark": dict(surface="#1a1a19", border="#383835", text="#ffffff", muted="#c3c2b7",
                 grid="#383835", ref="#8a8983", bar="#3987e5"),
}

FONT = "-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif"


def bar(x, y, w, h, fill):
    """Horizontal bar: square at the baseline, 4px rounded at the data end."""
    r = min(4, w)
    return (f'<path d="M{x},{y}h{w - r:.1f}a{r},{r} 0 0 1 {r},{r}v{h - 2 * r}'
            f'a{r},{r} 0 0 1 -{r},{r}h-{w - r:.1f}z" fill="{fill}"/>')


def panel(t, x0, width, title, subtitle, values, vmax, ticks, fmt, refs):
    top, pitch, h = 96, 34, 16
    scale = width / vmax
    bottom = top + pitch * len(values) - (pitch - h) + 10
    out = [f'<text x="{x0}" y="40" font-size="15" font-weight="600" fill="{t["text"]}">{title}</text>',
           f'<text x="{x0}" y="60" font-size="12.5" fill="{t["muted"]}">{subtitle}</text>']
    for v in ticks:
        x = x0 + v * scale
        out.append(f'<line x1="{x:.1f}" y1="{top - 8}" x2="{x:.1f}" y2="{bottom}" stroke="{t["grid"]}"/>')
        out.append(f'<text x="{x:.1f}" y="{bottom + 16}" font-size="11.5" text-anchor="middle" '
                   f'fill="{t["muted"]}">{fmt(v)}</text>')
    for v, label in refs:
        x = x0 + v * scale
        out.append(f'<line x1="{x:.1f}" y1="{top - 8}" x2="{x:.1f}" y2="{bottom}" stroke="{t["ref"]}" stroke-width="1.5"/>')
        out.append(f'<text x="{x:.1f}" y="{top - 14}" font-size="11.5" text-anchor="middle" '
                   f'fill="{t["muted"]}">{label}</text>')
    for i, v in enumerate(values):
        y = top + i * pitch
        out.append(bar(x0, y, v * scale, h, t["bar"]))
        # Value at the tip, over a surface-coloured halo so a reference line never cuts it.
        label = f'x="{x0 + v * scale + 7:.1f}" y="{y + 12.5}" font-size="12.5" font-weight="600"'
        out.append(f'<text {label} fill="{t["surface"]}" stroke="{t["surface"]}" stroke-width="5" '
                   f'stroke-linejoin="round">{fmt(v)}</text>')
        out.append(f'<text {label} fill="{t["text"]}">{fmt(v)}</text>')
    return "\n".join(out)


def performance(t):
    W, H, label_x, a_x, b_x, pw = 880, 300, 24, 232, 572, 270
    rows = []
    for i, (name, *_rest) in enumerate(BENCH):
        rows.append(f'<text x="{label_x}" y="{96 + i * 34 + 12.5}" font-size="12.5" fill="{t["text"]}">{name}</text>')
    speed = [m / r for _, m, r, _, _ in BENCH]
    mem = [100 * r / m for _, _, _, m, r in BENCH]
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="{FONT}" role="img" aria-label="macs3-rs is 3.4 to 4.0 times faster than MACS3 and uses 39 to 62 percent of its peak memory across five callpeak workloads">
<rect x="0.5" y="0.5" width="{W - 1}" height="{H - 1}" rx="10" fill="{t["surface"]}" stroke="{t["border"]}"/>
{chr(10).join(rows)}
{panel(t, a_x, pw, "Speedup over MACS3", "wall clock, higher is better", speed, 5, [0, 2, 4], lambda v: f"{v:.1f}×" if v % 1 else f"{v:.0f}×", [(1, "MACS3"), (3, "target")])}
{panel(t, b_x, pw, "Peak memory vs MACS3", "share of upstream RSS, lower is better", mem, 125, [0, 25, 75], lambda v: f"{v:.0f}%", [(50, "target"), (100, "MACS3")])}
</svg>
'''


def tiles(t):
    W, H, pad, gap = 880, 112, 0, 16
    tw = (W - gap * (len(TILES) - 1)) / len(TILES)
    out = []
    for i, (value, label) in enumerate(TILES):
        x = pad + i * (tw + gap)
        out.append(f'<rect x="{x + 0.5:.1f}" y="0.5" width="{tw - 1:.1f}" height="{H - 1}" rx="10" '
                   f'fill="{t["surface"]}" stroke="{t["border"]}"/>')
        out.append(f'<text x="{x + 20:.1f}" y="54" font-size="23" font-weight="600" fill="{t["text"]}">{value}</text>')
        out.append(f'<text x="{x + 20:.1f}" y="82" font-size="13" fill="{t["muted"]}">{label}</text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
            f'font-family="{FONT}" role="img" aria-label="'
            + "; ".join(f"{v} {l}" for v, l in TILES) + '">\n' + "\n".join(out) + "\n</svg>\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, t in THEMES.items():
        (OUT / f"performance-{name}.svg").write_text(performance(t))
        (OUT / f"parity-{name}.svg").write_text(tiles(t))
    print(f"wrote 4 figures to {OUT}")


if __name__ == "__main__":
    main()
