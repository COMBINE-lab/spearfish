#!/usr/bin/env python3
"""Generate the spearfish logo.

The fish is drawn from what spearfish does. Its spear-like bill is a long read,
studded with the seeds that anchor it. Its body carries the read's candidate
isoforms as colored lanes that run together where the isoforms share exons and
part in the forked tail where they diverge: the tail fork is the delta of delta
mapping, which aligns only where candidates differ.

Usage: generate-logo.py light|dark|favicon-light|favicon-dark > out.svg
"""

import sys

# Isoform colors: the COMBINE-lab palette shared with the cuttlefish logo.
ISO = ["#E4572E", "#F2A93B", "#17A398", "#4C6EF5"]  # coral, amber, teal, indigo
BODY = {"light": "#22303F", "dark": "#EAF0F7"}
EYE = {"light": "#FFFFFF", "dark": "#22303F"}

W, H = 512, 256


def body_path():
    # Fish facing left: snout at x=150, caudal peduncle at x=420.
    return (
        "M 150 128 "
        "C 172 104 222 90 286 92 "
        "C 344 94 396 114 422 123 "
        "L 422 133 "
        "C 396 142 344 162 286 164 "
        "C 222 166 172 152 150 128 Z"
    )


def dorsal_fin():
    # A long, rayed sail: highest just behind the head, sloping back.
    return (
        "M 190 101 "
        "C 196 80 204 62 216 56 "
        "C 252 62 310 80 364 100 "
        "C 320 95 258 94 190 101 Z"
    )


def fin_rays():
    """Thin rays across the sail."""
    out = []
    for x0, x1, y1 in [(204, 210, 66), (222, 230, 62), (242, 250, 66), (262, 270, 71),
                       (282, 290, 76), (302, 310, 82), (322, 330, 88)]:
        out.append((x0, 99, x1, y1))
    return out


def pectoral_fin():
    return "M 214 148 C 224 154 234 162 242 170 C 230 167 220 161 212 154 Z"


def lobe(sign):
    """One lobe of the forked tail; sign -1 = upper, +1 = lower."""
    def p(x, y):
        return f"{x} {128 + sign * (y - 128)}"
    return (
        f"M {p(414, 124)} "
        f"C {p(436, 108)} {p(458, 82)} {p(484, 48)} "
        f"C {p(482, 78)} {p(468, 106)} {p(446, 126)} Z"
    )


def tail():
    return lobe(-1) + " " + lobe(+1)


# Lane offsets (y) through the body, pinched through the peduncle, and each
# lane's path inside its tail lobe (outer lanes ride the lobes' leading edges).
BODY_Y = [119.0, 125.0, 131.0, 137.0]
PED_Y = [124.6, 126.9, 129.1, 131.4]
TAIL = [
    ((434, 118), (450, 96), (470, 66)),
    ((438, 122), (456, 102), (476, 76)),
    ((438, 134), (456, 154), (476, 180)),
    ((434, 138), (450, 160), (470, 190)),
]


def lanes():
    """Isoform lanes: shared through the body, splitting in the tail lobes."""
    out = []
    for i, c in enumerate(ISO):
        y = BODY_Y[i]
        (c1x, c1y), (c2x, c2y), (ex, ey) = TAIL[i]
        d = (
            f"M 196 128 C 214 {128 + (y - 128) * 0.35:.1f} 236 {y} 262 {y} "
            f"L 370 {y} "
            f"C 392 {y} 404 {PED_Y[i]} 420 {PED_Y[i]} "
            f"C {c1x} {c1y} {c2x} {c2y} {ex} {ey}"
        )
        out.append(
            f'<path d="{d}" fill="none" stroke="{c}" stroke-width="4.6" '
            f'stroke-linecap="round" stroke-linejoin="round"/>'
        )
    return out


def bill(theme):
    """The read: a long spear with seed anchors along it."""
    col = BODY[theme]
    out = [
        f'<path d="M 22 124.5 L 156 123 L 156 133 L 22 131.5 Z" fill="{col}"/>',
        f'<circle cx="22" cy="128" r="3.4" fill="{col}"/>',
    ]
    # Seeds: anchors shared by every isoform (two-tone would be too fine here),
    # drawn in the isoform colors in turn.
    for i, x in enumerate([44, 70, 96, 122]):
        out.append(f'<circle cx="{x}" cy="128" r="5.6" fill="{ISO[i]}"/>')
    return out


def logo(theme):
    col = BODY[theme]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="12 34 484 188" width="484" height="188" '
        'role="img" aria-label="A spearfish whose bill is a long read and whose body '
        'carries isoform lanes that split in its forked tail">',
        f'<path d="{dorsal_fin()}" fill="{col}" opacity="0.38"/>',
        *[f'<line x1="{a}" y1="{b}" x2="{c}" y2="{d}" stroke="{col}" stroke-width="1.6" opacity="0.6"/>'
          for a, b, c, d in fin_rays()],
        f'<path d="{tail()}" fill="{col}" opacity="0.38"/>',
        f'<path d="{body_path()}" fill="{col}"/>',
        f'<path d="{pectoral_fin()}" fill="{col}" opacity="0.55"/>',
        *lanes(),
        *bill(theme),
        f'<circle cx="176" cy="120" r="7.5" fill="{EYE[theme]}"/>',
        f'<circle cx="174" cy="120" r="3.6" fill="{col}"/>',
        "</svg>",
    ]
    return "\n".join(parts)


def favicon(theme):
    # A different drawing, not the logo scaled down: at 16 px the lanes, seeds
    # and fin rays turn to mush, and a wide fish is a sliver in a square tab.
    # Keep the spear, a solid body and two heavy lanes splitting into the tail
    # fork, and tilt the fish nose-up so it fills the square.
    col = BODY[theme]
    fish = [
        f'<path d="{tail()}" fill="{col}" opacity="0.45"/>',
        f'<path d="{body_path()}" fill="{col}"/>',
        f'<path d="M 230 121 L 392 121 C 420 121 440 96 472 64" fill="none" stroke="{ISO[0]}" '
        'stroke-width="16" stroke-linecap="round"/>',
        f'<path d="M 230 135 L 392 135 C 420 135 440 160 472 192" fill="none" stroke="{ISO[2]}" '
        'stroke-width="16" stroke-linecap="round"/>',
        f'<path d="M 24 121 L 156 118 L 156 138 L 24 135 Z" fill="{col}"/>',
    ]
    return "\n".join(
        [
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 256 256" width="256" height="256" '
            'role="img" aria-label="spearfish">',
            '<g transform="translate(128 128) rotate(38) scale(0.6) translate(-255 -128)">',
            *fish,
            "</g>",
            "</svg>",
        ]
    )


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "light"
    if which.startswith("favicon"):
        print(favicon(which.split("-")[1] if "-" in which else "light"))
    else:
        print(logo(which))
