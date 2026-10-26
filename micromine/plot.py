"""Synteny-style comparison plot (vector SVG; PDF optional via cairosvg).

Layout: one HORIZONTAL axis per species (file order, top to bottom), so the figure scales
to ~50 genomes by growing vertically. Position along the genome runs left to right.
"""
from __future__ import annotations

import math
from xml.sax.saxutils import escape

from .core import CLASS_NAMES, exclusivity

COLORS = {"mono": "#1f77b4", "di": "#d62728", "tri": "#2ca02c",
          "tetra": "#ff7f0e", "penta": "#9467bd", "hexa": "#8c564b"}
_STEPS = (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)


def _star(cx, cy, r):
    pts = []
    for i in range(10):
        ang = math.pi / 2 + i * math.pi / 5
        rad = r if i % 2 == 0 else r * 0.45
        pts.append(f"{cx + rad * math.cos(ang):.1f},{cy - rad * math.sin(ang):.1f}")
    return " ".join(pts)


def _cross(cx, cy, r, tip=""):
    """Red cross: SSR with N bases in its flanking region (compromised)."""
    t = f"<title>{escape(tip)}</title>" if tip else ""
    return (f'<g stroke="#c0392b" stroke-width="2.2" stroke-linecap="round">{t}'
            f'<line x1="{cx - r:.1f}" y1="{cy:.1f}" x2="{cx + r:.1f}" y2="{cy:.1f}"/>'
            f'<line x1="{cx:.1f}" y1="{cy - r:.1f}" x2="{cx:.1f}" y2="{cy + r:.1f}"/></g>')


def _nice_ge(v):
    for s in _STEPS:
        if s >= v:
            return s
    return _STEPS[-1]


def plot_svg(path, genomes, clusters, title="MICRO-MINE - SSR comparison", only=None, region=None):
    """Draw one horizontal axis per species and connect shared SSRs between consecutive rows.

    only   : optional predicate to draw a subset of loci
    region : (start, end) in reference (first genome) coordinates -> zoomed view
    """
    n = len(genomes)
    all_clusters = clusters
    drawn = [c for c in clusters
             if (only is None or only(c)) and (region is None or region[0] <= c.refpos <= region[1])]
    zoom = region is not None

    label_font = 11 if n <= 25 else 10
    longest = max(len(g.label) for g in genomes)
    ml = int(min(330, max(150, longest * label_font * 0.58 + 30)))
    plot_w = 900
    mr = 150 if zoom else 60
    top = 105
    row = max(46, min(130, 700 / max(n - 1, 1)))
    width = ml + plot_w + mr
    legend_h = 180
    height = int(top + (n - 1) * row + 50 + legend_h)
    ys = [top + i * row for i in range(n)]
    glen = {g.idx: g.length for g in genomes}

    # genomes mostly on the opposite strand relative to the reference are drawn flipped
    strand = {g.idx: 0 for g in genomes}
    for c in all_clusters:
        if c.ref.genome == 0:
            for m in c.members:
                if m.genome != 0:
                    strand[m.genome] += m.orient
    strand[0] = 0

    def opos(gi, mid):
        return glen[gi] - mid if strand[gi] < 0 else mid

    def copies(c, gi):
        cs = sorted(c.by_genome().get(gi, []), key=lambda m: m.mid)
        return cs[:1] if zoom else cs

    # axis windows (oriented coordinates)
    if zoom:
        win = {}
        for g in genomes:
            ps = []
            for c in drawn:
                ps += [opos(g.idx, m.mid) for m in copies(c, g.idx)]
                if g.idx in c.sub:
                    ps.append(opos(g.idx, c.sub[g.idx].mid))
                if g.idx == 0:
                    ps.append(c.refpos)
            win[g.idx] = (min(ps), max(ps)) if ps else (0, 1)
        span0 = max(h - l for l, h in win.values()) or 1
        pad = max(0.04 * span0, 150)
        win = {k: (l - pad, h + pad) for k, (l, h) in win.items()}
        span = max(h - l for l, h in win.values())
    else:
        win = {g.idx: (0, g.length) for g in genomes}
        span = max(g.length for g in genomes) or 1

    def xpos(gi, mid):
        return ml + (opos(gi, mid) - win[gi][0]) / span * plot_w

    ref_label = genomes[0].label
    sub1 = (f"Zoom: {region[0]:,}-{region[1]:,} bp of {ref_label} (axes scaled to the window)."
            if zoom else f"Axes proportional to genome size (longest = {span / 1000:.1f} kb).")
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height}" '
           f'viewBox="0 0 {width:.0f} {height}" font-family="Arial, Helvetica, sans-serif">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="20" y="32" font-size="18" font-weight="bold">{escape(title)}</text>',
           f'<text x="20" y="52" font-size="11" fill="#555">{escape(sub1)}</text>',
           '<text x="20" y="66" font-size="11" fill="#555">Lines connect SSRs with the same motif '
           'and region. The kb scale is shown on the first (reference) genome only.</text>']

    # connections between consecutive present species (vertical S-curves)
    out.append('<g fill="none" stroke-width="0.9" stroke-opacity="0.6">')
    for c in drawn:
        present = c.genomes
        if len(present) < 2:
            continue
        col = COLORS[c.cls]
        for a, b in zip(present, present[1:]):
            for x, y in zip(copies(c, a), copies(c, b)):
                x1, y1, x2, y2 = xpos(a, x.mid), ys[a], xpos(b, y.mid), ys[b]
                ym = (y1 + y2) / 2
                dash = ' stroke-dasharray="4,3"' if b != a + 1 else ""
                out.append(f'<path d="M{x1:.1f},{y1:.1f} C{x1:.1f},{ym:.1f} {x2:.1f},{ym:.1f} '
                           f'{x2:.1f},{y2:.1f}" stroke="{col}"{dash}/>')
    out.append("</g>")

    # axes, species labels and coordinate ticks (labels sit ABOVE the axis, away from SSR marks)
    base = 10000 if (not zoom and span >= 40000) else 0
    step = _nice_ge(max(base, span * 45 / plot_w))
    for g in genomes:
        lo, hi = win[g.idx]
        y = ys[g.idx]
        x0, x1 = ml, ml + (hi - lo) / span * plot_w
        out.append(f'<line x1="{x0:.1f}" y1="{y:.1f}" x2="{x1:.1f}" y2="{y:.1f}" '
                   f'stroke="#222" stroke-width="3" stroke-linecap="round"/>')
        out.append(f'<text x="{ml - 14}" y="{y + 4:.1f}" font-size="{label_font}" text-anchor="end">'
                   f'{escape(g.label)}</text>')
        a, b = (glen[g.idx] - hi, glen[g.idx] - lo) if strand[g.idx] < 0 else (lo, hi)
        r = math.ceil(max(a, 0) / step) * step
        while g.idx == 0 and r <= min(b, glen[g.idx]):   # scale only on the reference (first) genome
            x = xpos(g.idx, r)
            out.append(f'<line x1="{x:.1f}" y1="{y - 11:.1f}" x2="{x:.1f}" y2="{y - 7:.1f}" '
                       f'stroke="#aaa" stroke-width="0.7"/>'
                       f'<text x="{x:.1f}" y="{y - 14:.1f}" font-size="7" text-anchor="middle" '
                       f'fill="#a0a0a0">{"0" if r == 0 else f"{r / 1000:g} kb"}</text>')
            r += step
        if zoom:
            out.append(f'<text x="{x1 + 10:.1f}" y="{y + 3:.1f}" font-size="8" fill="#888">'
                       f'{max(a, 0):,.0f}-{min(b, glen[g.idx]):,.0f} bp</text>')

    # SSR marks on the axes (hover shows ID, motif and position)
    for c in drawn:
        for gi in c.genomes:
            for m in copies(c, gi):
                x = xpos(gi, m.mid)
                out.append(f'<line x1="{x:.1f}" y1="{ys[gi] - 5:.1f}" x2="{x:.1f}" y2="{ys[gi] + 5:.1f}" '
                           f'stroke="{COLORS[c.cls]}" stroke-width="1.6"><title>{c.id} ({c.canon}){m.reps}, '
                           f'{m.start:,}-{m.end:,}</title></line>')
        for gi, h in c.sub.items():
            x = xpos(gi, h.mid)
            out.append(f'<line x1="{x:.1f}" y1="{ys[gi] - 4:.1f}" x2="{x:.1f}" y2="{ys[gi] + 4:.1f}" '
                       f'stroke="#999" stroke-width="1.2"><title>{c.id}: ({c.canon}){h.reps} below '
                       f'threshold in {escape(genomes[gi].label)}</title></line>')

    # red cross = N bases in the SSR flanking region (compromised SSR), below the axis
    for c in drawn:
        for gi in c.genomes:
            for m in copies(c, gi):
                if m.compromised:
                    out.append(_cross(xpos(gi, m.mid), ys[gi] + 30, 5,
                                      f"{c.id}: {m.n_flank} N base(s) in the SSR flanking region (compromised)"))

    # gold star = real exclusive; black X = exclusive only because of the threshold (below the axis)
    for c in drawn:
        kind = exclusivity(c, n)
        if kind is None:
            continue
        ex = c.genomes[0]
        for m in copies(c, ex):
            cx, cy = xpos(ex, m.mid), ys[ex] + 15
            tip = (f'<title>{c.id} ({c.canon})n exclusive ' + ('(real)' if kind == "real" else
                   '(threshold-only)') + f' to {escape(genomes[ex].label)}</title>')
            if kind == "real":
                out.append(f'<polygon points="{_star(cx, cy, 7)}" fill="#FFD700" stroke="#B8860B" '
                           f'stroke-width="1">{tip}</polygon>')
            else:
                out.append(f'<g stroke="black" stroke-width="1.6" stroke-linecap="round">{tip}'
                           f'<line x1="{cx - 4:.1f}" y1="{cy - 4:.1f}" x2="{cx + 4:.1f}" y2="{cy + 4:.1f}"/>'
                           f'<line x1="{cx - 4:.1f}" y1="{cy + 4:.1f}" x2="{cx + 4:.1f}" y2="{cy - 4:.1f}"/></g>')

    # legend
    ly = ys[-1] + 60
    lx = 20
    for name in CLASS_NAMES.values():
        out.append(f'<line x1="{lx}" y1="{ly}" x2="{lx + 22}" y2="{ly}" stroke="{COLORS[name]}" '
                   f'stroke-width="3"/><text x="{lx + 28}" y="{ly + 4}" font-size="11">{name}</text>')
        lx += 85
    out.append(f'<polygon points="{_star(28, ly + 26, 7)}" fill="#FFD700" stroke="#B8860B" '
               f'stroke-width="1"/><text x="42" y="{ly + 30}" font-size="11">'
               f'exclusive marker (no homologous locus in the other species)</text>')
    out.append(f'<g stroke="black" stroke-width="1.6" stroke-linecap="round">'
               f'<line x1="24" y1="{ly + 44}" x2="32" y2="{ly + 52}"/>'
               f'<line x1="24" y1="{ly + 52}" x2="32" y2="{ly + 44}"/></g>'
               f'<text x="42" y="{ly + 52}" font-size="11">exclusive only because of the threshold '
               f'(locus present with fewer repeats in the others)</text>')
    out.append(f'<line x1="22" y1="{ly + 70}" x2="34" y2="{ly + 70}" stroke="#999" stroke-width="1.5"/>'
               f'<text x="42" y="{ly + 74}" font-size="11">SSR below threshold in that species</text>')
    out.append(_cross(28, ly + 92, 5) + f'<text x="42" y="{ly + 96}" font-size="11">SSR compromised: N bases in '
               f'its flanking region (unreliable boundaries/primers)</text>')
    out.append(f'<text x="20" y="{ly + 124}" font-size="10" fill="#666">Dashed line: SSR absent in the '
               f'intermediate species. Genomes on the reverse strand are drawn flipped.</text>')
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))


def plot_annotation_svg(path, genome, ssrs, title=None, max_rows=20):
    """Map of the SSRs along each contig (colored by class) + counts per class."""
    contigs = genome.contigs[:max_rows]
    longest = max(len(sq) for _, sq in contigs) or 1
    label_font = 11
    ml = int(min(330, max(110, max(len(c) for c, _ in contigs) * label_font * 0.58 + 30)))
    plot_w, top, row = 900, 125, 52
    width = ml + plot_w + 60
    bars_top = top + (len(contigs) - 1) * row + 70
    height = int(bars_top + 40 + 6 * 24 + 60)
    step = _nice_ge(max(10000 if longest >= 40000 else 0, longest * 45 / plot_w))
    by_class = {n: 0 for n in CLASS_NAMES.values()}
    for sr in ssrs:
        by_class[sr.cls] += 1
    dens = len(ssrs) / (genome.length / 1000) if genome.length else 0
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}" font-family="Arial, Helvetica, sans-serif">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="20" y="32" font-size="18" font-weight="bold">'
           f'{escape(title or "MICRO-MINE - SSR annotation")}</text>',
           f'<text x="20" y="52" font-size="11" fill="#555">{escape(genome.label)}: {genome.length:,} bp, '
           f'{len(ssrs):,} SSRs ({dens:.2f} SSRs per kb).</text>']
    if len(genome.contigs) > max_rows:
        out.append(f'<text x="20" y="84" font-size="11" fill="#b42318">Only the first {max_rows} of '
                   f'{len(genome.contigs)} contigs are drawn.</text>')
    out.append(_cross(26, 66, 5) + '<text x="38" y="70" font-size="11" fill="#555">SSR compromised: N bases in '
               'its flanking region</text>')
    pos = {}
    for i, (cid, sq) in enumerate(contigs):
        y = top + i * row
        pos[cid] = (y, len(sq))
        x1 = ml + len(sq) / longest * plot_w
        out.append(f'<line x1="{ml}" y1="{y}" x2="{x1:.1f}" y2="{y}" stroke="#222" stroke-width="3" '
                   f'stroke-linecap="round"/>')
        out.append(f'<text x="{ml - 14}" y="{y + 4}" font-size="{label_font}" text-anchor="end">{escape(cid)}</text>')
        r = 0
        while r <= len(sq):
            x = ml + r / longest * plot_w
            out.append(f'<line x1="{x:.1f}" y1="{y - 11}" x2="{x:.1f}" y2="{y - 7}" stroke="#aaa" stroke-width="0.7"/>'
                       f'<text x="{x:.1f}" y="{y - 14}" font-size="7" text-anchor="middle" fill="#a0a0a0">'
                       f'{"0" if r == 0 else f"{r / 1000:g} kb"}</text>')
            r += step
    for sr in ssrs:
        if sr.contig not in pos:
            continue
        y, _ = pos[sr.contig]
        x = ml + (sr.start + sr.end) / 2 / longest * plot_w
        out.append(f'<line x1="{x:.1f}" y1="{y - 7}" x2="{x:.1f}" y2="{y + 7}" stroke="{COLORS[sr.cls]}" '
                   f'stroke-width="1.6"><title>{sr.cluster} ({sr.canon}){sr.reps}, {sr.start:,}-{sr.end:,}</title></line>')
    for sr in ssrs:
        if sr.contig in pos and sr.compromised:
            y, _ = pos[sr.contig]
            out.append(_cross(ml + (sr.start + sr.end) / 2 / longest * plot_w, y + 17, 4.5,
                              f"{sr.cluster}: {sr.n_flank} N base(s) in the flanking region (compromised)"))
    out.append(f'<text x="20" y="{bars_top}" font-size="13" font-weight="bold">SSRs per class</text>')
    mx = max(by_class.values()) or 1
    for i, (name, n) in enumerate(by_class.items()):
        y = bars_top + 20 + i * 24
        w = n / mx * 520
        out.append(f'<text x="{ml - 14}" y="{y + 12}" font-size="11" text-anchor="end">{name}</text>'
                   f'<rect x="{ml}" y="{y}" width="{max(w, 1):.1f}" height="16" fill="{COLORS[name]}" rx="2"/>'
                   f'<text x="{ml + w + 8:.1f}" y="{y + 12}" font-size="11">{n:,}</text>')
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))


def svg_to_pdf(svg_path, pdf_path):
    try:
        import cairosvg
    except ImportError:
        return False
    cairosvg.svg2pdf(url=svg_path, write_to=pdf_path)
    return True
