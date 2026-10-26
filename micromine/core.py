"""MICRO-MINE core: FASTA parsing, SSR mining and locus clustering.

Scope: small organellar genomes (chloroplast of plants/algae, mitochondria of animals).
SSRs are grouped across species by the SAME MOTIF (canonical) and the SAME REGION, the
latter defined by the similarity of the flanking sequences (SNP tolerant and independent of
coordinates, rotation or strand of the genome).
"""
from __future__ import annotations

import bisect
import gzip
import math
import os
import re
from collections import Counter
from dataclasses import dataclass, field

CLASS_NAMES = {1: "mono", 2: "di", 3: "tri", 4: "tetra", 5: "penta", 6: "hexa"}
# Default minimum repeat numbers (same as MISA-Web)
DEFAULT_MIN_REPEATS = {1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5}

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_NAME_RE = re.compile(r"^(\d+)_(.+?)\.(?:fasta|fa|fna|fas)(?:\.gz)?$", re.I)


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def canonical_motif(motif: str) -> str:
    """Smallest rotation among the motif and its reverse complement (AT = TA, A = T...)."""
    cands = []
    for s in (motif, revcomp(motif)):
        cands.extend(s[i:] + s[:i] for i in range(len(s)))
    return min(cands)


def is_primitive(motif: str) -> bool:
    n = len(motif)
    return not any(n % d == 0 and motif == motif[:d] * (n // d) for d in range(1, n))


# --------------------------------------------------------------------------- #
# Input
# --------------------------------------------------------------------------- #
@dataclass
class Genome:
    idx: int
    tag: str          # numeric prefix as in the file name (e.g. "01")
    name: str         # species name (e.g. "speciesA")
    path: str
    contigs: list = field(default_factory=list)   # [(id, sequence)]

    @property
    def length(self) -> int:
        return sum(len(s) for _, s in self.contigs)

    @property
    def label(self) -> str:
        return f"{self.tag}_{self.name}"

    @property
    def n_count(self) -> int:
        return sum(seq.count("N") for _, seq in self.contigs)


def read_fasta(path: str):
    opener = gzip.open if path.lower().endswith(".gz") else open
    contigs, name, buf = [], None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    contigs.append((name, "".join(buf).upper()))
                name, buf = (line[1:].split() or ["seq"])[0], []
            else:
                if name is None:
                    raise ValueError(f"{path}: sequence found before the first '>' header")
                buf.append(line)
    if name is not None:
        contigs.append((name, "".join(buf).upper()))
    if not contigs:
        raise ValueError(f"{path}: no FASTA sequence found")
    return contigs


_VALID_BASES = frozenset("ACGTN")


def invalid_bases(contigs) -> Counter:
    """Count of every character that is not A, C, G, T or N (case-insensitive input is upper-cased)."""
    bad: Counter = Counter()
    for _, seq in contigs:
        if not set(seq) <= _VALID_BASES:
            bad.update(re.findall(r"[^ACGTN]", seq))
    return bad


def format_bad(bad: Counter, limit: int = 6) -> str:
    items = [f"{'space' if c == ' ' else c} x{n}" for c, n in bad.most_common(limit)]
    return ", ".join(items) + (", ..." if len(bad) > limit else "")


def _check_bases(genomes) -> None:
    problems = []
    for g in genomes:
        bad = invalid_bases(g.contigs)
        if bad:
            problems.append(f"{os.path.basename(g.path)} ({format_bad(bad)})")
    if problems:
        raise ValueError("Only the bases A, C, G, T and N are accepted. Files with other characters: "
                         + "; ".join(problems))


def discover_genomes(folder: str, max_files: int = 50) -> list[Genome]:
    gs = _discover(folder, max_files)
    _check_bases(gs)
    return gs


def _discover(folder: str, max_files: int = 50) -> list[Genome]:
    if os.path.isfile(folder):          # a single FASTA file: annotation mode
        base = os.path.basename(folder)
        m = _NAME_RE.match(base)
        tag, name = (m.group(1), m.group(2)) if m else (
            "01", re.sub(r"\.(fasta|fa|fna|fas)(\.gz)?$", "", base, flags=re.I))
        return [Genome(idx=0, tag=tag, name=name, path=folder, contigs=read_fasta(folder))]
    items = []
    for fn in os.listdir(folder):
        m = _NAME_RE.match(fn)
        if m:
            items.append((int(m.group(1)), m.group(1), m.group(2), os.path.join(folder, fn)))
    if not items:
        raise ValueError("No FASTA files named 'NN_species.fasta' (e.g. 01_speciesA.fasta) "
                         "found in: " + folder)
    items.sort()
    if len({i[0] for i in items}) != len(items):
        raise ValueError("Duplicated numeric prefixes in the file names.")
    if len(items) > max_files:
        raise ValueError(f"File limit exceeded: {len(items)} files (maximum {max_files}).")
    return [Genome(idx=pos, tag=tag, name=name, path=path, contigs=read_fasta(path))
            for pos, (_, tag, name, path) in enumerate(items)]


# --------------------------------------------------------------------------- #
# Mining
# --------------------------------------------------------------------------- #
@dataclass
class SSR:
    genome: int
    contig: str
    start: int        # 1-based, inclusive
    end: int
    motif: str        # as found
    canon: str
    reps: int
    seq: str
    left: str
    right: str
    offset: int = 0   # contig offset in the concatenated genome (used for plotting)
    cluster: str = ""
    orient: int = 1   # orientation relative to the cluster representative
    n_flank: int = 0  # N bases in the saved flanking regions (the SSR itself never contains N)

    @property
    def compromised(self) -> bool:
        """True when N bases occur in the SSR flanking region (primers/boundaries are unreliable)."""
        return self.n_flank > 0

    @property
    def k(self) -> int:
        return len(self.motif)

    @property
    def cls(self) -> str:
        return CLASS_NAMES[self.k]

    @property
    def mid(self) -> float:
        return self.offset + (self.start + self.end) / 2

    @property
    def flanked(self) -> str:
        """Left flank + SSR (lowercase) + right flank."""
        return self.left + self.seq.lower() + self.right


def case_flanked(s: SSR, ssr_case: str = "lower") -> str:
    """Flank + SSR + flank with the SSR highlighted by case.

    'lower': flanks UPPERCASE, SSR lowercase (default); 'upper': flanks lowercase, SSR UPPERCASE.
    """
    if ssr_case == "upper":
        return s.left.lower() + s.seq.upper() + s.right.lower()
    return s.left.upper() + s.seq.lower() + s.right.upper()


def _patterns(min_repeats: dict) -> dict:
    return {k: re.compile(r"([ACGT]{%d})\1{%d,}" % (k, m - 1))
            for k, m in min_repeats.items() if m and m >= 2}


def mine_genome(genome: Genome, min_repeats: dict, flank: int = 150) -> list[SSR]:
    pats = _patterns(min_repeats)
    out, offset = [], 0
    for cid, seq in genome.contigs:
        hits = []
        for k, pat in pats.items():
            for m in pat.finditer(seq):
                if k > 1 and not is_primitive(m.group(1)):
                    continue
                hits.append((m.start(), m.end(), m.group(1)))
        hits.sort(key=lambda h: (h[0], -(h[1] - h[0])))
        kept = []
        for h in hits:   # overlap between classes: keep the longest
            if kept and h[0] < kept[-1][1]:
                if h[1] - h[0] > kept[-1][1] - kept[-1][0]:
                    kept[-1] = h
                continue
            kept.append(h)
        for s, e, motif in kept:
            left, right = seq[max(0, s - flank):s], seq[e:e + flank]
            out.append(SSR(
                genome=genome.idx, contig=cid, start=s + 1, end=e, motif=motif,
                canon=canonical_motif(motif), reps=(e - s) // len(motif), seq=seq[s:e],
                left=left, right=right, offset=offset, n_flank=left.count("N") + right.count("N")))
        offset += len(seq)
    return out


# --------------------------------------------------------------------------- #
# Clustering by motif + region (flanking sequences)
# --------------------------------------------------------------------------- #
_K = 6


def _kset(s: str) -> set:
    return {s[i:i + _K] for i in range(len(s) - _K + 1)}


class _Anchor:
    __slots__ = ("f", "r")

    def __init__(self, left: str, right: str, alen: int):
        l, r = left[-alen:] if left else "", right[:alen]
        self.f = (_kset(l), _kset(r))
        self.r = (_kset(revcomp(r)), _kset(revcomp(l)))


def _sim(a: set, b: set):
    if not a or not b:
        return None
    return len(a & b) / min(len(a), len(b))


def _match(a: _Anchor, b: _Anchor, thr: float):
    best = None
    for orient, (bl, br) in ((1, b.f), (-1, b.r)):
        sims = [x for x in (_sim(a.f[0], bl), _sim(a.f[1], br)) if x is not None]
        if sims and all(x >= thr for x in sims):
            score = sum(sims) / len(sims)
            if best is None or score > best[1]:
                best = (orient, score)
    return best


@dataclass
class Cluster:
    id: str
    canon: str
    anchor: _Anchor
    members: list = field(default_factory=list)
    sub: dict = field(default_factory=dict)   # genome -> SSR below threshold at the homologous locus
    refpos: float = 0.0                       # position on the reference (first) genome axis

    @property
    def k(self) -> int:
        return len(self.canon)

    @property
    def cls(self) -> str:
        return CLASS_NAMES[self.k]

    @property
    def ref(self) -> SSR:
        return self.members[0]

    @property
    def number(self) -> int:
        return int(self.id.split(".")[1])

    def by_genome(self) -> dict:
        d: dict[int, list] = {}
        for m in self.members:
            d.setdefault(m.genome, []).append(m)
        return d

    @property
    def genomes(self) -> list:
        return sorted({m.genome for m in self.members})


def cluster_ssrs(ssrs: list[SSR], sim_thr: float = 0.6, anchor_len: int = 60) -> list[Cluster]:
    """Group SSRs with the same canonical motif and similar flanks (either strand)."""
    clusters: list[Cluster] = []
    by_motif: dict[str, list[Cluster]] = {}
    for s in ssrs:
        a = _Anchor(s.left, s.right, anchor_len)
        best = None
        for c in by_motif.get(s.canon, []):
            r = _match(c.anchor, a, sim_thr)
            if r and (best is None or r[1] > best[1][1]):
                best = (c, r)
        if best:
            s.orient = best[1][0]
            best[0].members.append(s)
        else:
            c = Cluster(id="", canon=s.canon, anchor=a, members=[s])
            clusters.append(c)
            by_motif.setdefault(s.canon, []).append(c)
    _assign_ids(clusters)
    return clusters


def _assign_ids(clusters: list[Cluster]) -> None:
    width = max(2, len(str(len(clusters))))
    for i, c in enumerate(clusters, 1):
        c.id = f"SSR.{i:0{width}d}"
        for m in c.members:
            m.cluster = c.id
        for h in c.sub.values():
            h.cluster = c.id


def exclusive_genome(c: Cluster, n_genomes: int):
    """Genome index if the SSR occurs in a single species (and there are >= 2 species)."""
    g = c.genomes
    return g[0] if n_genomes > 1 and len(g) == 1 else None


# --------------------------------------------------------------------------- #
# Homologous locus below threshold (separates real from threshold-only exclusives)
# --------------------------------------------------------------------------- #
def relaxed_mins(mins: dict, relax: float) -> dict:
    """Relaxed thresholds: ceil(min * relax), never below 3 repeats."""
    if not relax or relax <= 0:
        return {}
    return {k: max(3, math.ceil(m * relax)) for k, m in mins.items() if m and m >= 2}


def find_subthreshold(genomes, clusters, ssrs, mins, relax=0.6, anchor_len=60, sim_thr=0.6):
    """For every locus absent in a species, look for the same locus with fewer repeats."""
    rel = relaxed_mins(mins, relax)
    if not rel:
        return
    occupied: dict = {}
    for s in ssrs:
        occupied.setdefault((s.genome, s.contig), []).append((s.start, s.end))
    for g in genomes:
        need: dict = {}
        for c in clusters:
            if g.idx not in c.genomes:
                need.setdefault(c.canon, []).append(c)
        if not need:
            continue
        for h in mine_genome(g, rel, anchor_len):
            cands = need.get(h.canon)
            if not cands:
                continue
            if any(h.start <= e and h.end >= s0 for s0, e in occupied.get((g.idx, h.contig), [])):
                continue
            a = _Anchor(h.left, h.right, anchor_len)
            best = None
            for c in cands:
                r = _match(c.anchor, a, sim_thr)
                if r and (best is None or r[1] > best[1][1]):
                    best = (c, r)
            if best:
                c, (orient, _) = best
                h.orient, h.cluster = orient, c.id
                cur = c.sub.get(g.idx)
                if cur is None or h.reps > cur.reps:
                    c.sub[g.idx] = h


def exclusivity(c: Cluster, n_genomes: int):
    """None, 'real' (no homologous locus elsewhere) or 'threshold' (homolog below threshold)."""
    if exclusive_genome(c, n_genomes) is None:
        return None
    return "threshold" if c.sub else "real"


def is_polymorphic(c: Cluster) -> bool:
    return len({m.reps for m in c.members} | {h.reps for h in c.sub.values()}) > 1


def is_variable(c: Cluster, n_genomes: int) -> bool:
    """Marker candidate: absent (>= threshold) in some species or length polymorphic."""
    return len(c.genomes) < n_genomes or is_polymorphic(c)


# --------------------------------------------------------------------------- #
# Numbering: SSR.01, SSR.02... in top-to-bottom order of the reference (first) genome
# --------------------------------------------------------------------------- #
def _first_pos(c: Cluster, gi: int):
    ps = [m.mid for m in c.members if m.genome == gi]
    if gi in c.sub:
        ps.append(c.sub[gi].mid)
    return min(ps) if ps else None


def number_clusters(clusters: list[Cluster], ref: int = 0) -> list[Cluster]:
    """Order loci along the reference genome and renumber them (SSR.01 = top of the plot).

    Loci absent from the reference are placed by linear interpolation between the nearest
    loci shared with the reference in the first genome where they occur.
    """
    anchors: dict[int, tuple] = {}

    def anchors_for(gi):
        if gi not in anchors:
            pts = []
            for c in clusters:
                pr, pg = _first_pos(c, ref), _first_pos(c, gi)
                if pr is not None and pg is not None:
                    pts.append((pg, pr))
            pts.sort()
            anchors[gi] = (pts, [p[0] for p in pts])
        return anchors[gi]

    for c in clusters:
        p = _first_pos(c, ref)
        if p is not None:
            c.refpos = p
            continue
        gi = min(set(c.genomes) | set(c.sub))
        pg = _first_pos(c, gi)
        pts, keys = anchors_for(gi)
        i = bisect.bisect(keys, pg)
        lo = pts[i - 1] if i > 0 else None
        hi = pts[i] if i < len(pts) else None
        if lo and hi and abs(hi[1] - lo[1]) <= 3 * abs(hi[0] - lo[0]) + 5000:
            t = (pg - lo[0]) / (hi[0] - lo[0]) if hi[0] != lo[0] else 0.0
            c.refpos = lo[1] + t * (hi[1] - lo[1])
        elif lo or hi:
            near = lo if hi is None else hi if lo is None else (lo if pg - lo[0] <= hi[0] - pg else hi)
            c.refpos = near[1]
        else:
            c.refpos = pg
    ordered = sorted(clusters, key=lambda c: (c.refpos, c.canon))
    _assign_ids(ordered)
    return ordered


# --------------------------------------------------------------------------- #
# Variation of each locus across species (longer / shorter / absent relative to the others)
# --------------------------------------------------------------------------- #
VARIATION_TYPES = ("Conserved", "Length polymorphism", "Below threshold in some species",
                   "Absent in some species")


def locus_values(c: Cluster) -> dict:
    """genome -> (repeats, 'meets'|'below'); the longest copy is used when there are several."""
    vals = {gi: (max(m.reps for m in ms), "meets") for gi, ms in c.by_genome().items()}
    for gi, h in c.sub.items():
        vals[gi] = (h.reps, "below")
    return vals


def relative_status(c: Cluster, n_genomes: int) -> dict:
    """genome -> (relative, repeats, median of the other species, status).

    relative: 'Longer' / 'Shorter' / 'Same' (vs. the median of the OTHER species that carry the
    locus), 'Only species with SSR at locus' or 'No SSR at locus'.
    """
    vals = locus_values(c)
    out = {}
    for gi in range(n_genomes):
        if gi not in vals:
            out[gi] = ("No SSR at locus", None, None, None)
            continue
        r, st = vals[gi]
        others = sorted(v[0] for h, v in vals.items() if h != gi)
        if not others:
            out[gi] = ("Only species with SSR at locus", r, None, st)
            continue
        mid = len(others) // 2
        med = others[mid] if len(others) % 2 else (others[mid - 1] + others[mid]) / 2
        rel = "Longer" if r > med else "Shorter" if r < med else "Same"
        out[gi] = (rel, r, med, st)
    return out


def variation_type(c: Cluster, n_genomes: int) -> str:
    vals = locus_values(c)
    if len(vals) < n_genomes:
        return VARIATION_TYPES[3]
    if any(st == "below" for _, st in vals.values()):
        return VARIATION_TYPES[2]
    return VARIATION_TYPES[0] if len({r for r, _ in vals.values()}) == 1 else VARIATION_TYPES[1]


def variation_text(c: Cluster, genomes) -> str:
    """Compact description of the species that differ from the others."""
    parts = []
    for gi, (rel, r, med, st) in relative_status(c, len(genomes)).items():
        tag = genomes[gi].tag
        if rel == "Same":
            continue
        if rel == "No SSR at locus":
            parts.append(f"{tag}: no SSR at locus")
        elif rel == "Only species with SSR at locus":
            parts.append(f"{tag}: only species with SSR")
        else:
            parts.append(f"{tag}: {rel.lower()} ({r} vs {med:g}" + (", below threshold)" if st == "below" else ")"))
    return "; ".join(parts)


def compromised_genomes(c: Cluster) -> list:
    """Genome indexes in which at least one copy of the locus has N in its flanking region."""
    return sorted({m.genome for m in c.members if m.compromised})
