"""Code shared by the two deployments of the web app (server with jobs / serverless).

Validation of uploads and parameters, and the summary shown on the results page.
"""
from __future__ import annotations

import json
import os
import re
from collections import Counter

from ..core import CLASS_NAMES, DEFAULT_MIN_REPEATS, VARIATION_TYPES, compromised_genomes, exclusivity, \
    format_bad, is_polymorphic, variation_type

MAX_FILES = 50
MIN_FILES = 1
MAX_CONTIGS = 2000
MAX_GENOME_BP = int(os.environ.get("MICROMINE_MAX_GENOME_BP", "2000000"))
REPO_URL = os.environ.get("MICROMINE_REPO_URL", "")
FASTA_EXT = re.compile(r"\.(fasta|fa|fna|fas)$", re.I)


def sanitize_label(label, i: int) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(label)).strip("._")[:60] or f"genome{i}"


def validate_uploads(items, max_file_bytes: int, max_genome_bp: int = MAX_GENOME_BP):
    """items: [(file name, raw bytes, label)] -> (cleaned [(NN_label.fasta, text)], problems [dict]).

    Every file is checked so that all the offending files can be reported together.
    """
    cleaned, problems = [], []
    for i, (name, raw, lab) in enumerate(items, 1):
        name = name or f"file{i}"
        try:
            if not FASTA_EXT.search(name):
                raise ValueError(f"{name}: expected a .fasta/.fa/.fna/.fas file")
            if len(raw) > max_file_bytes:
                raise ValueError(f"{name}: larger than {max_file_bytes // (1024 * 1024)} MB")
            text, _, bad = clean_fasta(raw, name, max_genome_bp)
        except ValueError as e:
            problems.append({"file": name, "message": str(e)})
            continue
        if bad:
            problems.append({"file": name, "bases": dict(bad.most_common(10)), "count": sum(bad.values()),
                             "message": f"{name}: contains bases other than A, C, G, T or N ({format_bad(bad)}); "
                                        f"only A, C, G, T and N are accepted"})
            continue
        cleaned.append((f"{i:02d}_{sanitize_label(lab, i)}.fasta", text))
    return cleaned, problems


def clean_fasta(raw: bytes, name: str, max_genome_bp: int = MAX_GENOME_BP):
    """Validate a FASTA and return (sanitized copy with safe contig names, length, invalid bases).

    Only A, C, G, T and N (any case) are valid; other characters are counted in the returned Counter.
    """
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        raise ValueError(f"{name}: not a plain-text FASTA file")
    out, total, contigs, bad = [], 0, 0, Counter()
    for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n"):
        ln = ln.strip()
        if not ln:
            continue
        if ln.startswith(">"):
            contigs += 1
            if contigs > MAX_CONTIGS:
                raise ValueError(f"{name}: more than {MAX_CONTIGS} sequences")
            tok = (ln[1:].split() or ["seq"])[0]
            out.append(">" + (re.sub(r"[^A-Za-z0-9_.|:-]", "_", tok)[:80] or "seq"))
        else:
            if not out:
                raise ValueError(f"{name}: sequence found before the first '>' header")
            up = ln.upper()
            if not set(up) <= set("ACGTN"):
                bad.update(re.findall(r"[^ACGTN]", up))
            total += len(ln)
            out.append(up)
    if contigs == 0 or total == 0:
        raise ValueError(f"{name}: no FASTA sequence found")
    if total > max_genome_bp:
        raise ValueError(f"{name}: {total:,} bp exceeds the {max_genome_bp:,} bp limit. This server "
                         f"accepts organellar genomes only (chloroplast of plants/algae, "
                         f"mitochondria of animals).")
    return "\n".join(out) + "\n", total, bad


def parse_params(raw) -> dict:
    """Validate the analysis parameters (a JSON string or an already-parsed dict)."""
    try:
        p = json.loads(raw) if isinstance(raw, (str, bytes)) else dict(raw)
    except (ValueError, TypeError):
        raise ValueError("Invalid parameters.")

    def num(key, lo, hi, cast, default):
        try:
            v = cast(p.get(key, default))
        except (TypeError, ValueError):
            raise ValueError(f"Parameter '{key}' must be a number.")
        if not lo <= v <= hi:
            raise ValueError(f"Parameter '{key}' must be between {lo} and {hi}.")
        return v

    mins = {}
    for k, name in CLASS_NAMES.items():
        mins[k] = num(name, 0, 50, int, DEFAULT_MIN_REPEATS[k]) if name in p else DEFAULT_MIN_REPEATS[k]
        if mins[k] == 1:
            raise ValueError(f"Minimum repeats for {name} must be 0 (disabled) or at least 2.")
    regions = []
    text = str(p.get("regions", "") or "")
    for part in re.split(r"[,;\n]+", text):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"([\d_]+)\s*-\s*([\d_]+)", part)
        if not m:
            raise ValueError(f"Region '{part}' must look like START-END, e.g. 20000-40000.")
        a, b = (int(x.replace("_", "")) for x in m.groups())
        if a >= b:
            raise ValueError(f"Region '{part}': START must be smaller than END.")
        regions.append([a, b])
    if len(regions) > 5:
        raise ValueError("At most 5 zoom regions are allowed.")
    layout = p.get("fasta_layout", ["species", "locus", "consolidated"])
    if not isinstance(layout, list) or any(x not in ("species", "locus", "consolidated") for x in layout):
        raise ValueError("Invalid FASTA layout.")
    ssr_case = p.get("ssr_case", "lower")
    if ssr_case not in ("lower", "upper"):
        raise ValueError("Invalid SSR highlight option.")
    flank, anchor = num("flank", 20, 500, int, 150), num("anchor", 20, 200, int, 60)
    if anchor > flank:
        raise ValueError("The anchor length cannot be larger than the flank length.")
    return {"min_repeats": mins, "fasta_layout": layout, "ssr_case": ssr_case, "flank": flank, "anchor": anchor,
            "similarity": num("similarity", 0.3, 1.0, float, 0.6),
            "relax": num("relax", 0.0, 1.0, float, 0.6), "regions": regions}


def build_summary(genomes, clusters, ssrs, out_dir) -> dict:
    """Summary of one analysis (used by the results page)."""
    files = sorted(f for f in os.listdir(out_dir) if os.path.isfile(os.path.join(out_dir, f)))

    def count(sub):
        d = os.path.join(out_dir, sub)
        return len(os.listdir(d)) if os.path.isdir(d) else 0

    if len(genomes) == 1:
        g = genomes[0]
        by_class = Counter(s.cls for s in ssrs)
        motifs = Counter((s.canon, s.cls) for s in ssrs).most_common(10)
        return {"mode": "annotation",
                "genome": {"tag": g.tag, "name": g.name, "length": g.length, "contigs": len(g.contigs),
                           "n_bases": g.n_count},
                "compromised": sum(1 for s in ssrs if s.compromised), "ssrs": len(ssrs),
                "density_per_kb": round(len(ssrs) / (g.length / 1000), 2) if g.length else 0,
                "by_class": {n: by_class.get(n, 0) for n in CLASS_NAMES.values()},
                "top_motifs": [[m, cls, c] for (m, cls), c in motifs], "files": files,
                "fasta_per_species": 0, "fasta_per_locus": 0}
    n = len(genomes)
    kinds = [exclusivity(c, n) for c in clusters]
    counts = {g.idx: 0 for g in genomes}
    for s in ssrs:
        counts[s.genome] += 1
    var = {t: 0 for t in VARIATION_TYPES}
    for c in clusters:
        var[variation_type(c, n)] += 1
    return {"mode": "comparison",
            "genomes": [{"tag": g.tag, "name": g.name, "length": g.length, "ssrs": counts[g.idx],
                         "n_bases": g.n_count,
                         "compromised": sum(1 for s in ssrs if s.genome == g.idx and s.compromised)}
                        for g in genomes],
            "compromised_loci": sum(1 for c in clusters if compromised_genomes(c)),
            "loci": len(clusters), "exclusive_real": kinds.count("real"),
            "exclusive_threshold": kinds.count("threshold"),
            "length_polymorphic": sum(1 for c in clusters if is_polymorphic(c)),
            "variation": var, "files": files,
            "fasta_per_species": count("fasta_per_species"), "fasta_per_locus": count("fasta_per_locus")}
