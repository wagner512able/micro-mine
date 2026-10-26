"""Command line interface and orchestration of MICRO-MINE."""
from __future__ import annotations

import argparse
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

from . import __author__, __contact__, __version__
from .core import (CLASS_NAMES, DEFAULT_MIN_REPEATS, cluster_ssrs, discover_genomes, exclusivity,
                   find_subthreshold, is_variable, mine_genome, number_clusters, relaxed_mins)
from .annotate import run_annotation
from .plot import plot_svg, svg_to_pdf
from .report import write_fastas, write_xlsx


def _parse_region(text):
    m = re.fullmatch(r"\s*([\d,_]+)\s*-\s*([\d,_]+)\s*", text)
    if not m:
        raise argparse.ArgumentTypeError("region must look like START-END, e.g. 20000-40000")
    a, b = (int(x.replace(",", "").replace("_", "")) for x in m.groups())
    if a >= b:
        raise argparse.ArgumentTypeError("region START must be smaller than END")
    return a, b


LAYOUTS = ("species", "locus", "consolidated")


def _parse_layout(text):
    parts = [t.strip().lower() for t in text.split(",") if t.strip()]
    if parts == ["none"]:
        return ()
    bad = [t for t in parts if t not in LAYOUTS]
    if bad:
        raise argparse.ArgumentTypeError(f"unknown FASTA layout {bad}; choose from {', '.join(LAYOUTS)} or none")
    return tuple(dict.fromkeys(parts))


def run(input_dir, output_dir, min_repeats=None, flank=150, anchor=60, similarity=0.6,
        max_files=50, threads=1, pdf=False, relax=0.6, regions=None, log=print, input_label=None,
        fasta_layout=("species", "locus", "consolidated"), ssr_case="lower"):
    mins = dict(DEFAULT_MIN_REPEATS)
    mins.update(min_repeats or {})
    genomes = discover_genomes(input_dir, max_files=max_files)
    if len(genomes) == 1:       # single FASTA: SSR annotation only
        if regions:
            log("[warning] zoom regions are ignored in single-genome annotation mode")
        return run_annotation(genomes[0], output_dir, mins, flank, pdf, ssr_case, fasta_layout, log,
                              input_label, input_dir)
    if anchor > flank:
        raise ValueError("the anchor length cannot be larger than the flank length")
    log(f"[MICRO-MINE] {len(genomes)} genome(s): " + ", ".join(g.label for g in genomes))

    if threads > 1 and len(genomes) > 1:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            results = list(ex.map(mine_genome, genomes, [mins] * len(genomes), [flank] * len(genomes)))
    else:
        results = [mine_genome(g, mins, flank) for g in genomes]
    ssrs = [s for res in results for s in res]   # order: file 01, 02, ... then position
    for g, res in zip(genomes, results):
        log(f"  {g.label}: {g.length:,} bp, {len(res)} SSRs")

    clusters = cluster_ssrs(ssrs, sim_thr=similarity, anchor_len=anchor)
    find_subthreshold(genomes, clusters, ssrs, mins, relax, anchor, similarity)
    clusters = number_clusters(clusters)          # SSR.01 = top of the reference genome
    kinds = [exclusivity(c, len(genomes)) for c in clusters]
    log(f"[MICRO-MINE] {len(clusters)} SSR loci; exclusive: {kinds.count('real')} real, "
        f"{kinds.count('threshold')} threshold-only")

    os.makedirs(output_dir, exist_ok=True)
    params = {"Software": "MICRO-MINE", "Version": __version__,
              "Developed by": f"{__author__} ({__contact__})", "Input folder": input_label or os.path.abspath(input_dir),
              "Reference genome (numbering/zoom)": genomes[0].label,
              "Flank length (bp)": flank, "Anchor length for locus matching (bp)": anchor,
              "Minimum flank similarity": similarity,
              "Relaxation factor (below-threshold search)": relax,
              "FASTA layout": ",".join(fasta_layout) or "none",
              "SSR highlight": "SSR lowercase, flanks UPPERCASE" if ssr_case != "upper" else "SSR UPPERCASE, flanks lowercase"}
    params.update({f"Min. repeats {name}": mins[k] for k, name in CLASS_NAMES.items()})
    params.update({f"Relaxed min. repeats {CLASS_NAMES[k]}": v for k, v in relaxed_mins(mins, relax).items()})
    write_xlsx(os.path.join(output_dir, "micromine_report.xlsx"), genomes, clusters, ssrs, params, ssr_case)
    write_fastas(output_dir, genomes, clusters, fasta_layout, ssr_case)

    plots = [(os.path.join(output_dir, "micromine_comparison.svg"), {}),
             (os.path.join(output_dir, "micromine_variable_markers.svg"),
              {"title": "MICRO-MINE - variable markers",
               "only": lambda c: is_variable(c, len(genomes))})]
    for a, b in regions or []:
        if not any(a <= c.refpos <= b for c in clusters):
            log(f"[warning] no SSR in region {a}-{b} of {genomes[0].label}; zoom plot skipped")
            continue
        plots.append((os.path.join(output_dir, f"micromine_region_{a}-{b}.svg"),
                      {"title": f"MICRO-MINE - region {a:,}-{b:,} bp", "region": (a, b)}))
    for path, kw in plots:
        plot_svg(path, genomes, clusters, **kw)
        if pdf and not svg_to_pdf(path, path.replace(".svg", ".pdf")):
            log("[warning] PDF output requires 'pip install cairosvg'; SVG files were created.")
            pdf = False
    log(f"[MICRO-MINE] results in: {input_label and 'job output folder' or os.path.abspath(output_dir)}")
    return genomes, clusters, ssrs


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="micromine",
        description="Mine and compare microsatellites (SSRs) across organellar genomes "
                    "(chloroplast of plants/algae; mitochondria of animals). "
                    "Input files: NN_species.fasta (the number sets the order of species).")
    p.add_argument("-i", "--input", required=True, help="folder with the FASTA files (01_sp.fasta, 02_sp.fasta ...) to compare species, or a single FASTA file (or a folder with one file) for SSR annotation only")
    p.add_argument("-o", "--output", default="micromine_out", help="output folder")
    for k, name in CLASS_NAMES.items():
        p.add_argument(f"--{name}", type=int, default=DEFAULT_MIN_REPEATS[k], metavar="N",
                       help=f"minimum repeats for {name}nucleotide SSRs (0 disables; default {DEFAULT_MIN_REPEATS[k]})")
    p.add_argument("--flank", type=int, default=150, help="flanking bp saved on each side of every SSR, e.g. 75 or 150 (default 150)")
    p.add_argument("--anchor", type=int, default=60, help="flank bp used to define the region (default 60)")
    p.add_argument("--similarity", type=float, default=0.6, help="minimum flank similarity, 0-1 (default 0.6)")
    p.add_argument("--relax", type=float, default=0.6,
                   help="factor (0-1) applied to the thresholds to look for the homologous locus below "
                        "threshold; 0 disables (default 0.6)")
    p.add_argument("--region", type=_parse_region, action="append", metavar="START-END",
                   help="also draw a zoomed plot of this region, in coordinates of the first (reference) "
                        "genome; can be repeated, e.g. --region 20000-40000")
    p.add_argument("--fasta-layout", type=_parse_layout, default=LAYOUTS, metavar="LIST",
                   help="FASTA outputs, comma separated: species (one file per species), locus (one file per "
                        "SSR locus with the copy from each species), consolidated (one file with all loci), "
                        "or none (default: species,locus,consolidated)")
    p.add_argument("--ssr-case", choices=("lower", "upper"), default="lower",
                   help="highlight the SSR by case: lower = SSR lowercase / flanks UPPERCASE (default); "
                        "upper = SSR UPPERCASE / flanks lowercase")
    p.add_argument("--max-files", type=int, default=50, help="maximum number of FASTA files (default 50)")
    p.add_argument("--threads", type=int, default=1, help="processes used for mining (default 1)")
    p.add_argument("--pdf", action="store_true", help="also write PDF plots (requires cairosvg)")
    p.add_argument("--version", action="version", version=__version__)
    a = p.parse_args(argv)
    mins = {k: getattr(a, name) for k, name in CLASS_NAMES.items()}
    try:
        run(a.input, a.output, mins, a.flank, a.anchor, a.similarity, a.max_files, a.threads,
            a.pdf, a.relax, a.region,
            fasta_layout=a.fasta_layout, ssr_case=a.ssr_case)
    except (ValueError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
