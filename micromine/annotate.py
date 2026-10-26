"""Single-genome mode: SSR annotation of one FASTA file (no species comparison).

Outputs: Excel report, GFF3 annotation, FASTA of SSRs with flanks and a map of the SSRs along
the genome (SVG, PDF optional).
"""
from __future__ import annotations

import os
import urllib.parse
from collections import Counter

from openpyxl import Workbook

from . import __author__, __contact__, __version__
from .core import CLASS_NAMES, case_flanked, mine_genome
from .plot import plot_annotation_svg, svg_to_pdf
from .report import _finish, _seq_header, _wrap


def number_ssrs(ssrs):
    """SSR.01, SSR.02... in genome order (contigs in file order, then position)."""
    width = max(2, len(str(len(ssrs))))
    for i, s in enumerate(ssrs, 1):
        s.cluster = f"SSR.{i:0{width}d}"
    return ssrs


def write_annotation_xlsx(path, genome, ssrs, mins, params, ssr_case="lower"):
    from openpyxl.utils import get_column_letter  # noqa: F401  (kept local; used by _finish widths)
    wb = Workbook()
    ws = wb.active
    ws.title = "SSR_annotation"
    ws.append(["SSR_ID", "Contig", "Start", "End", "Length (bp)", "Class", "Motif (canonical)",
               "Observed motif", "Repeats", "Left flank (bp)", "Right flank (bp)", "N bases in flanks",
               "Compromised by N", _seq_header(ssr_case)])
    for s in ssrs:
        ws.append([s.cluster, s.contig, s.start, s.end, len(s.seq), s.cls, s.canon, s.motif, s.reps,
                   len(s.left), len(s.right), s.n_flank, "Yes" if s.compromised else "No",
                   case_flanked(s, ssr_case)])
    _finish(ws, {14: 60})

    wsum = wb.create_sheet("Summary")
    wsum.append(["Class", "Min. repeats", "SSRs", "% of SSRs", "SSRs per kb"])
    by_class = Counter(s.cls for s in ssrs)
    for k, name in CLASS_NAMES.items():
        r = wsum.max_row + 1
        wsum.append([name, mins.get(k, ""), by_class.get(name, 0),
                     f"=IF($C$8=0,0,C{r}/$C$8)", f"=IF($B$10=0,0,C{r}/($B$10/1000))"])
    wsum.append(["Total", "", "=SUM(C2:C7)", "=IF($C$8=0,0,C8/$C$8)", "=IF($B$10=0,0,C8/($B$10/1000))"])
    wsum.cell(10, 1, "Genome size (bp)")
    wsum.cell(10, 2, genome.length)
    wsum.cell(11, 1, "Contigs")
    wsum.cell(11, 2, len(genome.contigs))
    wsum.cell(12, 1, "N bases in genome")
    wsum.cell(12, 2, genome.n_count)
    wsum.cell(13, 1, "SSRs compromised by N")
    wsum.cell(13, 2, sum(1 for s in ssrs if s.compromised))
    _finish(wsum, {1: 18})
    for r in range(2, 9):
        wsum.cell(r, 4).number_format = "0.0%"
        wsum.cell(r, 5).number_format = "0.00"

    wm = wb.create_sheet("Motif_counts")
    wm.append(["Motif (canonical)", "Class", "SSRs", "Min. repeats observed", "Max. repeats observed",
               "Mean repeats"])
    motifs: dict[str, list] = {}
    for s in ssrs:
        motifs.setdefault(s.canon, []).append(s)
    for motif, ms in sorted(motifs.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        reps = [m.reps for m in ms]
        wm.append([motif, ms[0].cls, len(ms), min(reps), max(reps), round(sum(reps) / len(reps), 2)])
    _finish(wm)

    wq = wb.create_sheet("Parameters")
    wq.append(["Parameter", "Value"])
    for k, v in params.items():
        wq.append([k, str(v)])
    _finish(wq, {1: 46, 2: 40})
    wb.calculation.fullCalcOnLoad = True
    wb.save(path)


def write_gff3(path, genome, ssrs):
    q = lambda t: urllib.parse.quote(str(t), safe="|:._-")
    with open(path, "w") as fh:
        fh.write("##gff-version 3\n")
        for cid, seq in genome.contigs:
            fh.write(f"##sequence-region {q(cid)} 1 {len(seq)}\n")
        for s in ssrs:
            attrs = (f"ID={s.cluster};Name=({s.motif}){s.reps};motif={s.motif};canonical_motif={s.canon};"
                     f"repeats={s.reps};class={s.cls};length={len(s.seq)}"
                     + (f";compromised_by_N={s.n_flank}" if s.n_flank else ""))
            fh.write("\t".join([q(s.contig), "MicroMine", "microsatellite", str(s.start), str(s.end),
                                ".", "+", ".", attrs]) + "\n")


def write_annotation_fasta(path, genome, ssrs, ssr_case="lower"):
    with open(path, "w") as fh:
        for s in ssrs:
            fh.write(f">{s.cluster}|{genome.label}|{s.contig}:{s.start}-{s.end}|({s.motif}){s.reps}"
                     f"|flank={len(s.left)}/{len(s.right)}|ssr_in_seq={len(s.left) + 1}-{len(s.left) + len(s.seq)}"
                     + (f"|WARNING_N={s.n_flank}" if s.n_flank else "") + "\n"
                     f"{_wrap(case_flanked(s, ssr_case))}\n")


def run_annotation(genome, output_dir, mins, flank=150, pdf=False, ssr_case="lower",
                   fasta_layout=("consolidated",), log=print, input_label=None, input_path=None):
    log(f"[MICRO-MINE] single genome: annotation mode (no comparison) for {genome.label}")
    ssrs = number_ssrs(mine_genome(genome, mins, flank))
    by_class = Counter(s.cls for s in ssrs)
    log(f"  {genome.label}: {genome.length:,} bp, {len(ssrs)} SSRs ("
        + ", ".join(f"{n} {by_class.get(n, 0)}" for n in CLASS_NAMES.values()) + ")")
    os.makedirs(output_dir, exist_ok=True)
    params = {"Software": "MICRO-MINE", "Version": __version__,
              "Developed by": f"{__author__} ({__contact__})", "Mode": "single-genome annotation",
              "Input": input_label or os.path.abspath(input_path or genome.path),
              "Flank length (bp)": flank,
              "SSR highlight": "SSR lowercase, flanks UPPERCASE" if ssr_case != "upper" else "SSR UPPERCASE, flanks lowercase"}
    params.update({f"Min. repeats {name}": mins[k] for k, name in CLASS_NAMES.items()})
    write_annotation_xlsx(os.path.join(output_dir, "micromine_report.xlsx"), genome, ssrs, mins, params, ssr_case)
    write_gff3(os.path.join(output_dir, "micromine_SSRs.gff3"), genome, ssrs)
    if fasta_layout:
        write_annotation_fasta(os.path.join(output_dir, "SSRs_with_flanks.fasta"), genome, ssrs, ssr_case)
    svg = os.path.join(output_dir, "micromine_annotation_map.svg")
    plot_annotation_svg(svg, genome, ssrs)
    if pdf and not svg_to_pdf(svg, svg.replace(".svg", ".pdf")):
        log("[warning] PDF output requires 'pip install cairosvg'; the SVG file was created.")
    log(f"[MICRO-MINE] results in: {input_label and 'job output folder' or os.path.abspath(output_dir)}")
    return [genome], [], ssrs
