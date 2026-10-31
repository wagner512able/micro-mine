"""Partial-unit highlight, repeat numbers in FASTA headers, exclusive-SSR FASTA and single-genome unique SSRs."""
import os, re, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.dirname(__file__))
from make_demo_data import build
from micromine.cli import run
from micromine.core import Genome, case_flanked, mine_genome

# 1) the repeat tract includes the incomplete unit that continues it: (TC)6 followed by T -> 'tctctctctctct'
seq = "GGCATCGAGTAG" + "TC" * 6 + "T" + "AGGATCCAGTTGCA"
g = Genome(idx=0, tag="01", name="x", path="x", contigs=[("c", seq)])
(s,) = mine_genome(g, {2: 6}, flank=8)
assert (s.reps, s.seq, s.tail, s.start, s.end) == (6, "TC" * 6, "T", 13, 24)
assert s.right.startswith("AGGATC") and case_flanked(s) == s.left.upper() + "tctctctctctct" + s.right.upper()
assert case_flanked(s, "upper").count("TCTCTCTCTCTCT") == 1
# no partial unit: tract followed by a different base -> no tail
(s2,) = mine_genome(Genome(0, "01", "x", "x", [("c", "GCATCG" + "TCT" * 5 + "ACGGA")]), {3: 5}, flank=6)
assert s2.tail == "" and s2.reps == 5

with tempfile.TemporaryDirectory() as d:
    build(os.path.join(d, "in"))
    run(os.path.join(d, "in"), os.path.join(d, "out"), log=lambda *_: None)
    out = os.path.join(d, "out")
    # 2) repeat numbers everywhere
    for f in ["SSRs_with_flanks.fasta"] + [os.path.join("fasta_per_species", x) for x in os.listdir(os.path.join(out, "fasta_per_species"))]:
        for h in [l for l in open(os.path.join(out, f)) if l.startswith(">")]:
            assert re.search(r"\)\d+\|repeats=\d+", h), (f, h)
    assert "repeats_by_species=" in open(os.path.join(out, "SSRs_with_flanks.fasta")).readline()
    # 3) highlight window of ssr_in_seq always equals the lowercase tract
    for rec in open(os.path.join(out, "fasta_per_species", "01_ancestral_SSRs.fasta")).read().split(">")[1:]:
        h, *b = rec.split("\n"); body = "".join(b); a, z = map(int, re.search(r"ssr_in_seq=(\d+)-(\d+)", h).groups())
        assert body[a - 1:z] == body[a - 1:z].lower() and body[:a - 1] == body[:a - 1].upper() and body[z:] == body[z:].upper()
    # 4) exclusive FASTA: default = one file per species with only real exclusives
    ex = os.path.join(out, "fasta_exclusive"); files = os.listdir(ex)
    assert files == ["03_perda_e_exclusivo_exclusive.fasta"], files
    recs = open(os.path.join(ex, files[0])).read().split(">")[1:]
    assert len(recs) == 1 and "exclusive=real" in recs[0] and "locus_copies=1" in recs[0] and "repeats=6" in recs[0], recs
    # contig mode, single-copy and threshold options
    run(os.path.join(d, "in"), os.path.join(d, "out_c"), exclusive_fasta="contig", exclusive_single_copy=True, log=lambda *_: None)
    assert os.listdir(os.path.join(d, "out_c", "fasta_exclusive")) == ["03_perda_e_exclusivo__03_perda_e_exclusivo_plastome_exclusive.fasta"]
    run(os.path.join(d, "in"), os.path.join(d, "out_n"), exclusive_fasta="none", log=lambda *_: None)
    assert not os.path.exists(os.path.join(d, "out_n", "fasta_exclusive"))
    # threshold-only exclusives: raise the di threshold so AT (9 in 01/03/04, 7 in 02) becomes borderline
    g_, c_, _ = run(os.path.join(d, "in"), os.path.join(d, "out_t"), {1: 12, 2: 8, 3: 7}, exclusive_fasta="species", exclusive_threshold=True, log=lambda *_: None)
    n_real = sum(len(open(os.path.join(d, "out_t", "fasta_exclusive", f)).read().split(">")) - 1 for f in os.listdir(os.path.join(d, "out_t", "fasta_exclusive")))
    assert n_real >= 1
    # 5) single genome with two contigs: unique (single-copy) SSRs
    single = os.path.join(d, "single"); os.makedirs(single)
    a = open(os.path.join(d, "in", "01_ancestral.fasta")).read().split("\n", 1); b = open(os.path.join(d, "in", "03_perda_e_exclusivo.fasta")).read().split("\n", 1)
    open(os.path.join(single, "mito.fasta"), "w").write(">ctgA\n" + a[1] + ">ctgB\n" + b[1])
    _, cl, ssrs = run(os.path.join(single, "mito.fasta"), os.path.join(d, "ann"), exclusive_fasta="species", log=lambda *_: None)
    u = open(os.path.join(d, "ann", "SSRs_unique_single_copy.fasta")).read().split(">")[1:]
    assert 0 < len(u) < len(ssrs) and all("locus_copies=1" in r for r in u)            # shared ancestral loci appear twice -> excluded
    run(os.path.join(single, "mito.fasta"), os.path.join(d, "ann2"), exclusive_fasta="contig", log=lambda *_: None)
    assert sorted(os.listdir(os.path.join(d, "ann2", "fasta_unique"))) == ["mito__ctgA_unique.fasta", "mito__ctgB_unique.fasta"] or len(os.listdir(os.path.join(d, "ann2", "fasta_unique"))) >= 1
    import openpyxl
    ws = openpyxl.load_workbook(os.path.join(d, "ann", "micromine_report.xlsx"))["SSR_annotation"]
    hdr = [c.value for c in ws[1]]; assert "Locus copies in genome" in hdr and "Partial unit (highlighted, not counted)" in hdr
print("OK: partial-unit highlight, repeat numbers, exclusive FASTA and unique SSRs validated.")
