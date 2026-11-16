"""N handling: valid SSRs are retained, but N-containing flanks are never exported."""
import os
import sys
import tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from micromine.core import Genome, mine_genome
from micromine.cli import run

mins = {1: 99, 2: 6, 3: 99, 4: 99, 5: 99, 6: 99}

# A valid (AT)6 SSR with N in its left flanking window must remain annotated and flagged.
seq = "C" * 20 + "N" + "C" * 19 + "AT" * 6 + "G" * 40
g = Genome(idx=0, tag="01", name="test", path="memory", contigs=[("ctg1", seq)])
ssrs = mine_genome(g, mins, flank=30)
valid = [s for s in ssrs if s.motif.upper() == "AT" and s.reps >= 6]
assert len(valid) == 1, [(s.motif, s.reps, s.start, s.end) for s in ssrs]
assert valid[0].compromised and valid[0].n_flank == 1
assert "N" not in valid[0].seq.upper()

# An N interrupting a repeat tract must not be included inside any annotated SSR.
interrupted = "C" * 20 + "AT" * 3 + "N" + "AT" * 3 + "G" * 40
g2 = Genome(idx=0, tag="01", name="interrupted", path="memory", contigs=[("ctg1", interrupted)])
ssrs2 = mine_genome(g2, mins, flank=10)
assert not [s for s in ssrs2 if s.start <= 27 and s.end >= 20], ssrs2
assert all("N" not in s.seq.upper() for s in ssrs2)

# End-to-end annotation: sequence FASTA omits compromised SSRs and a TSV explains why.
with tempfile.TemporaryDirectory() as d:
    fasta = os.path.join(d, "01_test.fasta")
    with open(fasta, "w") as fh:
        fh.write(f">ctg1\n{seq}\n")
    out = os.path.join(d, "out")
    _, _, result = run(fasta, out, min_repeats=mins, flank=30, exclusive_fasta="none", log=lambda *_: None)
    assert any(s.compromised for s in result)
    fa = open(os.path.join(out, "SSRs_with_flanks.fasta")).read()
    assert "N" not in "".join(x for x in fa.splitlines() if x and not x.startswith(">"))
    warnings = open(os.path.join(out, "N_flank_warnings.tsv")).read()
    assert "SSR annotated; flanked sequence withheld" in warnings
    gff = open(os.path.join(out, "micromine_SSRs.gff3")).read()
    assert "compromised_by_N=1" in gff
print("OK: N handling validated.")
