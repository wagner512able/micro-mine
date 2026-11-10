"""Synthetic benchmark: how well does MICRO-MINE recognise the SAME SSR locus in a diverged copy?

Simulates a plastid-derived sequence integrated in another organellar genome (an MTPT-like situation) or, more
generally, an orthologous locus that diverged between genomes. Genome 1 ("plastome") carries 140 SSR loci
(8 motifs, 1.5 kb of random flank on each side). Genome 2 ("mitogenome") carries one copy of every locus with random
substitutions at 0-20 % in the flanks, either keeping the SSR tract intact or mutating it as well.
For every SSR of genome 1 the script records whether MICRO-MINE reports the locus as
  shared     : the SSR meets the threshold in both genomes (same motif + similar flanks)
  below      : the homologous locus is found below the threshold in genome 2 (relaxed search)
  exclusive  : no homologous SSR found (the SSR would be called organelle-specific)

Usage: python benchmarks/synthetic_divergence.py [SEED]   (prints a table; deterministic for a given seed and version)
"""
import collections
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from micromine.cli import run  # noqa: E402

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 11
rng = random.Random(SEED)
MOTIFS = ["AT", "AG", "AAT", "ACC", "AAAT", "AATG", "AAGAG", "ATTCG"]
LEVELS = [0.0, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20]
N = 20


def rnd(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate, protect=None):
    out = list(s)
    for i, ch in enumerate(out):
        if protect and protect[0] <= i < protect[1]:
            continue
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != ch])
    return "".join(out)


def write(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


segments = []
for lvl in LEVELS:
    for _ in range(N):
        m = rng.choice(MOTIFS)
        reps = {2: 9, 3: 7, 4: 6, 5: 6}[len(m)]
        segments.append(rnd(1500) + m * reps + rnd(1500))
cp = "".join(segments)


def make_mt(tract_mutated):
    parts, idx = [rnd(8000)], 0
    for lvl in LEVELS:
        for _ in range(N):
            seg = segments[idx]
            tract = (1500, len(seg) - 1500)
            parts.append(mutate(seg, lvl, None if tract_mutated else tract))
            parts.append(rnd(6000))
            idx += 1
    return "".join(parts)


def classify(sim, tract_mutated):
    d = tempfile.mkdtemp()
    os.makedirs(d + "/in")
    write(d + "/in/01_plastome.fasta", "cp", cp)
    write(d + "/in/02_mitogenome.fasta", "mt", make_mt(tract_mutated))
    _, clusters, _ = run(d + "/in", d + "/out", similarity=sim, log=lambda *_: None, exclusive_fasta="none", fasta_layout=())
    starts = [sum(len(x) for x in segments[:i]) for i in range(len(segments))]
    res = collections.defaultdict(collections.Counter)
    for c in clusters:
        for m in c.members:
            if m.genome != 0:
                continue
            i = max(k for k, s0 in enumerate(starts) if s0 <= m.start - 1)
            lvl = LEVELS[i // N]
            res[lvl]["shared" if 1 in c.genomes else "below" if 1 in c.sub else "exclusive"] += 1
    return res


print("tract\tsimilarity\tdivergence_%\tshared\tbelow_threshold\texclusive\tdetected_%\tn")
for tm in (False, True):
    for sim in (0.6, 0.4):
        r = classify(sim, tm)
        for lvl in LEVELS:
            c = r[lvl]
            tot = sum(c.values())
            det = c["shared"] + c["below"]
            print(f"{'mutated' if tm else 'intact'}\t{sim}\t{int(lvl * 100)}\t{c['shared']}\t{c['below']}\t{c['exclusive']}\t{100 * det / tot:.0f}\t{tot}")
