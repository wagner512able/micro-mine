"""Generate 4 synthetic plastomes to test MICRO-MINE.

01 ancestral | 02 flank SNPs + repeat-number changes
03 SSR #6 broken by a SNP + one exclusive SSR | 04 rotated and reverse strand
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from micromine.core import revcomp  # noqa: E402

DEFS = [("A", 12), ("T", 13), ("AT", 9), ("AG", 7), ("AAT", 7), ("ACC", 6),
        ("AAAT", 6), ("AATG", 5), ("AAGAG", 5), ("ATTCG", 5), ("AAGTCT", 5), ("C", 11)]


def rnd(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, pos):
    return s[:pos] + {"A": "C", "C": "G", "G": "T", "T": "A"}[s[pos]] + s[pos + 1:]


def build(outdir, seed=7):
    rng = random.Random(seed)
    flanks = [(rnd(rng, 600), rnd(rng, 600)) for _ in DEFS]
    uniq = (rnd(rng, 600), "CTTAGA" * 6, rnd(rng, 600))

    def make(counts, snp_flank=(), broken=(), extra=False):
        parts, bounds = [], []
        for i, ((m, _), (l, r)) in enumerate(zip(DEFS, flanks)):
            if i in snp_flank:
                l, r = mutate(l, 570), mutate(r, 25)
                l = mutate(l, 540)
            ssr = m * counts[i]
            if i in broken:
                mid = len(ssr) // 2
                ssr = ssr[:mid] + mutate(ssr, mid)[mid] + ssr[mid + 1:]
            bounds.append(sum(len(p) for p in parts))
            parts.append(l + ssr + r)
        if extra:
            parts.append("".join(uniq))
        return "".join(parts), bounds

    base = [c for _, c in DEFS]
    g1, bounds = make(base)
    c2 = list(base); c2[2] = 7; c2[4] = 9
    g2, _ = make(c2, snp_flank=(0, 3, 6, 9))
    g3, _ = make(base, broken=(5,), extra=True)
    cut = bounds[6]
    g4 = revcomp(g1[cut:] + g1[:cut])

    os.makedirs(outdir, exist_ok=True)
    for name, g in [("01_ancestral", g1), ("02_snps_flanco", g2), ("03_perda_e_exclusivo", g3),
                    ("04_rotacionado_rc", g4)]:
        with open(os.path.join(outdir, name + ".fasta"), "w") as fh:
            fh.write(f">{name}_plastome\n")
            for i in range(0, len(g), 70):
                fh.write(g[i:i + 70] + "\n")


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "demo_data")
