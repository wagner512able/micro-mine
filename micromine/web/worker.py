"""Job worker: runs one MICRO-MINE analysis inside its own process.

Usage: python -m micromine.web.worker JOB_DIR
JOB_DIR contains input/ (NN_species.fasta) and params.json; writes output/ and summary.json.
"""
from __future__ import annotations

import json
import os
import sys

from micromine.cli import run
from micromine.web.common import build_summary


def main(job_dir: str) -> int:
    with open(os.path.join(job_dir, "params.json")) as fh:
        p = json.load(fh)
    out = os.path.join(job_dir, "output")
    try:
        import cairosvg  # noqa: F401
        pdf = True
    except ImportError:
        pdf = False
    try:
        genomes, clusters, ssrs = run(
            os.path.join(job_dir, "input"), out,
            {int(k): int(v) for k, v in p["min_repeats"].items()},
            flank=p["flank"], anchor=p["anchor"], similarity=p["similarity"], max_files=50,
            threads=1, pdf=pdf, relax=p["relax"], regions=[tuple(r) for r in p["regions"]],
            input_label="uploaded files (web)",
            fasta_layout=tuple(p.get("fasta_layout", ())), ssr_case=p.get("ssr_case", "lower"))
    except (ValueError, OSError) as e:
        print(f"Error: {e}", flush=True)
        return 1
    with open(os.path.join(job_dir, "summary.json"), "w") as fh:
        json.dump(build_summary(genomes, clusters, ssrs, out), fh)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
