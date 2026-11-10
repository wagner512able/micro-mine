"""Scalability benchmark of the REAL MICRO-MINE pipeline (reading, mining, locus matching, reports, figures).

Self-contained: one file, no other script needed. Every measurement runs `micromine.cli.run` in a separate process
(this same file, started with --child), so that wall time and peak memory belong to MICRO-MINE alone.

A) number of genomes: N copies (N = 2..50) of a real chloroplast genome, each with 0.5 % random substitutions
B) genome size: synthetic AT-rich genomes (0.5-5 Mb) with ~0.14 SSRs/kb; annotation of one genome and
   comparison of two genomes

Usage (from the repository root):
    python benchmarks/scaling.py PLASTOME.fasta [REPEATS] [--quick] [--save results.tsv]

    PLASTOME.fasta  a real chloroplast genome (single record), e.g. one of your Syagrus plastomes
    REPEATS         repetitions of each measurement (default 3); the table shows the MEDIAN time and the MAX memory
    --quick         smaller sweep (2/5/10 genomes; 0.5/1 Mb) to check that everything works
    --save FILE     also write the tables to FILE

The first output line records the machine (OS, CPU, logical CPUs, Python version): report it together with the results.
Run times depend on the load of the machine, so always use REPEATS >= 3 and close other programs.
Peak memory is the peak resident/working set of the whole process (interpreter and libraries included), as shown by
the task manager. On Windows it is read with the Win32 API (PeakWorkingSetSize); if it cannot be read, 'nan' is printed.
"""
import os
import platform
import random
import shutil
import statistics
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))


# --------------------------------------------------------------------------- child: one real analysis
def peak_memory_mb():
    try:
        import resource                                   # Linux / macOS
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return peak / (1024 * 1024) if sys.platform == "darwin" else peak / 1024
    except ImportError:                                   # Windows
        try:
            import ctypes
            from ctypes import wintypes

            class PMC(ctypes.Structure):
                _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                            ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                            ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

            kernel32, psapi = ctypes.windll.kernel32, ctypes.windll.psapi
            kernel32.GetCurrentProcess.restype = wintypes.HANDLE
            psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(PMC), wintypes.DWORD]
            pmc = PMC()
            pmc.cb = ctypes.sizeof(PMC)
            psapi.GetProcessMemoryInfo(kernel32.GetCurrentProcess(), ctypes.byref(pmc), pmc.cb)
            return pmc.PeakWorkingSetSize / (1024 * 1024)
        except Exception:                                 # memory could not be read: do not invent a number
            return float("nan")


def child(folder, out):
    from micromine.cli import run
    t0 = time.perf_counter()
    _, clusters, ssrs = run(folder, out, log=lambda *_: None, fasta_layout=(), exclusive_fasta="none")
    print(f"{time.perf_counter() - t0:.2f}\t{peak_memory_mb():.0f}\t{len(ssrs)}\t{len(clusters)}")


# --------------------------------------------------------------------------- parent: build data and measure
def read_seq(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">")).upper()


def write(path, name, seq):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")


def mutate(seq, rate, rng):
    s = list(seq)
    for i in rng.sample(range(len(s)), int(len(s) * rate)):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def synthetic(size, rng):
    motifs = ["A", "AT", "AAT", "AATAT", "AAAT", "AG", "AATC"]
    parts, total = [], 0
    while total < size:
        seg = "".join(rng.choices("ACGT", weights=[0.32, 0.18, 0.18, 0.32], k=rng.randint(4000, 12000)))
        m = rng.choice(motifs)
        reps = 10 if len(m) == 1 else rng.randint(max(5, 10 // len(m)), 14)
        parts.append(seg + m * reps)
        total += len(parts[-1])
    return "".join(parts)[:size]


def measure(folder, repeats):
    runs = []
    for _ in range(repeats):
        out = tempfile.mkdtemp()
        r = subprocess.run([sys.executable, os.path.abspath(__file__), "--child", folder, out],
                           capture_output=True, text=True)
        shutil.rmtree(out, ignore_errors=True)
        if r.returncode:
            raise SystemExit("MICRO-MINE failed:\n" + r.stderr)
        runs.append([float(x) for x in r.stdout.split()])
    return (statistics.median(x[0] for x in runs), max(x[1] for x in runs),
            int(runs[0][2]), int(runs[0][3]))


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    quick = "--quick" in argv
    save = argv[argv.index("--save") + 1] if "--save" in argv else None
    if save in args:
        args.remove(save)
    if not args:
        raise SystemExit(__doc__)
    plastome = read_seq(args[0])
    repeats = int(args[1]) if len(args) > 1 else 3
    rng = random.Random(7)
    lines = []

    def emit(text=""):
        print(text, flush=True)
        lines.append(text)

    emit(f"# machine: {platform.platform()} | {platform.processor() or platform.machine()} | "
         f"{os.cpu_count()} logical CPUs | python {platform.python_version()} | repeats={repeats} "
         f"(median time, max memory)")
    emit(f"\n# A) number of genomes ({len(plastome)} bp each)")
    emit("genomes\tseconds\tpeak_MB\tSSRs\tloci")
    for n in ((2, 5, 10) if quick else (2, 5, 10, 25, 50)):
        d = tempfile.mkdtemp()
        for i in range(n):
            write(f"{d}/{i + 1:02d}_g{i + 1}.fasta", f"g{i + 1}", mutate(plastome, 0.005, rng))
        t, mb, s, c = measure(d, repeats)
        shutil.rmtree(d, ignore_errors=True)
        emit(f"{n}\t{t:.1f}\t{mb:.0f}\t{s}\t{c}")
    emit("\n# B) genome size (synthetic): annotation of one genome / comparison of two genomes")
    emit("size_Mb\tmode\tseconds\tpeak_MB\tSSRs\tloci")
    for size_mb in ((0.5, 1) if quick else (0.5, 1, 2, 3, 5)):
        base = synthetic(int(size_mb * 1e6), rng)
        d = tempfile.mkdtemp()
        write(f"{d}/01_genome.fasta", "g1", base)
        t, mb, s, c = measure(d, repeats)
        shutil.rmtree(d, ignore_errors=True)
        emit(f"{size_mb}\tannotation\t{t:.1f}\t{mb:.0f}\t{s}\t{c}")
        d = tempfile.mkdtemp()
        write(f"{d}/01_genome.fasta", "g1", base)
        write(f"{d}/02_genome.fasta", "g2", mutate(base, 0.01, rng))
        t, mb, s, c = measure(d, repeats)
        shutil.rmtree(d, ignore_errors=True)
        emit(f"{size_mb}\tcomparison(2)\t{t:.1f}\t{mb:.0f}\t{s}\t{c}")
    if save:
        with open(save, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"\nsaved to {save}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        child(sys.argv[2], sys.argv[3])
    else:
        main(sys.argv[1:])
