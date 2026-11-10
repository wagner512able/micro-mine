# Benchmarks

Both scripts run the real MICRO-MINE code (`micromine.cli.run`) and print tab-separated tables.

* `scaling.py PLASTOME.fasta [REPEATS] [--quick] [--save FILE]` - self-contained script: run time and peak memory
  against the number of genomes (2-50 copies of a real chloroplast genome with 0.5 % substitutions) and against genome
  size (synthetic genomes of 0.5-5 Mb). The first output line records the machine (OS, CPU, Python): report it with the
  results.
* `synthetic_divergence.py [SEED]` - fraction of SSR loci that are still recognised as the same locus when the flanks
  of the second genome diverge by 0-20 % (default similarity 0.6 and 0.4). Deterministic for a given seed and version.
* `results/` - outputs obtained in a Linux sandbox (1 vCPU); rerun on your own hardware before citing numbers.

Run time varies with the load of the machine: use several repeats and report the median.
