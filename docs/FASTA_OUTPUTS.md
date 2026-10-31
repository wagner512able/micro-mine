# FASTA outputs (MICRO-MINE)

Every record is `flank + repeat tract + flank`. The **repeat tract is highlighted by case** (`--ssr-case`, default: tract in
lowercase, flanks in UPPERCASE). The tract is the full repeat units **plus the incomplete unit that continues the tract**
(for example `(TC)6` followed by `T` is highlighted as `tctctctctctct`). Repeat number, length and coordinates (table,
GFF3) count full units only, as MISA does, so results stay comparable; the extra base(s) are listed in the report column
*Partial unit* and in the header field `partial=`. The right flank starts after the partial unit.

## Header fields
`>SSR.012|01_species|contig:start-end|(AAG)6|repeats=6|flank=150/150|ssr_in_seq=151-168|partial=T|WARNING_N=3`

| Field | Meaning |
|---|---|
| `(AAG)6` and `repeats=6` | motif as found and number of full repeat units |
| `ssr_in_seq=151-168` | position of the whole highlighted tract inside the sequence (1-based); use it as the primer3 target |
| `partial=T` | incomplete unit included in the highlighted tract (only when present) |
| `WARNING_N=3` | N bases in the flanks (unreliable boundaries/primers) |
| `repeats_by_species=01:6,02:7,03:3(below_threshold)` | consolidated file: repeat number of the locus in every species |
| `exclusive=real` / `threshold_only` | exclusive-SSR files: kind of exclusivity |
| `locus_copies=1` | copies of the locus in that species (genome in single-genome mode) |
| `strand=-` | per-locus files: the copy was reverse-complemented to match the first copy |

## Files
* `fasta_per_species/` - all SSRs of each species.
* `fasta_per_locus/` - one file per locus, copies from all species oriented alike.
* `SSRs_with_flanks.fasta` - one record per locus.
* `fasta_exclusive/` (`--exclusive-fasta species|contig`) - **only the exclusive SSRs**, one file per species (or per
  species and contig). Real exclusive = present in one species and with no homologous locus, even below the threshold,
  in the others. `--exclusive-single-copy` keeps only loci with one copy; `--exclusive-threshold` also adds SSRs that are
  exclusive only because of the threshold.
* Single-genome mode: `SSRs_unique_single_copy.fasta` (or `fasta_unique/` per contig) - SSRs whose locus has a single
  copy in the whole genome (no similar locus elsewhere, e.g. outside repeated regions): the safest primer candidates.
