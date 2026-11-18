# deliverables/ — what to open and how to read it

Tables and figures for the wet lab. Everything here is derived from an
nf-core/smrnaseq 2.4.x run; see `../run_info/` for provenance.

## Start here

1. `report/smrna_report.html` — the whole run in one page.
2. `qc/qc_table.tsv` and `qc/flagged_samples.tsv` — per-sample QC. Samples are
   **flagged, never excluded**; flags are advisory.
3. `de/de_summary.tsv` — how many miRNAs changed per contrast.
4. `de/<contrast_id>/de_significant.tsv` — the hits.

## Folders

- `counts/counts_mirtop.tsv` — miRNA x sample matrix. Multi-assigned reads are
  split 1/k across their miRNAs and rounded; first column `mirna`, then one
  column per sample in `samplesheet.csv` order.
- `qc/` — `qc_table.tsv` (all metrics), `qc_flags.tsv` / `flagged_samples.txt`
  (which sample failed which check), `composition_by_sample.tsv` and
  `composition_delta.tsv` (mirTrace RNA-type composition), `plots/`.
- `explore/<norm_group>/` — PCA coordinates and sample-distance heatmap.
- `de/<contrast_id>/` — `de_full.tsv` (every tested miRNA) and
  `de_significant.tsv` (`padj < alpha`), plus `plots/volcano.pdf`, `.png`.
- `probes/` — optional (read probes): `probe_summary.tsv` (reads derived from given
  probe sequences), `probe_hits.tsv` and `top_unknown_by_sample.tsv`.

## Column dictionary — `de_*.tsv`

| column | meaning |
| --- | --- |
| `mirna` | mature miRNA name |
| `baseMean` | mean of normalised counts across the contrast samples |
| `log2FC_mle` | DESeq2 MLE log2 fold change |
| `log2FC_shrunk` | shrunken log2FC (used for ranking and `sig_lfc1`) |
| `lfcSE`, `stat`, `pvalue`, `padj` | Wald test on the MLE |
| `direction` | `up`/`down`/`ns` (`padj < alpha` and sign of `log2FC_shrunk`) |
| `sig_lfc1` | `padj < alpha` **and** `|log2FC_shrunk| >= lfc_flag` (significance + effect) |
| `n_ref`, `n_trt` | samples per level |
| `mean_norm_count_ref`, `mean_norm_count_trt` | normalised means per level |

**Sign convention:** `log2FC = treatment / reference`. Positive = more abundant
in the `treatment` level of `contrasts.csv`.

## Column dictionary — `qc_table.tsv`

| column | meaning |
| --- | --- |
| `reads_raw`, `reads_passed_fastp`, `pct_passed` | fastp input/output |
| `mirtrace_pct_<type>` | mirTrace RNA-type composition (percent, sums to 100) |
| `mirtop_raw_total` | mirtop reads before the isomiR filter |
| `mirna_reads_mirtop` | mirtop reads after the isomiR filter (used by `LOW_DEPTH`) |
| `mature_total` | reads in the mature-miRNA (edgeR) counts |
| `canonical_frac`, `multi_assigned_frac` | fraction canonical / split across >1 miRNA |
| `genome_mapped_pct` | reads mapped to the genome (on `*_mature_hairpin_genome`) |

## Notes

- Composition percentages are relative to the small-RNA content of each sample
  (compositional, not absolute).
- DE is relative to the total miRNA fraction, not absolute abundance.
- Small n: interpret each contrast with care; see `de/concordance.tsv`.
