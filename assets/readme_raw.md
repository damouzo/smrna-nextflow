# raw/ — native dumps for internal audit

Intermediate and per-tool dumps kept for traceability. **Not needed to read the
results** — use `../deliverables/` for that.

- `counts/counts_mature.tsv` — mature-miRNA counts from the nf-core edgeR table
  (sensitivity source, `count_source: mature`).
- `counts/isomir_filtered.tsv.gz` — the mirtop isomiR table after the isomiR
  filter, with the 1/k weight per (UID, miRNA). Input of the isomiR phase.
- `deseq2/<norm_group>/` — `dds.rds` (DESeqDataSet), `size_factors.tsv`,
  `filtered_counts.tsv` after the expression filter.
- `deseq2/sens_*/<contrast_id>/` — sensitivity runs (`sens_condition`,
  `sens_mature`, `sens_no_flagged`) as `de_full.tsv` / `de_significant.tsv`.
- `biotype/` — phase 2 alignment stats and feature counts (only with
  `run_biotype`).
