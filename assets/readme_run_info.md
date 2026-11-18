# run_info/ — provenance and audit trail

Everything needed to reproduce and interpret a run. Written before and during
the run; nothing here is an analysis result.

| file | contents |
| --- | --- |
| `audit.md` | audit report: nf-core version, adapter check, 3' leftover check, miRBase references, sample/group counts, design, warnings |
| `samplesheet_validated.csv`, `contrasts_validated.csv` | the inputs after validation |
| `sample_index.tsv` | sample order as used (guards against silent column swaps) + samplesheet |
| `design_check/<norm_group>.tsv` | levels of the contrast variable and the batch column, plus the design chosen per group (`main_design`) |
| `params_used.tsv` | every parameter with its effective value |
| `thresholds_used.tsv` | the analysis thresholds (QC, expression filter, DE) |
| `software_versions.tsv` | R/packages, nf-core/smrnaseq, Nextflow and the container image |
| `run_manifest.tsv` | run date, outdir, nf-core dir and version, Nextflow version, container |
| `qa_report.tsv` | per-sample QC pass/flag summary |
| `pipeline_info/` | `params.json`, Nextflow `execution_report.html` and `execution_trace.txt` |

No threshold or path is hardcoded in the code: everything in `params_used.tsv` /
`thresholds_used.tsv` comes from `nextflow.config` (defaults) and the run's
`params.yaml` (overrides).
