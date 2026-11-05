# smrna-nextflow

Downstream analysis of small RNA-seq on top of [nf-core/smrnaseq](https://nf-co.re/smrnaseq) 2.4.x results: count matrix, QC flags, biotype composition, PCA, differential expression (DESeq2) and an HTML report. Optional modules: full biotype composition, isomiRs and miRNA target enrichment.

## Requirements

- Nextflow >= 24.10
- Singularity/Apptainer or Docker (the image is built by GitHub Actions and published to GHCR)

## Quick start

```bash
nextflow run . -profile singularity,slurm -params-file params.yaml
```

```yaml
# params.yaml
nfcore_dir: /path/to/nfcore_smrnaseq/results
metadata:   metadata.csv
contrasts:  contrasts.csv
outdir:     results_smrna
```

Generic templates are in `assets/`; a working example is in `demo/batch_2026_04/`.

## Inputs

**metadata.csv** — one row per sample.

| column | required | description |
| --- | --- | --- |
| `sample_id` | yes | same name as in nf-core |
| `condition` | yes | label used in `contrasts.csv` |
| `norm_group` | yes | samples sharing a model and normalisation |
| `exp`, `subject`, `library_id` | no | batch, donor, library id |

**contrasts.csv** — one row per contrast; log2FC = treatment / reference.

```csv
contrast_id,variable,reference,treatment,norm_group
treated_vs_control,condition,control,treated,groupA
```

## Main parameters

| parameter | default | description |
| --- | --- | --- |
| `batch_col` | null | metadata column for an extra `~ batch + condition` run |
| `count_source` | `mirtop` | `mirtop` or `mature` |
| `alpha`, `lfc_flag` | 0.05, 1 | FDR and strong-effect log2FC threshold |
| `run_biotype`, `run_isomir`, `run_functional` | false | enable optional modules |

All defaults live in `nextflow.config`.

## Optional resources

Large files are never downloaded by the pipeline. Download them once and set the path in `params.yaml` or `nextflow.config`.

```bash
# ncrna_bed (run_biotype): tRNA, snoRNA, rRNA, miRNA, snRNA... coordinates (skip if you already have it)
wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/genome_coordinates/bed/homo_sapiens.GRCh38.bed.gz
gunzip homo_sapiens.GRCh38.bed.gz

# genome_fasta and gtf (run_biotype): same Ensembl release for both
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

# bowtie_index (run_biotype, optional): built from genome_fasta if left empty
bowtie-build --threads 8 <(zcat Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz) GRCh38
```

`mirtarbase_file` (`run_functional`): download the Homo sapiens MTI file from <https://mirtarbase.cuhk.edu.cn/miRTarBase/downloads/>.

## Upstream note

Set `three_prime_adapter` in the nf-core run to the adapter of your kit. A wrong adapter leaves one adapter base on most reads and distorts isomiR and canonical counts; the pipeline warns about it. Pass `mature`, `hairpin` and `mirna_gtf` from the same miRBase release, pinned by path (e.g. `.../mirbase/v23.1/`), never by a `CURRENT` URL. To use a newer release, download it into a new `mirbase/vX.Y/` folder and change the three paths; results will differ because miRNA names and sequences change between releases.

## Outputs

```
results_smrna/
├── 00_audit/     validation report
├── 01_counts/    miRNA count matrices
├── 02_qc/        QC table, flags, composition
├── 03_explore/   PCA and sample distances
├── 04_de/        DESeq2 results per run and contrast
├── 05_biotype/   optional
├── 06_isomir/    optional
├── 07_functional/ optional
└── report/       smrna_report.html
```
