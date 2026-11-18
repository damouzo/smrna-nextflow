#!/usr/bin/env Rscript
# Build the miRNA count matrix and the QC table from an nf-core/smrnaseq run.
# Multi-assigned reads are split 1/k across their miRNAs, then rounded.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
})

option_list <- list(
  make_option("--nfcore_dir", type = "character",
              help = "nf-core/smrnaseq results directory [required]"),
  make_option("--samplesheet", type = "character",
              help = "samplesheet_validated.csv [required]"),
  make_option("--isomir_exclude_variants", type = "character",
              help = "comma-separated Variant tags to drop [required]"),
  make_option("--isomir_max_shift_5p", type = "integer",
              help = "max absolute 5' shift kept [required]"),
  make_option("--isomir_max_shift_3p", type = "integer",
              help = "max absolute 3' shift kept [required]"),
  make_option("--isomir_min_len", type = "integer",
              help = "min read length kept [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Build the miRNA count matrix and QC table."))

for (a in c("nfcore_dir", "samplesheet", "isomir_exclude_variants",
            "isomir_max_shift_5p", "isomir_max_shift_3p", "isomir_min_len")) {
  if (is.null(opt[[a]]) || !nzchar(as.character(opt[[a]]))) {
    stop(sprintf("Missing required argument --%s. Run with --help for usage.", a), call. = FALSE)
  }
}
if (!dir.exists(opt$nfcore_dir)) {
  stop(sprintf("--nfcore_dir does not exist: %s", opt$nfcore_dir), call. = FALSE)
}
if (!file.exists(opt$samplesheet)) {
  stop(sprintf("--samplesheet does not exist: %s", opt$samplesheet), call. = FALSE)
}

die <- function(msg) stop(msg, call. = FALSE)
excl_tags <- trimws(strsplit(opt$isomir_exclude_variants, ",", fixed = TRUE)[[1]])

nfcore_dir    <- opt$nfcore_dir
joined_tsv    <- file.path(nfcore_dir, "mirna_quant/mirtop/joined_samples_mirtop.tsv")
mature_csv    <- file.path(nfcore_dir, "mirna_quant/edger_qc/mature_counts.csv")
samtools_tsv  <- file.path(nfcore_dir, "multiqc/multiqc_data/multiqc_samtools_stats.txt")
fastp_dir     <- file.path(nfcore_dir, "fastp/on_raw")
mirtrace_dir  <- file.path(nfcore_dir, "mirtrace")

for (p in c(joined_tsv, mature_csv, samtools_tsv)) {
  if (!file.exists(p)) die(sprintf("Missing nf-core file: %s", p))
}
if (!dir.exists(fastp_dir))    die(sprintf("Missing fastp directory: %s", fastp_dir))
if (!dir.exists(mirtrace_dir)) die(sprintf("Missing mirtrace directory: %s", mirtrace_dir))

meta    <- fread(opt$samplesheet)
samples <- as.character(meta$sample_id)

# ------------------------------------------------------------------------
# mirtop counts: filter isomiRs, split 1/k, round
# ------------------------------------------------------------------------
dt <- fread(joined_tsv, sep = "\t", header = TRUE)
ann_cols <- c("UID", "Read", "miRNA", "Variant",
              "iso_5p", "iso_3p", "iso_add3p", "iso_snp")
missing_cols <- setdiff(c(ann_cols, samples), names(dt))
if (length(missing_cols) > 0) {
  die(paste0("joined_samples_mirtop.tsv is missing columns: ",
             paste(missing_cols, collapse = ", "),
             "\nCheck that samplesheet sample_id matches the nf-core sample names."))
}

m_full <- as.matrix(dt[, ..samples])
mirtop_raw_total <- colSums(m_full)

has_excl <- vapply(strsplit(dt$Variant, ",", fixed = TRUE),
                   function(tags) any(tags %in% excl_tags), logical(1))
read_len <- nchar(dt$Read)
keep <- !has_excl &
  abs(as.integer(dt$iso_5p)) <= opt$isomir_max_shift_5p &
  abs(as.integer(dt$iso_3p)) <= opt$isomir_max_shift_3p &
  read_len >= opt$isomir_min_len
dt_f <- dt[keep]

removed_pct <- 100 * (mirtop_raw_total - colSums(as.matrix(dt_f[, ..samples]))) / mirtop_raw_total
message(sprintf("isomiR filter removed %.1f%%-%.1f%% of mirtop reads per sample.",
                min(removed_pct), max(removed_pct)))

dt_f[, k := uniqueN(miRNA), by = UID]
contrib <- as.matrix(dt_f[, ..samples]) * (1 / dt_f$k)

agg_raw        <- rowsum(contrib, group = dt_f$miRNA)
counts_mirtop  <- round(agg_raw)                       # miRNA x sample
# rowsum keeps first-appearance order; sort by name so reruns diff cleanly.
counts_mirtop  <- counts_mirtop[order(rownames(counts_mirtop)), , drop = FALSE]
read_tot       <- colSums(contrib)

canon <- as.integer(dt_f$iso_5p) == 0 & as.integer(dt_f$iso_3p) == 0 &
  as.integer(dt_f$iso_add3p) == 0 & as.integer(dt_f$iso_snp) == 0
canonical_frac <- colSums(contrib[canon, , drop = FALSE]) / read_tot
multi_frac     <- colSums(contrib[dt_f$k > 1, , drop = FALSE]) / read_tot

# ------------------------------------------------------------------------
# mature source (sensitivity): mature_counts.csv is samples x miRNA
# ------------------------------------------------------------------------
mature_dt   <- fread(mature_csv)
mature_samp <- as.character(mature_dt[[1]])
mature_mat  <- as.matrix(mature_dt[, -1, with = FALSE])
rownames(mature_mat) <- mature_samp
missing_mature <- setdiff(samples, mature_samp)
if (length(missing_mature) > 0) {
  die(paste0("mature_counts.csv is missing samples: ", paste(missing_mature, collapse = ", ")))
}
mature_mat     <- mature_mat[samples, , drop = FALSE]
mature_total   <- rowSums(mature_mat)
counts_mature  <- t(mature_mat)                        # miRNA x sample

unmatched <- setdiff(rownames(counts_mature), rownames(counts_mirtop))
message(sprintf("mature miRNAs not present in mirtop: %d (kept, sensitivity only).",
                length(unmatched)))

# ------------------------------------------------------------------------
# QC table: fastp, mirTrace, samtools + mirtop/mature totals
# ------------------------------------------------------------------------
fastp_reads <- function(s) {
  p <- file.path(fastp_dir, paste0(s, ".fastp.json"))
  if (!file.exists(p)) die(sprintf("Missing fastp report: %s", p))
  j <- fromJSON(p)
  c(raw = j$summary$before_filtering$total_reads,
    passed = j$summary$after_filtering$total_reads)
}
fp <- vapply(samples, fastp_reads, numeric(2))         # 2 x n
reads_raw    <- fp["raw", ]
reads_passed <- fp["passed", ]

st <- fread(samtools_tsv)
st_gen <- st[grepl("_mature_hairpin_genome$", Sample)]
genome_pct <- setNames(as.numeric(st_gen$reads_mapped_percent),
                       sub("_mature_hairpin_genome$", "", st_gen$Sample))

rnatype_pct <- function(s) {
  p <- file.path(mirtrace_dir, s, "mirtrace-stats-rnatype.tsv")
  if (!file.exists(p)) die(sprintf("Missing mirTrace rnatype file: %s", p))
  rt <- fread(p)
  setnames(rt, names(rt), c("RNA_TYPE", "value"))
  v <- setNames(as.numeric(rt$value), rt$RNA_TYPE)
  cats <- c("miRNA", "tRNA", "rRNA", "artifact", "unknown")
  # percentage of each category; missing categories become 0
  sapply(cats, function(k) if (k %in% names(v)) 100 * v[[k]] / sum(v) else 0)
}
mt <- vapply(samples, rnatype_pct, numeric(5))         # 5 x n

qc <- data.table(
  sample_id             = samples,
  condition             = as.character(meta$condition),
  norm_group            = as.character(meta$norm_group),
  batch                 = if ("batch" %in% names(meta)) as.character(meta$batch) else NA_character_,
  reads_raw             = as.numeric(reads_raw),
  reads_passed_fastp    = as.numeric(reads_passed),
  pct_passed            = 100 * reads_passed / reads_raw,
  mirtrace_pct_miRNA    = as.numeric(mt["miRNA", ]),
  mirtrace_pct_tRNA     = as.numeric(mt["tRNA", ]),
  mirtrace_pct_rRNA     = as.numeric(mt["rRNA", ]),
  mirtrace_pct_artifact = as.numeric(mt["artifact", ]),
  mirtrace_pct_unknown  = as.numeric(mt["unknown", ]),
  mirtop_raw_total      = as.numeric(mirtop_raw_total),
  mirna_reads_mirtop    = as.numeric(read_tot),
  mature_total          = as.numeric(mature_total),
  canonical_frac        = as.numeric(canonical_frac),
  multi_assigned_frac   = as.numeric(multi_frac),
  genome_mapped_pct     = as.numeric(genome_pct[samples])
)

# ------------------------------------------------------------------------
# Outputs (published by Nextflow: deliverables/ and raw/ subtrees)
# ------------------------------------------------------------------------
dir.create("deliverables/counts", recursive = TRUE, showWarnings = FALSE)
dir.create("deliverables/qc",     recursive = TRUE, showWarnings = FALSE)
dir.create("raw/counts",          recursive = TRUE, showWarnings = FALSE)

fwrite(qc, "deliverables/qc/qc_table.tsv", sep = "\t", quote = FALSE, na = "")

counts_out <- data.table(mirna = rownames(counts_mirtop), as.data.table(counts_mirtop))
setnames(counts_out, c("mirna", samples))
fwrite(counts_out, "deliverables/counts/counts_mirtop.tsv", sep = "\t", quote = FALSE)

mature_out <- data.table(mirna = rownames(counts_mature), as.data.table(counts_mature))
setnames(mature_out, c("mirna", samples))
fwrite(mature_out, "raw/counts/counts_mature.tsv", sep = "\t", quote = FALSE)

isomir_out <- data.table(
  UID = dt_f$UID, Read = dt_f$Read, miRNA = dt_f$miRNA, Variant = dt_f$Variant,
  iso_5p = dt_f$iso_5p, iso_3p = dt_f$iso_3p, iso_add3p = dt_f$iso_add3p, iso_snp = dt_f$iso_snp,
  k = dt_f$k, weight = 1 / dt_f$k
)
isomir_out <- cbind(isomir_out, dt_f[, ..samples])
fwrite(isomir_out, "raw/counts/isomir_filtered.tsv.gz",
       sep = "\t", quote = FALSE, compress = "gzip")

message(sprintf("build_counts ok: %d miRNAs x %d samples; %d mature miRNAs in raw/.",
                nrow(counts_mirtop), length(samples), nrow(counts_mature)))
