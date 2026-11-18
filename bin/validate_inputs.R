#!/usr/bin/env Rscript
# Validate smrna-nextflow inputs against an nf-core/smrnaseq 2.4.x run.
# Hard failures stop the pipeline; suspicious-but-tolerable states are warnings.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

option_list <- list(
  make_option("--nfcore_dir", type = "character",
              help = "nf-core/smrnaseq results directory [required]"),
  make_option("--samplesheet", type = "character",
              help = "samplesheet CSV, one row per sample [required]"),
  make_option("--contrasts", type = "character",
              help = "contrasts CSV [required]"),
  make_option("--batch_col", type = "character", default = NULL,
              help = "samplesheet column with the batch (e.g. Exp4) [optional]"),
  make_option("--expected_adapter", type = "character",
              help = "adapter the nf-core run should have declared [required]"),
  make_option("--source_dir", type = "character", default = NULL,
              help = "original nfcore_dir path recorded in the audit [optional]"),
  make_option("--outdir", type = "character", default = "results",
              help = "pipeline outdir recorded in run_manifest [optional]"),
  make_option("--container", type = "character", default = NA,
              help = "container image recorded in software_versions/run_manifest [optional]"),
  make_option("--nextflow_version", type = "character", default = NA,
              help = "Nextflow version recorded in software_versions/run_manifest [optional]")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Validate smrna-nextflow inputs."))

for (a in c("nfcore_dir", "samplesheet", "contrasts", "expected_adapter")) {
  if (is.null(opt[[a]]) || !nzchar(opt[[a]])) {
    stop(sprintf("Missing required argument --%s. Run with --help for usage.", a), call. = FALSE)
  }
}
if (!dir.exists(opt$nfcore_dir)) {
  stop(sprintf("--nfcore_dir does not exist: %s", opt$nfcore_dir), call. = FALSE)
}
if (!file.exists(opt$samplesheet)) {
  stop(sprintf("--samplesheet does not exist: %s", opt$samplesheet), call. = FALSE)
}
if (!file.exists(opt$contrasts)) {
  stop(sprintf("--contrasts does not exist: %s", opt$contrasts), call. = FALSE)
}

warn_msgs <- character(0)
warn <- function(msg) {
  warn_msgs <<- c(warn_msgs, msg)
  warning(msg, call. = FALSE, immediate. = TRUE)
}
die <- function(msg) stop(msg, call. = FALSE)

# NULL / length-0 / NA all collapse to a fallback.
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

nfcore_dir <- opt$nfcore_dir
# The audit records the path the user passed; Nextflow stages nfcore_dir under
# work/ and its basename alone is not traceable.
source_dir <- if (is.null(opt$source_dir) || !nzchar(as.character(opt$source_dir))) nfcore_dir else opt$source_dir
joined_tsv <- file.path(nfcore_dir, "mirna_quant/mirtop/joined_samples_mirtop.tsv")

# ------------------------------------------------------------------------
# 1. Required nf-core files
# ------------------------------------------------------------------------
required_files <- c(
  "mirna_quant/mirtop/joined_samples_mirtop.tsv",
  "mirna_quant/mirtop/mirna.tsv",
  "mirna_quant/edger_qc/mature_counts.csv",
  "multiqc/multiqc_data/multiqc_general_stats.txt",
  "multiqc/multiqc_data/multiqc_samtools_stats.txt"
)
missing_files <- required_files[!file.exists(file.path(nfcore_dir, required_files))]
if (length(missing_files) > 0) {
  die(paste0("Missing required nf-core files under '", nfcore_dir, "':\n  ",
             paste(missing_files, collapse = "\n  "),
             "\nExpected the results/ directory of an nf-core/smrnaseq 2.4.x run; check --nfcore_dir."))
}

params_files <- list.files(file.path(nfcore_dir, "pipeline_info"),
                           pattern = "^params_.*\\.json$", full.names = TRUE)
if (length(params_files) == 0) {
  die(paste0("No params_*.json in ", file.path(nfcore_dir, "pipeline_info"),
             ". Expected pipeline_info/params_<timestamp>.json from the nf-core run."))
}
params_file <- params_files[which.max(file.mtime(params_files))]

# ------------------------------------------------------------------------
# 2. Run parameters: adapter and ncRNA-reference coherence (warnings only)
# ------------------------------------------------------------------------
nf <- fromJSON(params_file, simplifyVector = TRUE)

genome <- as.character(nf$fasta %||% nf$genome %||% NA_character_)

# The pipeline version is not in the params dump (its 'version' key is nf-core's
# boolean param); read it from the software-versions YAML written by the run.
versions_yml <- file.path(nfcore_dir,
                          "pipeline_info/nf_core_smrnaseq_software_mqc_versions.yml")
nf_version <- NA_character_
if (file.exists(versions_yml)) {
  vlines <- grep("nf-core/smrnaseq:", readLines(versions_yml, warn = FALSE), value = TRUE)
  if (length(vlines) > 0) nf_version <- trimws(sub(".*nf-core/smrnaseq:", "", vlines[1]))
}
if (is.na(nf_version)) {
  warn(sprintf("Could not read the nf-core/smrnaseq version from %s", versions_yml))
}

declared_adapter <- as.character(nf$three_prime_adapter %||% NA_character_)
if (is.na(declared_adapter) || !identical(declared_adapter, opt$expected_adapter)) {
  warn(sprintf("Adapter mismatch: nf-core declared '%s', expected '%s'. Reads may keep an adapter base.",
               declared_adapter, opt$expected_adapter))
}
if (isTRUE(nf$with_umi)) {
  warn("nf-core ran with with_umi = true; this pipeline expects false.")
}
if (isTRUE(nf$filter_contamination)) {
  warn("nf-core ran with filter_contamination = true; rRNA/tRNA were removed before quantification, so composition is not meaningful.")
}

mir_paths <- c(
  mature    = as.character(nf$mature    %||% NA_character_),
  hairpin   = as.character(nf$hairpin   %||% NA_character_),
  mirna_gtf = as.character(nf$mirna_gtf %||% NA_character_)
)
is_url      <- grepl("^https?://", mir_paths)
is_testdata <- grepl("test-datasets", mir_paths)
if (any(is_url, na.rm = TRUE) || any(is_testdata, na.rm = TRUE) || is.na(mir_paths[["mirna_gtf"]])) {
  warn("miRBase references are not pinned by local path (URL, test-datasets or empty mirna_gtf); mature/hairpin/mirna_gtf may mix releases.")
}
mir_dirs <- dirname(mir_paths[!is.na(mir_paths)])
if (length(unique(mir_dirs)) > 1) {
  warn(sprintf("mature/hairpin/mirna_gtf are not in the same directory: %s",
               paste(unique(mir_dirs), collapse = " | ")))
}

# ------------------------------------------------------------------------
# 3. 3' adapter leftover: last base of collapsed reads
# ------------------------------------------------------------------------
fasta_files <- list.files(file.path(nfcore_dir, "mirtrace"),
                          pattern = "\\.fastp\\.fasta$", recursive = TRUE, full.names = TRUE)
last_base_note <- "not checked (no mirtrace collapsed FASTA found)"
if (length(fasta_files) == 0) {
  warn("No mirtrace collapsed FASTA found; skipped the 3' adapter-leftover check.")
} else {
  lines <- readLines(fasta_files[1], warn = FALSE)
  seqs  <- lines[!startsWith(lines, ">")]
  seqs  <- seqs[nzchar(seqs)]
  if (length(seqs) > 0) {
    tab  <- sort(table(substr(seqs, nchar(seqs), nchar(seqs))), decreasing = TRUE)
    frac <- as.numeric(tab[1]) / sum(tab)
    last_base_note <- sprintf("%s in %.1f%% of reads (%s)", names(tab)[1], 100 * frac, basename(fasta_files[1]))
    # 90% cutoff fixed; not an analysis parameter.
    if (frac > 0.9) {
      warn(sprintf("Last base is '%s' in %.1f%% of reads; likely adapter leftover.",
                   names(tab)[1], 100 * frac))
    }
  }
}

# ------------------------------------------------------------------------
# 4a. samplesheet.csv
# ------------------------------------------------------------------------
meta <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
required_meta_cols <- c("sample_id", "condition", "norm_group")
missing_meta_cols  <- setdiff(required_meta_cols, names(meta))
if (length(missing_meta_cols) > 0) {
  die(sprintf("samplesheet.csv is missing required columns: %s",
              paste(missing_meta_cols, collapse = ", ")))
}
if (anyDuplicated(meta$sample_id) > 0) {
  dups <- unique(meta$sample_id[duplicated(meta$sample_id)])
  die(sprintf("samplesheet.csv has duplicated sample_id: %s", paste(dups, collapse = ", ")))
}

hdr <- strsplit(readLines(joined_tsv, n = 1), "\t")[[1]]
fixed_cols <- c("UID", "Read", "miRNA", "Variant", "iso_5p", "iso_3p", "iso_add3p",
                "iso_snp", "iso_5p_nt", "iso_3p_nt", "iso_add3p_nt", "iso_snp_nt")
sample_cols <- setdiff(hdr, fixed_cols)

missing_in_counts <- setdiff(meta$sample_id, sample_cols)
if (length(missing_in_counts) > 0) {
  die(paste0("samplesheet samples are not in the mirtop columns:",
             "\n  in samplesheet but not in counts: ", paste(missing_in_counts, collapse = ", "),
             "\nFix samplesheet.csv so sample_id matches the nf-core sample names exactly."))
}
# nf-core samples absent from samplesheet.csv are simply not analysed: the
# samplesheet selects the subset of interest, so this is not reported.

missing_evidence <- character(0)
for (s in meta$sample_id) {
  rt <- file.path(nfcore_dir, "mirtrace", s, "mirtrace-stats-rnatype.tsv")
  fa <- file.path(nfcore_dir, "mirtrace", s, "qc_passed_reads.all.collapsed", paste0(s, ".fastp.fasta"))
  if (!file.exists(rt)) missing_evidence <- c(missing_evidence, rt)
  if (!file.exists(fa)) missing_evidence <- c(missing_evidence, fa)
}
if (length(missing_evidence) > 0) {
  die(paste0("Missing per-sample mirtrace files:\n  ",
             paste(missing_evidence, collapse = "\n  "),
             "\nExpected mirtrace/<sample>/... for every sample in samplesheet.csv."))
}

# ------------------------------------------------------------------------
# 4b. contrasts.csv
# ------------------------------------------------------------------------
con <- read.csv(opt$contrasts, stringsAsFactors = FALSE, check.names = FALSE)
required_con_cols <- c("contrast_id", "variable", "reference", "treatment", "norm_group")
missing_con_cols  <- setdiff(required_con_cols, names(con))
if (length(missing_con_cols) > 0) {
  die(sprintf("contrasts.csv is missing required columns: %s",
              paste(missing_con_cols, collapse = ", ")))
}
if (anyDuplicated(con$contrast_id) > 0) {
  dups <- unique(con$contrast_id[duplicated(con$contrast_id)])
  die(sprintf("contrasts.csv has duplicated contrast_id: %s", paste(dups, collapse = ", ")))
}

contrast_errors <- character(0)
keep_con <- rep(TRUE, nrow(con))
for (i in seq_len(nrow(con))) {
  r <- con[i, ]
  if (!r$variable %in% names(meta)) {
    contrast_errors <- c(contrast_errors,
      sprintf("'%s': variable '%s' is not a column of samplesheet.csv", r$contrast_id, r$variable))
    next
  }
  if (identical(as.character(r$reference), as.character(r$treatment))) {
    contrast_errors <- c(contrast_errors,
      sprintf("'%s': reference and treatment are identical", r$contrast_id))
    next
  }

  sub <- meta[meta$norm_group == r$norm_group, , drop = FALSE]
  if (nrow(sub) == 0) {
    keep_con[i] <- FALSE
    warn(sprintf("Contrast '%s' dropped: norm_group '%s' has no samples in samplesheet.csv.",
                 r$contrast_id, r$norm_group))
    next
  }
  vals  <- sub[[r$variable]]
  n_ref <- sum(vals == r$reference)
  n_trt <- sum(vals == r$treatment)
  if (n_ref < 2 || n_trt < 2) {
    # A contrast whose levels are not in the (possibly subset) samplesheet is
    # dropped with a warning, not a hard error, so a subset run is possible.
    keep_con[i] <- FALSE
    warn(sprintf("Contrast '%s' dropped: needs >= 2 samples per level in '%s' (ref=%d, trt=%d).",
                 r$contrast_id, r$norm_group, n_ref, n_trt))
  } else if (n_ref < 3 || n_trt < 3) {
    warn(sprintf("Contrast '%s': small group (ref=%d, trt=%d); interpret with care.",
                 r$contrast_id, n_ref, n_trt))
  }
}
if (length(contrast_errors) > 0) {
  die(paste0("Invalid contrasts.csv:\n  ", paste(contrast_errors, collapse = "\n  ")))
}
if (!all(keep_con)) {
  con <- con[keep_con, , drop = FALSE]
  if (nrow(con) == 0) {
    die("No usable contrast remains after applying samplesheet.csv (all levels missing or < 2 samples).")
  }
}

# ------------------------------------------------------------------------
# 4c. batch column: complete or empty per group; design decided per group
# ------------------------------------------------------------------------
# Design is chosen per norm_group: ~ <batch_col> + condition when the group has
# >= 2 non-empty batch levels and the design matrix is full rank, else
# ~ condition. A partially filled batch column inside a group is an error.
group_design <- function(sub) {
  if (is.null(opt$batch_col) || !nzchar(opt$batch_col) || !opt$batch_col %in% names(sub)) {
    return(list(usable = FALSE, formula = "~ condition"))
  }
  v <- as.character(sub[[opt$batch_col]])
  v[is.na(v)] <- ""
  n_empty <- sum(!nzchar(v))
  if (n_empty > 0 && n_empty < length(v)) {
    die(sprintf("batch column '%s' is partially empty in norm_group '%s' (%d of %d). Fill every sample in the group or leave the whole group empty.",
                opt$batch_col, as.character(sub$norm_group[1]), n_empty, length(v)))
  }
  if (n_empty == length(v) || length(unique(v)) < 2) {
    return(list(usable = FALSE, formula = "~ condition"))
  }
  mm <- model.matrix(~ v + condition,
                     data = data.frame(v = factor(v), condition = factor(sub$condition)))
  if (qr(mm)$rank < ncol(mm)) {
    warn(sprintf("batch column '%s' is confounded with condition in norm_group '%s' (rank-deficient); using ~ condition.",
                 opt$batch_col, as.character(sub$norm_group[1])))
    return(list(usable = FALSE, formula = "~ condition"))
  }
  list(usable = TRUE, formula = sprintf("~ %s + condition", opt$batch_col))
}
group_levels <- unique(as.character(meta$norm_group))
group_designs <- setNames(
  lapply(group_levels, function(g) group_design(meta[meta$norm_group == g, , drop = FALSE])),
  group_levels)

# ------------------------------------------------------------------------
# Outputs (published to run_info/ by Nextflow)
# ------------------------------------------------------------------------
as_str <- function(x) if (is.na(x)) "NA" else x
count_lines <- function(x) {
  tab <- table(x)
  paste0("  - ", names(tab), ": ", as.integer(tab))
}

audit <- c(
  "# Audit",
  "",
  sprintf("- nf-core/smrnaseq version: %s", as_str(nf_version)),
  sprintf("- nf-core results dir: %s", source_dir),
  sprintf("- params file used: %s", params_file),
  sprintf("- genome FASTA: %s", as_str(genome)),
  sprintf("- declared adapter: %s (expected %s)", as_str(declared_adapter), opt$expected_adapter),
  sprintf("- 3' end check: %s", last_base_note),
  "",
  "## miRBase references",
  sprintf("- mature: %s", as_str(mir_paths[["mature"]])),
  sprintf("- hairpin: %s", as_str(mir_paths[["hairpin"]])),
  sprintf("- mirna_gtf: %s", as_str(mir_paths[["mirna_gtf"]])),
  "",
  "## Samples",
  sprintf("- total: %d", nrow(meta)),
  "- by condition:",
  count_lines(meta$condition),
  "- by norm_group:",
  count_lines(meta$norm_group),
  "",
  "## Design",
  if (is.null(opt$batch_col) || !nzchar(opt$batch_col)) "- batch_col not set: ~ condition for every group"
  else paste0("- ", names(group_designs), ": ", vapply(group_designs, function(d) d$formula, character(1))),
  "",
  "## Warnings",
  if (length(warn_msgs) == 0) "- none" else paste0("- ", warn_msgs),
  ""
)
writeLines(audit, "audit.md")

write.csv(meta, "samplesheet_validated.csv", row.names = FALSE, quote = FALSE, na = "")
write.csv(con,  "contrasts_validated.csv", row.names = FALSE, quote = FALSE, na = "")

# ------------------------------------------------------------------------
# run_info extras
# ------------------------------------------------------------------------
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE,
                                           row.names = FALSE, na = "")

# sample_index.tsv: the sample order used by the pipeline (= samplesheet.csv
# order), plus the full samplesheet. nf-core samples absent from samplesheet.csv
# are ignored.
sample_index <- data.frame(index = seq_len(nrow(meta)), meta, stringsAsFactors = FALSE)
write_tsv(sample_index, "sample_index.tsv")

# design_check/<group>.tsv: one row per level of the contrast variable(s) and
# the batch column, plus the design chosen for that group.
dir.create("design_check", recursive = TRUE, showWarnings = FALSE)
for (g in unique(meta$norm_group)) {
  sub    <- meta[meta$norm_group == g, , drop = FALSE]
  design <- group_designs[[as.character(g)]]$formula
  vars   <- unique(as.character(con$variable[con$norm_group == g]))
  if (length(vars) == 0) vars <- "condition"
  if (!is.null(opt$batch_col) && nzchar(opt$batch_col) && opt$batch_col %in% names(sub)) {
    vars <- unique(c(vars, opt$batch_col))
  }
  rows <- list()
  for (v in vars) {
    if (!v %in% names(sub)) next
    vals <- as.character(sub[[v]])
    vals[is.na(vals)] <- ""
    lvls <- unique(vals)
    if (identical(v, opt$batch_col)) lvls <- lvls[nzchar(lvls)]
    if (length(lvls) == 0) lvls <- "(empty)"
    for (lv in lvls) {
      in_lv <- if (identical(lv, "(empty)")) !nzchar(vals) else vals == lv
      rows[[length(rows) + 1]] <- data.frame(
        variable    = v,
        level       = lv,
        n_samples   = sum(in_lv),
        sample_ids  = paste(sub$sample_id[in_lv], collapse = ","),
        main_design = design,
        stringsAsFactors = FALSE
      )
    }
  }
  dd <- do.call(rbind, rows)
  write_tsv(dd, file.path("design_check", paste0(g, ".tsv")))
}

# software_versions.tsv: R and key packages + the nf-core/Nextflow/container.
pkgs <- c("DESeq2", "data.table", "ggplot2", "ggrepel", "pheatmap", "optparse",
          "jsonlite", "yaml", "rmarkdown", "knitr", "clusterProfiler")
pkg_ver <- vapply(pkgs, function(p)
  tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_),
  character(1))
software_versions <- data.frame(
  software = c("R", pkgs, "nf-core/smrnaseq", "nextflow", "container"),
  version  = c(as.character(getRversion()), unname(pkg_ver),
               as_str(nf_version), as_str(opt$nextflow_version), as_str(opt$container)),
  stringsAsFactors = FALSE
)
write_tsv(software_versions, "software_versions.tsv")

# run_manifest.tsv: what produced this run and where.
run_manifest <- data.frame(
  key   = c("run_date", "outdir", "nfcore_dir", "nfcore_version",
            "nextflow_version", "container"),
  value = c(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), as_str(opt$outdir), source_dir,
            as_str(nf_version), as_str(opt$nextflow_version), as_str(opt$container)),
  stringsAsFactors = FALSE
)
write_tsv(run_manifest, "run_manifest.tsv")

message(sprintf("validate_inputs ok: %d samples, %d contrasts, %d warning(s).",
                nrow(meta), nrow(con), length(warn_msgs)))
