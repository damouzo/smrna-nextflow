#!/usr/bin/env Rscript
# QC flags and mirTrace composition. Samples are flagged, never excluded.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

option_list <- list(
  make_option("--qc_table", type = "character", help = "qc_table.tsv [required]"),
  make_option("--samplesheet", type = "character", help = "samplesheet_validated.csv [required]"),
  make_option("--contrasts", type = "character", help = "contrasts_validated.csv [required]"),
  make_option("--qc_min_mirna_frac", type = "double", help = "LOW_MIRNA_FRAC threshold, fraction [required]"),
  make_option("--qc_max_trna_frac", type = "double", help = "HIGH_TRNA_FRAC threshold, fraction [required]"),
  make_option("--qc_max_rrna_frac", type = "double", help = "HIGH_RRNA_FRAC threshold, fraction [required]"),
  make_option("--qc_min_mirna_reads", type = "double", help = "LOW_DEPTH threshold, miRNA reads [required]"),
  make_option("--seed", type = "integer", help = "random seed [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "QC flags and mirTrace composition."))

for (a in c("qc_table", "samplesheet", "contrasts", "qc_min_mirna_frac",
            "qc_max_trna_frac", "qc_max_rrna_frac", "qc_min_mirna_reads", "seed")) {
  if (is.null(opt[[a]])) stop(sprintf("Missing required argument --%s. Run with --help.", a), call. = FALSE)
}
for (p in c(opt$qc_table, opt$samplesheet, opt$contrasts)) {
  if (!file.exists(p)) stop(sprintf("File not found: %s", p), call. = FALSE)
}
set.seed(opt$seed)

qc   <- fread(opt$qc_table)
meta <- fread(opt$samplesheet)
con  <- fread(opt$contrasts)

categories <- c("miRNA", "tRNA", "rRNA", "artifact", "unknown")
pct_cols   <- paste0("mirtrace_pct_", categories)
stopifnot(all(c("sample_id", pct_cols, "mirna_reads_mirtop") %in% names(qc)))

qc[, sample_id := factor(sample_id, levels = sample_id)]

# ------------------------------------------------------------------------
# Flags (marked, not excluded)
# ------------------------------------------------------------------------
flag_specs <- list(
  LOW_MIRNA_FRAC = list(value = qc$mirtrace_pct_miRNA / 100, threshold = opt$qc_min_mirna_frac, lower = TRUE),
  HIGH_TRNA_FRAC = list(value = qc$mirtrace_pct_tRNA  / 100, threshold = opt$qc_max_trna_frac,  lower = FALSE),
  HIGH_RRNA_FRAC = list(value = qc$mirtrace_pct_rRNA  / 100, threshold = opt$qc_max_rrna_frac,  lower = FALSE),
  LOW_DEPTH      = list(value = qc$mirna_reads_mirtop,       threshold = opt$qc_min_mirna_reads, lower = TRUE)
)
flag_rows <- rbindlist(lapply(names(flag_specs), function(f) {
  s <- flag_specs[[f]]
  hit <- if (s$lower) which(s$value < s$threshold) else which(s$value > s$threshold)
  if (length(hit) == 0) return(NULL)
  data.table(sample_id = as.character(qc$sample_id)[hit], flag = f,
             value = s$value[hit], threshold = s$threshold)
}))
if (is.null(flag_rows)) flag_rows <- data.table(sample_id = character(), flag = character(),
                                                value = numeric(), threshold = numeric())
flagged <- sort(unique(flag_rows$sample_id))

# ------------------------------------------------------------------------
# Composition per sample and per-contrast delta
# ------------------------------------------------------------------------
comp <- data.table(
  sample_id = as.character(qc$sample_id),
  condition = qc$condition,
  norm_group = qc$norm_group
)
for (cat in categories) comp[[cat]] <- qc[[paste0("mirtrace_pct_", cat)]]

comp_long <- melt(comp, id.vars = c("sample_id", "condition", "norm_group"),
                  measure.vars = categories, variable.name = "category", value.name = "pct")

delta_rows <- list()
for (i in seq_len(nrow(con))) {
  r        <- con[i]
  variable <- as.character(r$variable)
  if (!variable %in% names(meta)) {
    warning(sprintf("contrast '%s': variable '%s' not in samplesheet; skipped.", r$contrast_id, variable), call. = FALSE)
    next
  }
  grp <- meta[norm_group == r$norm_group]
  ref <- as.character(grp$sample_id[grp[[variable]] == r$reference])
  trt <- as.character(grp$sample_id[grp[[variable]] == r$treatment])
  if (length(ref) == 0 || length(trt) == 0) {
    warning(sprintf("contrast '%s': empty reference/treatment in '%s'; skipped.", r$contrast_id, r$norm_group), call. = FALSE)
    next
  }
  for (cat in categories) {
    mr <- mean(comp[[cat]][comp$sample_id %in% ref])
    mt <- mean(comp[[cat]][comp$sample_id %in% trt])
    delta_rows[[length(delta_rows) + 1]] <- data.table(
      contrast_id = r$contrast_id, category = cat,
      mean_ref = mr, mean_trt = mt, delta_pp = mt - mr
    )
  }
}
comp_delta <- if (length(delta_rows) > 0) rbindlist(delta_rows) else
  data.table(contrast_id = character(), category = character(),
             mean_ref = numeric(), mean_trt = numeric(), delta_pp = numeric())

# ------------------------------------------------------------------------
# Figures (PDF + PNG)
# ------------------------------------------------------------------------
dir.create("deliverables/qc/plots", recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, stem, width, height) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = "pdf")
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 150)
}

ov <- melt(qc, id.vars = "sample_id",
           measure.vars = c("mirtrace_pct_miRNA", "mirtrace_pct_tRNA", "mirtrace_pct_rRNA",
                            "mirtrace_pct_unknown", "mirna_reads_mirtop"))
ov[, flagged := sample_id %in% flagged]
p_overview <- ggplot(ov, aes(sample_id, value, fill = flagged)) +
  geom_col() +
  facet_wrap(~ variable, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
  labs(x = NULL, y = NULL, fill = "flagged") +
  theme_bw()
save_plot(p_overview, "deliverables/qc/plots/qc_overview", width = 9, height = 10)

comp_long[, sample_id := factor(sample_id, levels = levels(qc$sample_id))]
comp_long[, category := factor(category, levels = categories)]
p_stack <- ggplot(comp_long, aes(sample_id, pct, fill = category)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "percent of reads", fill = NULL) +
  theme_bw()
save_plot(p_stack, "deliverables/qc/plots/composition_stacked_bar", width = 9, height = 10)

# ------------------------------------------------------------------------
# Outputs (published by Nextflow into deliverables/qc/)
# ------------------------------------------------------------------------
fwrite(flag_rows, "deliverables/qc/qc_flags.tsv", sep = "\t", quote = FALSE, na = "")
writeLines(flagged, "deliverables/qc/flagged_samples.txt")

comp_out <- copy(comp)
setcolorder(comp_out, c("sample_id", "condition", "norm_group", categories))
fwrite(comp_out, "deliverables/qc/composition_by_sample.tsv", sep = "\t", quote = FALSE, na = "")
fwrite(comp_delta, "deliverables/qc/composition_delta.tsv", sep = "\t", quote = FALSE, na = "")

# run_info/qa_report.tsv: per-sample pass/flag summary.
dir.create("run_info", recursive = TRUE, showWarnings = FALSE)
flag_by_sample <- split(as.character(flag_rows$flag), as.character(flag_rows$sample_id))
sid <- as.character(qc$sample_id)
qa <- data.frame(
  sample_id = sid,
  condition = as.character(qc$condition),
  norm_group = as.character(qc$norm_group),
  n_flags   = as.integer(vapply(sid, function(s) length(flag_by_sample[[s]]), integer(1))),
  flags     = vapply(sid, function(s) paste(flag_by_sample[[s]], collapse = ","), character(1)),
  pass_qc   = vapply(sid, function(s) length(flag_by_sample[[s]]) == 0, logical(1)),
  stringsAsFactors = FALSE
)
write.table(qa, "run_info/qa_report.tsv", sep = "\t", quote = FALSE, row.names = FALSE, na = "")

message(sprintf("qc_composition ok: %d flagged samples, %d flag rows, %d contrasts.",
                length(flagged), nrow(flag_rows), uniqueN(comp_delta$contrast_id)))
