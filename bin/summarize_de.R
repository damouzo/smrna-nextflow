#!/usr/bin/env Rscript
# Summarize DE runs and compare main against each sensitivity.
# Reads the published outdir tree (deliverables/de and raw/deseq2).

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--outdir", type = "character", help = "pipeline outdir [required]"),
  make_option("--alpha", type = "double", help = "FDR threshold [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Summarize DE runs and concordance."))

for (a in c("outdir", "alpha")) {
  if (is.null(opt[[a]])) stop(sprintf("Missing required argument --%s. Run with --help.", a), call. = FALSE)
}
if (!dir.exists(opt$outdir)) stop(sprintf("--outdir does not exist: %s", opt$outdir), call. = FALSE)

alpha     <- opt$alpha
main_root <- file.path(opt$outdir, "deliverables/de")
sens_root <- file.path(opt$outdir, "raw/deseq2")

de_path <- function(run, cid) {
  if (identical(run, "main")) file.path(main_root, cid, "de_full.tsv")
  else file.path(sens_root, run, cid, "de_full.tsv")
}
contrast_ids <- function(dir) {
  if (!dir.exists(dir)) return(character(0))
  sub <- list.dirs(dir, recursive = FALSE)
  sub[file.exists(file.path(sub, "de_full.tsv"))]
}

# main contrasts + sensitivity runs
main_dirs   <- contrast_ids(main_root)
main_cids   <- basename(main_dirs)
sens_runs   <- list.dirs(sens_root, recursive = FALSE)
sens_runs   <- basename(sens_runs[startsWith(basename(sens_runs), "sens_")])
sens_runs   <- sens_runs[order(sens_runs)]

read_de <- function(run, cid) {
  p <- de_path(run, cid)
  if (file.exists(p)) fread(p) else NULL
}

# ------------------------------------------------------------------------
# de_summary.tsv
# ------------------------------------------------------------------------
summarize <- function(dt, run, cid) {
  data.table(
    run         = run,
    contrast_id = cid,
    n_tested    = sum(!is.na(dt$pvalue)),
    n_up        = sum(dt$direction == "up", na.rm = TRUE),
    n_down      = sum(dt$direction == "down", na.rm = TRUE),
    n_sig_lfc1  = sum(dt$sig_lfc1, na.rm = TRUE)
  )
}

sens_cids <- unlist(lapply(sens_runs, function(r) basename(contrast_ids(file.path(sens_root, r)))))
all_cids  <- sort(unique(c(main_cids, sens_cids)))

summary_rows <- list()
for (cid in all_cids) {
  for (run in c("main", sens_runs)) {
    dt <- read_de(run, cid)
    if (!is.null(dt) && nrow(dt) > 0) summary_rows[[length(summary_rows) + 1]] <- summarize(dt, run, cid)
  }
}
de_summary <- if (length(summary_rows) > 0) rbindlist(summary_rows) else
  data.table(run = character(), contrast_id = character(), n_tested = integer(),
             n_up = integer(), n_down = integer(), n_sig_lfc1 = integer())

# ------------------------------------------------------------------------
# concordance.tsv — main vs each sensitivity, per contrast
# ------------------------------------------------------------------------
conc_rows <- list()
for (cid in main_cids) {
  m <- read_de("main", cid)
  if (is.null(m)) next
  for (run in sens_runs) {
    s <- read_de(run, cid)
    if (is.null(s)) next
    j <- merge(m[, .(mirna, lfc_main = log2FC_shrunk, padj_main = padj)],
               s[, .(mirna, lfc_sens = log2FC_shrunk, padj_sens = padj)],
               by = "mirna")
    ok <- !is.na(j$lfc_main) & !is.na(j$lfc_sens)
    rho <- if (sum(ok) >= 3) stats::cor(j$lfc_main[ok], j$lfc_sens[ok], method = "spearman") else NA_real_
    conc_rows[[length(conc_rows) + 1]] <- data.table(
      contrast_id    = cid,
      run_a          = "main",
      run_b          = run,
      n_common       = sum(ok),
      spearman_log2fc = rho,
      n_sig_shared   = sum(!is.na(j$padj_main) & j$padj_main < alpha &
                             !is.na(j$padj_sens) & j$padj_sens < alpha)
    )
  }
}
concordance <- if (length(conc_rows) > 0) rbindlist(conc_rows) else
  data.table(contrast_id = character(), run_a = character(), run_b = character(),
             n_common = integer(), spearman_log2fc = numeric(), n_sig_shared = integer())

# ------------------------------------------------------------------------
# Outputs (published to deliverables/de/)
# ------------------------------------------------------------------------
dir.create(file.path("deliverables/de"), recursive = TRUE, showWarnings = FALSE)
fwrite(de_summary,  "deliverables/de/de_summary.tsv",  sep = "\t", quote = FALSE, na = "")
fwrite(concordance, "deliverables/de/concordance.tsv", sep = "\t", quote = FALSE, na = "")

message(sprintf("summarize_de ok: %d summary rows, %d concordance rows, %d sensitivity runs.",
                nrow(de_summary), nrow(concordance), length(sens_runs)))
