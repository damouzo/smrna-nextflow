#!/usr/bin/env Rscript
# DESeq2 per norm_group. One invocation handles one group (or "all") and
# every run (main + sensitivities). main writes deliverables; sensitivities raw.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

option_list <- list(
  make_option("--counts_mirtop", type = "character", help = "counts_mirtop.tsv [required]"),
  make_option("--counts_mature", type = "character", help = "counts_mature.tsv [required]"),
  make_option("--samplesheet", type = "character", help = "samplesheet_validated.csv [required]"),
  make_option("--contrasts", type = "character", help = "contrasts_validated.csv [required]"),
  make_option("--flagged", type = "character", help = "flagged_samples.txt [required]"),
  make_option("--group", type = "character", help = "norm_group value, or 'all' [required]"),
  make_option("--count_source", type = "character", help = "mirtop or mature [required]"),
  make_option("--min_cpm", type = "double", help = "min CPM [required]"),
  make_option("--min_total_count", type = "double", help = "min total count [required]"),
  make_option("--alpha", type = "double", help = "FDR threshold [required]"),
  make_option("--lfc_flag", type = "double", help = "sig_lfc1 |log2FC| threshold [required]"),
  make_option("--shrink_type", type = "character", help = "lfcShrink type [required]"),
  make_option("--n_top_var", type = "integer", help = "miRNAs used for PCA [required]"),
  make_option("--batch_col", type = "character", default = NULL, help = "samplesheet batch column (e.g. batch) [optional]"),
  make_option("--seed", type = "integer", help = "random seed [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "DESeq2 per norm_group."))

for (a in c("counts_mirtop", "counts_mature", "samplesheet", "contrasts", "flagged",
            "group", "count_source", "min_cpm", "min_total_count", "alpha",
            "lfc_flag", "shrink_type", "n_top_var", "seed")) {
  if (is.null(opt[[a]])) stop(sprintf("Missing required argument --%s. Run with --help.", a), call. = FALSE)
}
for (p in c(opt$counts_mirtop, opt$counts_mature, opt$samplesheet, opt$contrasts)) {
  if (!file.exists(p)) stop(sprintf("File not found: %s", p), call. = FALSE)
}
set.seed(opt$seed)

group <- opt$group
meta  <- fread(opt$samplesheet)
con   <- fread(opt$contrasts)
flagged <- if (file.exists(opt$flagged)) readLines(opt$flagged, warn = FALSE) else character(0)

read_counts <- function(path) {
  m <- fread(path)
  mat <- as.matrix(m[, -1, with = FALSE])
  rownames(mat) <- m[[1]]
  mat
}
counts_mirtop <- read_counts(opt$counts_mirtop)
counts_mature <- read_counts(opt$counts_mature)
main_counts   <- if (identical(opt$count_source, "mature")) counts_mature else counts_mirtop

grp_samples <- if (identical(group, "all")) as.character(meta$sample_id) else
  as.character(meta$sample_id[meta$norm_group == group])
if (length(grp_samples) == 0) stop(sprintf("norm_group '%s' has no samples.", group), call. = FALSE)

grp_con <- if (identical(group, "all")) con[0] else con[norm_group == group]
if (nrow(grp_con) > 0 && any(as.character(grp_con$variable) != "condition")) {
  warning("Some contrasts use a variable other than 'condition'; condition is always in the model.", call. = FALSE)
}

meta_idx <- setNames(seq_len(nrow(meta)), as.character(meta$sample_id))

# ------------------------------------------------------------------------
# Design: per group, ~ <batch_col> + condition when the batch column has >= 2
# non-empty levels and the design matrix is full rank; else ~ condition.
# Donor groups with an empty batch column fall back to ~ condition and do not
# block the cell-line groups.
# ------------------------------------------------------------------------
batch_col <- if (!is.null(opt$batch_col) && nzchar(opt$batch_col)) opt$batch_col else NULL
batch_usable <- FALSE
if (!is.null(batch_col)) {
  if (!batch_col %in% names(meta)) {
    message(sprintf("batch_col '%s' is not a samplesheet column; using ~ condition.", batch_col))
  } else {
    v <- as.character(meta[[batch_col]][meta_idx[grp_samples]])
    v[is.na(v)] <- ""
    if (all(nzchar(v)) && length(unique(v)) >= 2) {
      col_tmp <- data.frame(batch = factor(v),
                            condition = factor(as.character(meta$condition[meta_idx[grp_samples]])))
      mm <- model.matrix(~ batch + condition, data = col_tmp)
      batch_usable <- qr(mm)$rank == ncol(mm)
      if (!batch_usable) {
        message(sprintf("group '%s': batch '%s' is confounded with condition; using ~ condition.", group, batch_col))
      }
    }
  }
}
main_design <- if (batch_usable) sprintf("~ %s + condition", batch_col) else "~ condition"
message(sprintf("group '%s': main design %s", group, main_design))

# ------------------------------------------------------------------------
# DESeq2 helpers
# ------------------------------------------------------------------------
build_dds <- function(counts, sel, design_str) {
  sel <- sel[sel %in% colnames(counts)]
  col <- data.frame(row.names = sel,
                    condition = factor(as.character(meta$condition[meta_idx[sel]])))
  if (!is.null(batch_col) && batch_col %in% names(meta)) {
    col[[batch_col]] <- factor(as.character(meta[[batch_col]][meta_idx[sel]]))
  }
  mat <- counts[, sel, drop = FALSE]
  cpm <- sweep(mat, 2, pmax(colSums(mat), 1), "/") * 1e6
  min_samples <- min(table(col$condition))
  keep <- rowSums(cpm >= opt$min_cpm) >= min_samples & rowSums(mat) >= opt$min_total_count
  dds <- DESeqDataSetFromMatrix(mat[keep, , drop = FALSE], col, design = as.formula(design_str))
  dds <- tryCatch(estimateSizeFactors(dds, type = "ratio"),
                  error = function(e) {
                    message("size factors: ratio failed, using poscounts")
                    estimateSizeFactors(dds, type = "poscounts")
                  })
  dds
}

fit_dds <- function(dds) {
  tryCatch(DESeq(dds, fitType = "parametric", quiet = TRUE),
           error = function(e) {
             message("fitType parametric failed, trying local")
             tryCatch(DESeq(dds, fitType = "local", quiet = TRUE),
                      error = function(e2) {
                        message("fitType local failed, using mean")
                        DESeq(dds, fitType = "mean", quiet = TRUE)
                      })
           })
}

de_table <- function(dds, contrast_id, reference, treatment) {
  res <- results(dds, contrast = c("condition", treatment, reference))
  shr <- lfcShrink(dds, contrast = c("condition", treatment, reference),
                   type = opt$shrink_type, quiet = TRUE)
  norm <- counts(dds, normalized = TRUE)
  is_ref <- colData(dds)$condition == reference
  is_trt <- colData(dds)$condition == treatment
  dt <- data.table(
    mirna          = rownames(res),
    baseMean       = res$baseMean,
    log2FC_mle     = res$log2FoldChange,
    log2FC_shrunk  = shr$log2FoldChange,
    lfcSE          = res$lfcSE,
    stat           = res$stat,
    pvalue         = res$pvalue,
    padj           = res$padj
  )
  dt[, direction := fifelse(!is.na(padj) & padj < opt$alpha,
                            fifelse(log2FC_shrunk > 0, "up", "down"), "ns")]
  dt[, sig_lfc1 := !is.na(padj) & padj < opt$alpha & abs(log2FC_shrunk) >= opt$lfc_flag]
  dt[, n_ref := sum(is_ref)]
  dt[, n_trt := sum(is_trt)]
  dt[, mean_norm_count_ref := rowMeans(norm[, is_ref, drop = FALSE])]
  dt[, mean_norm_count_trt := rowMeans(norm[, is_trt, drop = FALSE])]
  dt[]
}

write_de <- function(dt, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(dt, file.path(dir, "de_full.tsv"), sep = "\t", quote = FALSE, na = "")
  fwrite(dt[!is.na(padj) & padj < opt$alpha], file.path(dir, "de_significant.tsv"),
         sep = "\t", quote = FALSE, na = "")
}

# ------------------------------------------------------------------------
# Exploration (main only): VST, PCA, sample distances
# ------------------------------------------------------------------------
explore <- function(dds, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # vst() subsamples nsub=1000 rows and errors when there are fewer; fall back
  # to the full transformation for small designs.
  vsd <- if (nrow(dds) < 1000) varianceStabilizingTransformation(dds, blind = TRUE) else vst(dds, blind = TRUE)
  v <- assay(vsd)
  ntop <- min(opt$n_top_var, nrow(v))
  v <- v[order(apply(v, 1, stats::var), decreasing = TRUE)[seq_len(ntop)], , drop = FALSE]
  pca <- prcomp(t(v))

  coords <- data.table(sample_id = rownames(pca$x),
                       PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                       condition = colData(dds)$condition,
                       batch = if (!is.null(batch_col) && batch_col %in% names(colData(dds)))
                                 as.character(colData(dds)[[batch_col]]) else NA_character_)
  coords[, flagged := sample_id %in% flagged]
  fwrite(coords, file.path(out_dir, "pca_coords.tsv"), sep = "\t", quote = FALSE, na = "")

  # Shape by batch only when every sample has one: mapping shape to an
  # all-NA column drops every point (the donor groups have no batch).
  pca_aes <- if (!anyNA(coords$batch)) {
    aes(PC1, PC2, color = condition, shape = batch)
  } else {
    aes(PC1, PC2, color = condition)
  }

  pca_plot <- ggplot(coords, pca_aes) +
    geom_point(size = 3, stroke = 0.5) +
    geom_point(data = coords[flagged == TRUE], color = "red", size = 3, shape = 1) +
    ggrepel::geom_text_repel(aes(label = sample_id), size = 2.5, show.legend = FALSE) +
    labs(color = "condition", shape = "batch") +
    theme_bw()
  ggsave(file.path(out_dir, "pca.pdf"), pca_plot, width = 7, height = 6, device = "pdf")
  ggsave(file.path(out_dir, "pca.png"), pca_plot, width = 7, height = 6, dpi = 150)

  samp_dists <- dist(t(assay(vsd)))
  m <- as.matrix(samp_dists)
  pheatmap(m, filename = file.path(out_dir, "sample_distance_heatmap.pdf"))
  grDevices::png(file.path(out_dir, "sample_distance_heatmap.png"), width = 1400, height = 1200, res = 150)
  pheatmap(m)
  grDevices::dev.off()
}

# ------------------------------------------------------------------------
# Runs
# ------------------------------------------------------------------------
runs <- list(main = list(counts = main_counts, samples = grp_samples, design = main_design))
if (!identical(group, "all")) {
  if (!identical(opt$count_source, "mature")) {
    runs$sens_mature <- list(counts = counts_mature, samples = grp_samples, design = main_design)
  }
  nf_samples <- setdiff(grp_samples, flagged)
  if (length(nf_samples) < length(grp_samples)) {
    runs$sens_no_flagged <- list(counts = main_counts, samples = nf_samples, design = main_design)
  }
  if (batch_usable) {
    runs$sens_condition <- list(counts = main_counts, samples = grp_samples, design = "~ condition")
  }
}

for (run_name in names(runs)) {
  r <- runs[[run_name]]
  dds <- build_dds(r$counts, r$samples, r$design)

  if (identical(run_name, "main")) {
    dir.create(file.path("raw/deseq2", group), recursive = TRUE, showWarnings = FALSE)
    sf <- data.table(sample_id = colnames(dds), size_factor = sizeFactors(dds))
    fwrite(sf, file.path("raw/deseq2", group, "size_factors.tsv"), sep = "\t", quote = FALSE)
    fc <- data.table(mirna = rownames(counts(dds)), as.data.table(counts(dds)))
    fwrite(fc, file.path("raw/deseq2", group, "filtered_counts.tsv"), sep = "\t", quote = FALSE)
  }

  # All group: exploration only, no DE
  if (identical(group, "all")) {
    if (identical(run_name, "main")) {
      dds <- fit_dds(dds)
      explore(dds, file.path("deliverables/explore", group))
      saveRDS(dds, file.path("raw/deseq2", group, "dds.rds"))
    }
    next
  }

  dds <- fit_dds(dds)
  if (identical(run_name, "main")) {
    saveRDS(dds, file.path("raw/deseq2", group, "dds.rds"))
    explore(dds, file.path("deliverables/explore", group))
  }

  for (i in seq_len(nrow(grp_con))) {
    cid <- as.character(grp_con$contrast_id[i])
    ref <- as.character(grp_con$reference[i])
    trt <- as.character(grp_con$treatment[i])
    if (sum(colData(dds)$condition == ref) < 2 || sum(colData(dds)$condition == trt) < 2) {
      message(sprintf("%s / %s: fewer than 2 samples per level, skipped.", run_name, cid))
      next
    }
    dt <- de_table(dds, cid, ref, trt)
    if (identical(run_name, "main")) {
      write_de(dt, file.path("deliverables/de", cid))
      dir.create(file.path("deliverables/de", cid, "plots"), recursive = TRUE, showWarnings = FALSE)
      sig <- dt[!is.na(padj) & padj < opt$alpha & abs(log2FC_shrunk) >= opt$lfc_flag]
      lab <- sig[order(padj)][seq_len(min(15, nrow(sig)))]
      volcano <- ggplot(dt, aes(log2FC_shrunk, -log10(padj), color = direction)) +
        geom_point(alpha = 0.6) +
        geom_text_repel(data = lab, aes(label = mirna), size = 2.5, show.legend = FALSE) +
        labs(x = "log2FC (shrunken)", y = "-log10(padj)", color = NULL) +
        theme_bw()
      ggsave(file.path("deliverables/de", cid, "plots", "volcano.pdf"), volcano, width = 7, height = 6, device = "pdf")
      ggsave(file.path("deliverables/de", cid, "plots", "volcano.png"), volcano, width = 7, height = 6, dpi = 150)
    } else {
      write_de(dt, file.path("raw/deseq2", run_name, cid))
    }
  }
}

# Completion stamp for downstream ordering
dir.create("raw/deseq2", recursive = TRUE, showWarnings = FALSE)
writeLines(group, "deseq_group.done")
message(sprintf("deseq2_group ok: group '%s', runs: %s.", group, paste(names(runs), collapse = ", ")))
