#!/usr/bin/env Rscript
# isomiR profile (phase 3). Requires reads with a clean 3' end.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

option_list <- list(
  make_option("--isomir_filtered", type = "character", help = "raw/counts/isomir_filtered.tsv.gz [required]"),
  make_option("--samplesheet", type = "character", help = "samplesheet_validated.csv [required]"),
  make_option("--contrasts", type = "character", help = "contrasts_validated.csv [required]"),
  make_option("--isomir_min_reads", type = "double", help = "min weighted reads per miRNA/sample [required]"),
  make_option("--min_cpm", type = "double", help = "min CPM for isomiR DE [required]"),
  make_option("--min_total_count", type = "double", help = "min total count for isomiR DE [required]"),
  make_option("--alpha", type = "double", help = "FDR threshold [required]"),
  make_option("--lfc_flag", type = "double", help = "sig_lfc1 |log2FC| [required]"),
  make_option("--shrink_type", type = "character", help = "lfcShrink type [required]"),
  make_option("--seed", type = "integer", help = "random seed [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "isomiR profile."))

for (a in c("isomir_filtered", "samplesheet", "contrasts", "isomir_min_reads",
            "min_cpm", "min_total_count", "alpha", "lfc_flag", "shrink_type", "seed")) {
  if (is.null(opt[[a]])) stop(sprintf("Missing required argument --%s.", a), call. = FALSE)
}
for (p in c(opt$isomir_filtered, opt$samplesheet, opt$contrasts)) {
  if (!file.exists(p)) stop(sprintf("File not found: %s", p), call. = FALSE)
}
set.seed(opt$seed)

meta <- fread(opt$samplesheet)
con  <- fread(opt$contrasts)
samples <- as.character(meta$sample_id)

dt <- fread(opt$isomir_filtered)
sample_cols <- intersect(samples, names(dt))
if (length(sample_cols) != length(samples)) stop("isomir_filtered.tsv.gz is missing some samplesheet samples.", call. = FALSE)

m <- as.matrix(dt[, ..sample_cols]) * dt$weight        # weighted counts (1/k)
miRNA <- dt$miRNA

iso5 <- as.integer(dt$iso_5p)
iso3 <- as.integer(dt$iso_3p)
add3 <- as.integer(dt$iso_add3p)
snp  <- as.integer(dt$iso_snp)

tot <- rowsum(m, group = miRNA)                         # miRNA x sample
sum_if <- function(mask) {
  if (sum(mask) == 0) return(matrix(0, nrow = nrow(tot), ncol = ncol(tot),
                                    dimnames = dimnames(tot)))
  rr <- rowsum(m[mask, , drop = FALSE], group = miRNA[mask])
  rr[rownames(tot), , drop = FALSE]
}
canon <- sum_if(iso5 == 0 & iso3 == 0 & add3 == 0 & snp == 0)
s5    <- sum_if(iso5 != 0)
s3    <- sum_if(iso3 != 0)
a3    <- sum_if(add3 != 0)

profile <- rbindlist(lapply(seq_along(samples), function(j) {
  smp <- samples[j]
  keep <- tot[, j] >= opt$isomir_min_reads
  if (!any(keep)) return(NULL)
  mirs <- rownames(tot)[keep]
  data.table(sample_id = smp, mirna = mirs,
             total_reads = tot[mirs, j],
             frac_canonical   = canon[mirs, j] / tot[mirs, j],
             frac_5p_shifted  = s5[mirs, j] / tot[mirs, j],
             frac_3p_shifted  = s3[mirs, j] / tot[mirs, j],
             frac_3p_added    = a3[mirs, j] / tot[mirs, j])
}))

# dominant isoform per miRNA/sample
isomir_id <- paste(miRNA, iso5, iso3, add3, sep = ":")
agg <- rowsum(m, group = isomir_id)
agg_dt <- data.table(isomir = rownames(agg), as.data.table(agg))
agg_dt[, mirna := sub(":.*", "", isomir)]
long_agg <- melt(agg_dt, id.vars = c("isomir", "mirna"), variable.name = "sample", value.name = "w")
dominant <- long_agg[, .(dominant_isomir = isomir[which.max(w)]), by = .(sample, mirna)]
profile <- merge(profile, dominant, by.x = c("sample_id", "mirna"), by.y = c("sample", "mirna"), all.x = TRUE)
setorder(profile, sample_id, mirna)

# ------------------------------------------------------------------------
# delta 5' shift per contrast (Welch t on logit, indicative)
# ------------------------------------------------------------------------
logit <- function(p) log((p + 1e-3) / (1 - p + 1e-3))
wide5 <- dcast(profile, mirna ~ sample_id, value.var = "frac_5p_shifted")

delta_rows <- list()
for (i in seq_len(nrow(con))) {
  r <- con[i]
  grp <- meta[norm_group == r$norm_group]
  ref <- intersect(as.character(grp$sample_id[grp[[as.character(r$variable)]] == r$reference]), colnames(wide5))
  trt <- intersect(as.character(grp$sample_id[grp[[as.character(r$variable)]] == r$treatment]), colnames(wide5))
  if (length(ref) < 2 || length(trt) < 2) next
  rows <- list()
  for (mir in wide5$mirna) {
    a <- logit(as.numeric(wide5[mirna == mir, ..ref]))
    b <- logit(as.numeric(wide5[mirna == mir, ..trt]))
    a <- a[is.finite(a)]; b <- b[is.finite(b)]
    if (length(a) < 2 || length(b) < 2) next
    tt <- tryCatch(t.test(a, b), error = function(e) NULL)
    rows[[length(rows) + 1]] <- data.table(
      contrast_id = r$contrast_id, mirna = mir,
      n_ref = length(a), n_trt = length(b),
      mean_ref = mean(a), mean_trt = mean(b), delta = mean(b) - mean(a),
      t_stat = if (is.null(tt)) NA_real_ else unname(tt$statistic),
      p_value = if (is.null(tt)) NA_real_ else tt$p.value)
  }
  if (length(rows) > 0) {
    d <- rbindlist(rows)
    dir.create(file.path("deliverables/isomir"), recursive = TRUE, showWarnings = FALSE)
    fwrite(d, file.path("deliverables/isomir", paste0("delta_5p_shift_", r$contrast_id, ".tsv")),
           sep = "\t", quote = FALSE, na = "")
    delta_rows[[length(delta_rows) + 1]] <- d
  }
}

# ------------------------------------------------------------------------
# isomiR DE per contrast (feature = isomir)
# ------------------------------------------------------------------------
counts_iso <- round(agg)
for (i in seq_len(nrow(con))) {
  r <- con[i]
  grp <- as.character(r$norm_group)
  sel <- samples[meta$norm_group == grp & as.character(meta$condition) %in%
                   c(as.character(r$reference), as.character(r$treatment))]
  if (length(unique(as.character(meta$condition[meta$sample_id %in% sel]))) < 2) next
  col <- data.frame(row.names = sel,
                    condition = factor(as.character(meta$condition[match(sel, meta$sample_id)])))
  msel <- counts_iso[, sel, drop = FALSE]
  cpm <- sweep(msel, 2, pmax(colSums(msel), 1), "/") * 1e6
  min_samples <- min(table(col$condition))
  keep <- rowSums(cpm >= opt$min_cpm) >= min_samples & rowSums(msel) >= opt$min_total_count
  if (sum(keep) < 2) next
  dds <- DESeqDataSetFromMatrix(msel[keep, , drop = FALSE], col, design = ~ condition)
  dds <- tryCatch(DESeq(dds, fitType = "parametric", quiet = TRUE),
                  error = function(e) DESeq(dds, fitType = "mean", quiet = TRUE))
  res <- results(dds, contrast = c("condition", as.character(r$treatment), as.character(r$reference)))
  shr <- lfcShrink(dds, contrast = c("condition", as.character(r$treatment), as.character(r$reference)),
                   type = opt$shrink_type, quiet = TRUE)
  d <- data.table(isomir = rownames(res), baseMean = res$baseMean,
                  log2FC_mle = res$log2FoldChange, log2FC_shrunk = shr$log2FoldChange,
                  lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
  d[, direction := fifelse(!is.na(padj) & padj < opt$alpha, fifelse(log2FC_shrunk > 0, "up", "down"), "ns")]
  d[, sig_lfc1 := !is.na(padj) & padj < opt$alpha & abs(log2FC_shrunk) >= opt$lfc_flag]
  od <- file.path("deliverables/isomir/de_isomir", as.character(r$contrast_id))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  fwrite(d, file.path(od, "de_full.tsv"), sep = "\t", quote = FALSE, na = "")
  fwrite(d[!is.na(padj) & padj < opt$alpha], file.path(od, "de_significant.tsv"),
         sep = "\t", quote = FALSE, na = "")
  sig <- d[!is.na(padj) & padj < opt$alpha & abs(log2FC_shrunk) >= opt$lfc_flag]
  lab <- sig[order(padj)][seq_len(min(15, nrow(sig)))]
  v <- ggplot(d, aes(log2FC_shrunk, -log10(padj), color = direction)) +
    geom_point(alpha = 0.6) + ggrepel::geom_text_repel(data = lab, aes(label = isomir), size = 2.5) +
    labs(x = "log2FC (shrunken)", y = "-log10(padj)", color = NULL) + theme_bw()
  ggsave(file.path(od, "volcano.pdf"), v, width = 7, height = 6, device = "pdf")
  ggsave(file.path(od, "volcano.png"), v, width = 7, height = 6, dpi = 150)
}

# ------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------
dir.create("deliverables/isomir", recursive = TRUE, showWarnings = FALSE)
fwrite(profile, "deliverables/isomir/isomir_profile_by_sample.tsv", sep = "\t", quote = FALSE, na = "")

iso_counts_out <- data.table(isomir = rownames(counts_iso), as.data.table(counts_iso))
setnames(iso_counts_out, c("isomir", samples))
fwrite(iso_counts_out, "deliverables/isomir/isomir_counts.tsv", sep = "\t", quote = FALSE)

# overall fraction per sample (weighted across miRNAs)
overall <- data.table(
  sample_id = samples,
  canonical  = colSums(canon) / colSums(tot),
  shift_5p   = colSums(s5) / colSums(tot),
  shift_3p   = colSums(s3) / colSums(tot),
  added_3p   = colSums(a3) / colSums(tot)
)
ol <- melt(overall, id.vars = "sample_id", variable.name = "category", value.name = "frac")
ol[, sample_id := factor(sample_id, levels = samples)]
p <- ggplot(ol, aes(sample_id, frac, fill = category)) + geom_col() + coord_flip() +
  labs(x = NULL, y = "fraction of reads", fill = NULL) + theme_bw()
ggsave("deliverables/isomir/isomir_fraction_stacked_bar.pdf", p, width = 9, height = 10, device = "pdf")
ggsave("deliverables/isomir/isomir_fraction_stacked_bar.png", p, width = 9, height = 10, dpi = 150)

message(sprintf("isomir_profile ok: %d isomiRs, %d profile rows.", nrow(counts_iso), nrow(profile)))
