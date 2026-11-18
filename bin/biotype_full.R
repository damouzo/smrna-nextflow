#!/usr/bin/env Rscript
# Full biotype composition from the mirtrace collapsed FASTA (phase 2).
# Builds annotation once, aligns every sample with bowtie, intersects with bedtools,
# assigns each alignment to the highest-priority overlapped feature, then counts.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

option_list <- list(
  make_option("--nfcore_dir", type = "character", help = "nf-core results dir [required]"),
  make_option("--genome_fasta", type = "character", default = NULL, help = "genome FASTA [optional if index given]"),
  make_option("--bowtie_index", type = "character", default = NULL, help = "bowtie1 index prefix [optional]"),
  make_option("--gtf", type = "character", help = "Ensembl GTF [required]"),
  make_option("--ncrna_bed", type = "character", help = "RNAcentral BED12 [required]"),
  make_option("--samplesheet", type = "character", help = "samplesheet_validated.csv [required]"),
  make_option("--contrasts", type = "character", help = "contrasts_validated.csv [required]"),
  make_option("--trim3p", type = "integer", help = "nt trimmed from 3' end [required]"),
  make_option("--bowtie_mismatches", type = "integer", help = "bowtie -v [required]"),
  make_option("--bowtie_max_hits", type = "integer", help = "max alignments kept [required]"),
  make_option("--min_overlap_frac", type = "double", help = "bedtools -f [required]"),
  make_option("--biotype_priority", type = "character", help = "comma-separated class priority [required]"),
  make_option("--de_classes", type = "character", help = "comma-separated classes with feature-level DE [required]"),
  make_option("--min_cpm", type = "double", help = "min CPM [required]"),
  make_option("--min_total_count", type = "double", help = "min total count [required]"),
  make_option("--alpha", type = "double", help = "FDR threshold [required]"),
  make_option("--lfc_flag", type = "double", help = "sig_lfc1 |log2FC| [required]"),
  make_option("--shrink_type", type = "character", help = "lfcShrink type [required]"),
  make_option("--batch_col", type = "character", default = NULL, help = "samplesheet batch column [optional]"),
  make_option("--threads", type = "integer", help = "bowtie threads [required]"),
  make_option("--seed", type = "integer", help = "random seed [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Full biotype composition."))

req <- c("nfcore_dir", "gtf", "ncrna_bed", "samplesheet", "contrasts", "trim3p",
         "bowtie_mismatches", "bowtie_max_hits", "min_overlap_frac", "biotype_priority",
         "de_classes", "min_cpm", "min_total_count", "alpha", "lfc_flag",
         "shrink_type", "threads", "seed")
for (a in req) if (is.null(opt[[a]])) stop(sprintf("Missing required argument --%s.", a), call. = FALSE)
if (is.null(opt$bowtie_index) && is.null(opt$genome_fasta)) {
  stop("Provide --bowtie_index or --genome_fasta.", call. = FALSE)
}
for (p in c(opt$gtf, opt$ncrna_bed, opt$samplesheet, opt$contrasts)) {
  if (!file.exists(p)) stop(sprintf("File not found: %s", p), call. = FALSE)
}
set.seed(opt$seed)

priority <- trimws(strsplit(opt$biotype_priority, ",", fixed = TRUE)[[1]])
de_classes <- trimws(strsplit(opt$de_classes, ",", fixed = TRUE)[[1]])
meta <- fread(opt$samplesheet)
con  <- fread(opt$contrasts)
samples <- as.character(meta$sample_id)

harmonize_chr <- function(x) {
  x <- sub("^chr", "", x)
  x[x == "M"] <- "MT"
  x
}

# ------------------------------------------------------------------------
# 1. Annotation once: ncRNA from BED12 (drop circRNA/piRNA/lncRNA) + GTF gene lines
# ------------------------------------------------------------------------
build_annotation <- function(path) {
  bed <- fread(opt$ncrna_bed, header = FALSE, sep = "\t", select = c(1, 2, 3, 4, 6, 14))
  setnames(bed, c("chr", "start", "end", "feature_id", "strand", "class"))
  bed <- bed[class %in% priority & !class %in% c("lncRNA", "protein_coding")]
  bed[, chr := harmonize_chr(chr)]
  bed[, `:=`(score = ".")]

  gtf <- fread(opt$gtf, header = FALSE, sep = "\t", quote = "", skip = "##")
  setnames(gtf, c("chr", "src", "feature", "start", "end", "score", "strand", "frame", "attr"))
  gtf <- gtf[feature == "gene"]
  gtf[, class := sub('.*gene_biotype "([^"]+)".*', "\\1", attr)]
  gtf <- gtf[class %in% c("protein_coding", "lncRNA")]
  gtf[, feature_id := sub('.*gene_id "([^"]+)".*', "\\1", attr)]
  gtf[, chr := harmonize_chr(chr)]
  gtf[, start := start - 1L]
  gtf[, score := "."]

  ann <- rbind(bed[, .(chr, start, end, feature_id, score, strand, class)],
               gtf[, .(chr, start, end, feature_id, score, strand, class)])
  ann[, start := as.integer(start)]
  ann[, end := as.integer(end)]
  ann
}
ann <- build_annotation(opt$ncrna_bed)
ann[, prio := match(class, priority)]
fwrite(ann[, .(chr, start, end, feature_id, score, strand, class)], "annotation.bed",
       sep = "\t", col.names = FALSE, quote = FALSE)
message(sprintf("annotation: %d features, %d classes", nrow(ann), uniqueN(ann$class)))

# ------------------------------------------------------------------------
# 2. bowtie index (build once if not provided)
# ------------------------------------------------------------------------
index_prefix <- opt$bowtie_index
if (is.null(index_prefix) || !file.exists(paste0(index_prefix, ".1.ebwt"))) {
  index_prefix <- file.path(getwd(), "bowtie_index")
  fa <- opt$genome_fasta
  if (grepl("\\.gz$", fa)) {
    system2("bash", c("-c", sprintf("zcat '%s' | bowtie-build -f - '%s' >/dev/null 2>bowtie-build.log", fa, index_prefix)), stdout = TRUE)
  } else {
    system2("bowtie-build", c("-f", shQuote(fa), index_prefix),
            stdout = TRUE, stderr = "bowtie-build.log")
  }
  if (!file.exists(paste0(index_prefix, ".1.ebwt"))) stop("bowtie-build failed; see bowtie-build.log", call. = FALSE)
}

# ------------------------------------------------------------------------
# 3. Per sample: align, intersect, assign, count
# ------------------------------------------------------------------------
read_fasta_reads <- function(fa) {
  lines <- readLines(fa, warn = FALSE)
  hdr <- lines[startsWith(lines, ">")]
  name  <- sub("^>", "", sub(" .*$", "", hdr))
  count <- as.numeric(sub(".*_x", "", name))
  data.table(qname = name, count = count)
}

assign_sample <- function(sample) {
  fa <- file.path(opt$nfcore_dir, "mirtrace", sample,
                  "qc_passed_reads.all.collapsed", paste0(sample, ".fastp.fasta"))
  if (!file.exists(fa)) stop(sprintf("Missing collapsed FASTA: %s", fa), call. = FALSE)

  reads <- read_fasta_reads(fa)
  total_reads <- sum(reads$count)

  sam_file <- paste0(sample, ".sam")
  system2("bowtie",
          c("-S", "-f", "-v", opt$bowtie_mismatches, "-k", opt$bowtie_max_hits + 1,
            "--best", "--strata", "-3", opt$trim3p, "-p", opt$threads,
            shQuote(index_prefix), shQuote(fa)),
          stdout = sam_file, stderr = paste0(sample, ".bowtie.log"))

  sam <- fread(sam_file, header = FALSE, skip = "@", select = 1:11)
  if (nrow(sam) == 0) {
    return(list(counts = data.table(class = "unaligned", count = total_reads),
                features = data.table(), stats = data.table(
                  sample = sample, total_reads = total_reads, unaligned = total_reads,
                  multimap_gt_max = 0, annotated = 0, unannotated = 0)))
  }
  setnames(sam, c("qname", "flag", "rname", "pos", "mapq", "cigar",
                  "rnext", "pnext", "tlen", "seq", "qual"),
           skip_absent = TRUE)
  sam <- sam[bitwAnd(flag, 4L) == 0L]
  sam[, read_len := nchar(seq)]
  sam[, strand := fifelse(bitwAnd(flag, 16L) != 0L, "-", "+")]
  sam[, start := pos - 1L]
  sam[, end := pos - 1L + read_len]

  n_by_read <- sam[, .(n_align = .N), by = qname]
  reads <- merge(reads, n_by_read, by = "qname", all.x = TRUE)
  reads[is.na(n_align), n_align := 0L]

  unaligned <- sum(reads$count[reads$n_align == 0L])
  multimap  <- sum(reads$count[reads$n_align > opt$bowtie_max_hits])

  assignable <- reads[n_align >= 1L & n_align <= opt$bowtie_max_hits]
  aln <- sam[qname %in% assignable$qname]
  aln <- merge(aln, assignable[, .(qname, count, n_align)], by = "qname")

  records <- 0; unannotated <- 0; feature_dt <- data.table()
  class_counts <- data.table(class = character(), count = numeric())
  if (nrow(aln) > 0) {
    reads_bed <- paste0(sample, ".reads.bed")
    fwrite(aln[, .(rname, start, end, qname, score = ".", strand)], reads_bed,
           sep = "\t", col.names = FALSE, quote = FALSE)
    inter_file <- paste0(sample, ".intersect.tsv")
    system2("bedtools",
            c("intersect", "-wo", "-s", "-f", opt$min_overlap_frac,
              "-a", shQuote(reads_bed), "-b", "annotation.bed"),
            stdout = inter_file, stderr = paste0(sample, ".bedtools.log"))

    hits <- fread(inter_file, header = FALSE)
    if (nrow(hits) > 0) {
      setnames(hits, c("r_chr", "r_start", "r_end", "r_name", "r_score", "r_strand",
                       "a_chr", "a_start", "a_end", "a_feature", "a_score", "a_strand", "a_class"))
      hits <- merge(hits, ann[, .(feature_id, class, prio)], by.x = "a_feature", by.y = "feature_id")
      setorder(hits, r_chr, r_start, r_end, r_name, prio, a_feature)
      best <- hits[, .SD[1L], by = .(r_chr, r_start, r_end, r_name, r_strand)]
      best <- merge(best, assignable[, .(qname, count, n_align)], by.x = "r_name", by.y = "qname")
      best[, w := count / n_align]

      class_counts <- best[, .(count = sum(w)), by = .(class)]
      feature_dt <- best[class %in% de_classes, .(count = sum(w)), by = .(feature = a_feature, class)]
      annotated_w <- sum(best$w)
      unannotated <- sum(assignable$count) - annotated_w
    } else {
      annotated_w <- 0
      unannotated <- sum(assignable$count)
    }
  }

  counts <- rbindlist(list(
    class_counts[, .(class, count)],
    data.table(class = "unaligned", count = unaligned),
    data.table(class = "multimap_gt_max", count = multimap),
    data.table(class = "unannotated", count = unannotated)
  ), fill = TRUE)

  stats <- data.table(sample = sample, total_reads = total_reads, unaligned = unaligned,
                      multimap_gt_max = multimap,
                      annotated = total_reads - unaligned - multimap - unannotated,
                      unannotated = unannotated)
  file.remove(sam_file)
  list(counts = counts, features = feature_dt, stats = stats)
}

all_counts <- list(); all_features <- list(); all_stats <- list()
for (s in samples) {
  res <- assign_sample(s)
  ct <- res$counts; ct[, sample := s]
  all_counts[[s]] <- ct
  if (nrow(res$features) > 0) { res$features[, sample := s]; all_features[[s]] <- res$features }
  all_stats[[s]] <- res$stats
  message(sprintf("biotype %s done", s))
}

counts_long <- rbindlist(all_counts, fill = TRUE)
counts_long <- counts_long[, .(count = sum(count)), by = .(class, sample)]

class_order <- c(priority, "unannotated", "multimap_gt_max", "unaligned")
counts_wide <- dcast(counts_long, class ~ sample, value.var = "count", fill = 0)
counts_wide[, class := factor(class, levels = intersect(class_order, class))]
setorder(counts_wide, class)
counts_wide[, class := as.character(class)]

# composition per sample (percent of total assigned per sample)
comp <- copy(counts_wide)
tot <- colSums(as.matrix(comp[, -1, with = FALSE]))
comp[, (samples) := lapply(samples, function(s) 100 * comp[[s]] / max(tot[[s]], 1))]

comp_wide <- data.table(sample_id = samples)
for (cl in comp$class) comp_wide[[cl]] <- as.numeric(unlist(comp[class == cl, ..samples]))
setnames(comp_wide, c("sample_id", comp$class))
comp_wide <- merge(data.table(sample_id = samples,
                              condition = as.character(meta$condition),
                              norm_group = as.character(meta$norm_group)),
                   comp_wide, by = "sample_id")

# per-contrast delta on composition
delta_rows <- list()
for (i in seq_len(nrow(con))) {
  r <- con[i]
  variable <- as.character(r$variable)
  grp <- meta[norm_group == r$norm_group]
  ref <- as.character(grp$sample_id[grp[[variable]] == r$reference])
  trt <- as.character(grp$sample_id[grp[[variable]] == r$treatment])
  if (length(ref) == 0 || length(trt) == 0) next
  for (cl in comp$class) {
    mr <- mean(comp_wide[[cl]][comp_wide$sample_id %in% ref])
    mt <- mean(comp_wide[[cl]][comp_wide$sample_id %in% trt])
    delta_rows[[length(delta_rows) + 1]] <- data.table(
      contrast_id = r$contrast_id, category = cl, mean_ref = mr, mean_trt = mt, delta_pp = mt - mr)
  }
}
comp_delta <- if (length(delta_rows) > 0) rbindlist(delta_rows) else
  data.table(contrast_id = character(), category = character(), mean_ref = numeric(),
             mean_trt = numeric(), delta_pp = numeric())

# ------------------------------------------------------------------------
# feature-level DE per de_class
# ------------------------------------------------------------------------
feat_long <- if (length(all_features) > 0) rbindlist(all_features, fill = TRUE) else data.table()

run_de_class <- function(cl) {
  f <- feat_long[class == cl]
  if (nrow(f) == 0) return(invisible(NULL))
  mat <- dcast(f, feature ~ sample, value.var = "count", fill = 0)
  feats <- mat$feature
  m <- as.matrix(mat[, -1, with = FALSE]); rownames(m) <- feats
  for (cid in unique(con$contrast_id)) {
    cr <- con[contrast_id == cid]
    grp <- as.character(cr$norm_group)
    sel <- samples[meta$norm_group == grp & as.character(meta$condition) %in% c(cr$reference, cr$treatment)]
    if (length(unique(as.character(meta$condition[meta$sample_id %in% sel]))) < 2) next
    col <- data.frame(row.names = sel,
                      condition = factor(as.character(meta$condition[match(sel, meta$sample_id)])))
    msel <- m[, sel, drop = FALSE]
    cpm <- sweep(msel, 2, pmax(colSums(msel), 1), "/") * 1e6
    min_samples <- min(table(col$condition))
    keep <- rowSums(cpm >= opt$min_cpm) >= min_samples & rowSums(msel) >= opt$min_total_count
    if (sum(keep) < 2) next
    dds <- DESeqDataSetFromMatrix(msel[keep, , drop = FALSE], col, design = ~ condition)
    dds <- tryCatch(DESeq(dds, fitType = "parametric", quiet = TRUE),
                    error = function(e) DESeq(dds, fitType = "mean", quiet = TRUE))
    res <- results(dds, contrast = c("condition", as.character(cr$treatment), as.character(cr$reference)))
    shr <- lfcShrink(dds, contrast = c("condition", as.character(cr$treatment), as.character(cr$reference)),
                     type = opt$shrink_type, quiet = TRUE)
    dt <- data.table(feature = rownames(res), baseMean = res$baseMean,
                     log2FC_mle = res$log2FoldChange, log2FC_shrunk = shr$log2FoldChange,
                     lfcSE = res$lfcSE, stat = res$stat, pvalue = res$pvalue, padj = res$padj)
    dt[, direction := fifelse(!is.na(padj) & padj < opt$alpha, fifelse(log2FC_shrunk > 0, "up", "down"), "ns")]
    dt[, sig_lfc1 := !is.na(padj) & padj < opt$alpha & abs(log2FC_shrunk) >= opt$lfc_flag]
    dir.create(file.path("deliverables/biotype/de", cl, cid), recursive = TRUE, showWarnings = FALSE)
    fwrite(dt, file.path("deliverables/biotype/de", cl, cid, "de_full.tsv"),
           sep = "\t", quote = FALSE, na = "")
    fwrite(dt[!is.na(padj) & padj < opt$alpha],
           file.path("deliverables/biotype/de", cl, cid, "de_significant.tsv"),
           sep = "\t", quote = FALSE, na = "")
  }
}
for (cl in de_classes) run_de_class(cl)

# ------------------------------------------------------------------------
# Outputs
# ------------------------------------------------------------------------
dir.create("deliverables/biotype", recursive = TRUE, showWarnings = FALSE)
dir.create("raw/biotype/feature_counts", recursive = TRUE, showWarnings = FALSE)

fwrite(counts_wide, "deliverables/biotype/biotype_counts.tsv", sep = "\t", quote = FALSE)
fwrite(comp_wide,  "deliverables/biotype/biotype_composition.tsv", sep = "\t", quote = FALSE, na = "")
fwrite(comp_delta, "deliverables/biotype/composition_delta.tsv", sep = "\t", quote = FALSE, na = "")
fwrite(rbindlist(all_stats, fill = TRUE), "raw/biotype/alignment_stats.tsv", sep = "\t", quote = FALSE)

if (nrow(feat_long) > 0) {
  for (cl in unique(feat_long$class)) {
    w <- dcast(feat_long[class == cl], feature ~ sample, value.var = "count", fill = 0)
    fwrite(w, file.path("raw/biotype/feature_counts", paste0(cl, ".tsv")), sep = "\t", quote = FALSE)
  }
}

comp_long <- melt(comp_wide, id.vars = c("sample_id", "condition", "norm_group"),
                  variable.name = "category", value.name = "pct")
p <- ggplot(comp_long, aes(factor(sample_id, levels = samples), pct, fill = category)) +
  geom_col() + coord_flip() + labs(x = NULL, y = "percent", fill = NULL) + theme_bw()
ggsave("deliverables/biotype/composition_stacked_bar.pdf", p, width = 9, height = 10, device = "pdf")
ggsave("deliverables/biotype/composition_stacked_bar.png", p, width = 9, height = 10, dpi = 150)

message(sprintf("biotype_full ok: %d classes, %d samples.", nrow(counts_wide), length(samples)))
