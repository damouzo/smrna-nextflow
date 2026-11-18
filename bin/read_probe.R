#!/usr/bin/env Rscript
# Count reads derived from given probe sequences in the mirtrace collapsed
# FASTA. Generic: works for shRNA constructs or any sequence of interest. A read
# counts if it equals a substring of a probe (>= probe_min_len nt) in either
# orientation; optional 1-mismatch matching. Also dumps the top 'unknown'
# sequences per sample so unannotated constructs can be spotted without a probe.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option("--nfcore_dir", type = "character",
              help = "nf-core/smrnaseq results directory [required]"),
  make_option("--samplesheet", type = "character",
              help = "samplesheet_validated.csv [required]"),
  make_option("--probe_file", type = "character", default = NULL,
              help = "probes.tsv with columns probe_id, sequence [optional]"),
  make_option("--probe_mismatches", type = "integer", default = 0,
              help = "0 or 1 mismatches allowed [default 0]"),
  make_option("--probe_min_len", type = "integer", default = 18,
              help = "minimum match length [default 18]"),
  make_option("--top_unknown_n", type = "integer", default = 0,
              help = "top N 'unknown' sequences per sample, 0 disables [default 0]"),
  make_option("--top_hits_n", type = "integer", default = 20,
              help = "top N reads per sample and probe in probe_hits.tsv [default 20]"),
  make_option("--seed", type = "integer", help = "random seed [required]")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Count reads derived from probe sequences."))

for (a in c("nfcore_dir", "samplesheet", "seed")) {
  if (is.null(opt[[a]])) stop(sprintf("Missing required argument --%s. Run with --help.", a), call. = FALSE)
}
if (!dir.exists(opt$nfcore_dir)) stop(sprintf("--nfcore_dir does not exist: %s", opt$nfcore_dir), call. = FALSE)
if (!file.exists(opt$samplesheet)) stop(sprintf("--samplesheet does not exist: %s", opt$samplesheet), call. = FALSE)
if (!opt$probe_mismatches %in% c(0L, 1L)) stop("--probe_mismatches must be 0 or 1.", call. = FALSE)
if (opt$probe_min_len < 1L) stop("--probe_min_len must be >= 1.", call. = FALSE)
set.seed(opt$seed)

samples <- as.character(fread(opt$samplesheet)$sample_id)

# Uppercase, drop anything that is not a letter, read U as T (RNA input).
norm_seq <- function(x) gsub("U", "T", gsub("[^A-Za-z]", "", toupper(x)))

read_probes <- function(path) {
  empty <- data.table(probe_id = character(), sequence = character())
  if (is.null(path) || !nzchar(path) || !file.exists(path)) return(empty)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!startsWith(trimws(lines), "#") & nzchar(trimws(lines))]
  if (length(lines) < 1) return(empty)
  hdr <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  if (!all(c("probe_id", "sequence") %in% hdr)) {
    stop("probe_file must have columns 'probe_id' and 'sequence'.", call. = FALSE)
  }
  dt <- fread(text = paste(lines, collapse = "\n"), sep = "\t", header = TRUE)
  dt <- dt[!is.na(probe_id) & nzchar(probe_id)]
  dt[, sequence := norm_seq(sequence)]
  bad <- dt[!grepl("^[ACGT]+$", sequence)]
  if (nrow(bad) > 0) {
    stop(sprintf("probe_file has invalid sequences for: %s (only A/C/G/T/U allowed).",
                 paste(bad$probe_id, collapse = ", ")), call. = FALSE)
  }
  unique(dt[, .(probe_id = as.character(probe_id), sequence = sequence)])
}
probes <- read_probes(opt$probe_file)

revcomp1 <- function(s) {
  s <- chartr("ACGT", "TGCA", s)
  paste(rev(strsplit(s, "")[[1]]), collapse = "")
}

# All substrings of length >= min_len of one sequence (and, with mismatches = 1,
# their 1-mismatch variants), each with its 0-based offset in the probe.
make_variants <- function(seqs, min_len, mismatches) {
  bases <- c("A", "C", "G", "T")
  one <- function(s, probe_id, orientation) {
    n <- nchar(s)
    if (n < min_len) return(NULL)
    out <- list()
    for (L in min_len:n) {
      starts <- seq_len(n - L + 1)
      subs <- substring(s, starts, starts + L - 1)
      offs <- starts - 1L
      out[[length(out) + 1]] <- data.table(variant = subs, offset = offs, length = L)
      if (mismatches >= 1) {
        for (p in seq_len(L)) {
          orig <- substr(subs, p, p)
          for (b in bases) {
            idx <- which(orig != b)
            if (length(idx) == 0) next
            vsub <- subs[idx]
            substr(vsub, p, p) <- b
            out[[length(out) + 1]] <- data.table(variant = vsub, offset = offs[idx], length = L)
          }
        }
      }
    }
    if (length(out) == 0) return(NULL)
    dt <- rbindlist(out)
    dt[, `:=`(probe_id = probe_id, orientation = orientation)]
    dt
  }
  res <- list()
  for (i in seq_len(nrow(seqs))) {
    res[[length(res) + 1]] <- one(seqs$sequence[i], seqs$probe_id[i], "as_given")
    res[[length(res) + 1]] <- one(revcomp1(seqs$sequence[i]), seqs$probe_id[i], "reverse_complement")
  }
  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) == 0) {
    return(data.table(variant = character(), offset = integer(), length = integer(),
                      probe_id = character(), orientation = character()))
  }
  rbindlist(res)
}

variants <- if (nrow(probes) > 0) make_variants(probes, opt$probe_min_len, opt$probe_mismatches) else
  data.table(variant = character(), offset = integer(), length = integer(),
             probe_id = character(), orientation = character())
message(sprintf("read_probe: %d probes, %d variant strings.", nrow(probes), nrow(variants)))

# Collapsed FASTA: two lines per record, header '>seq_N_xCOUNT rnatype:CLASS'.
parse_fasta <- function(fa) {
  lines <- readLines(fa, warn = FALSE)
  is_hdr <- startsWith(lines, ">")
  hdr  <- lines[is_hdr]
  seqs <- toupper(lines[!is_hdr])
  name  <- sub("^>", "", hdr)
  count <- suppressWarnings(as.numeric(sub(".*_x", "", sub(" .*$", "", name))))
  rt    <- ifelse(grepl("rnatype:", hdr), trimws(sub(".*rnatype:", "", hdr)), NA_character_)
  dt <- data.table(seq = seqs, count = count, rnatype = rt)
  dt[!is.na(seq) & nzchar(seq)]
}

summary_rows <- list(); hits_rows <- list(); unknown_rows <- list()
total_by_sample   <- setNames(numeric(0), character(0))
unknown_by_sample <- setNames(numeric(0), character(0))

for (s in samples) {
  fa <- file.path(opt$nfcore_dir, "mirtrace", s,
                  "qc_passed_reads.all.collapsed", paste0(s, ".fastp.fasta"))
  if (!file.exists(fa)) stop(sprintf("Missing collapsed FASTA: %s", fa), call. = FALSE)
  reads <- parse_fasta(fa)
  reads <- reads[, .(count = sum(count), rnatype = rnatype[1]), by = seq]
  total_reads   <- sum(reads$count)
  unknown_reads <- sum(reads$count[!is.na(reads$rnatype) & reads$rnatype == "unknown"])
  total_by_sample[s]   <- total_reads
  unknown_by_sample[s] <- unknown_reads

  if (opt$top_unknown_n > 0 && unknown_reads > 0) {
    unk <- head(reads[!is.na(rnatype) & rnatype == "unknown"][order(-count)], opt$top_unknown_n)
    unknown_rows[[s]] <- data.table(
      sample_id = s, rank = seq_len(nrow(unk)), read_sequence = unk$seq,
      length = nchar(unk$seq), count = unk$count, frac_unknown = unk$count / unknown_reads)
  }

  if (nrow(variants) == 0) next
  matched <- reads[seq %in% variants$variant]
  if (nrow(matched) == 0) next

  hit <- variants[matched, on = .(variant = seq), nomatch = 0, allow.cartesian = TRUE]
  # A read counts once per (probe, orientation) even if it matches at several offsets.
  uniq <- unique(hit, by = c("probe_id", "orientation", "variant"))
  summary_rows[[s]] <- uniq[, .(sample_id = s, n_reads = sum(count)), by = .(probe_id, orientation)]

  top <- uniq[order(-count), .SD[seq_len(min(opt$top_hits_n, .N))],
              by = .(probe_id, orientation)]
  top[, rank := seq_len(.N), by = .(probe_id, orientation)]
  hits_rows[[s]] <- top[, .(sample_id = s, probe_id, orientation, rank,
                            read_sequence = variant, length = nchar(variant),
                            count, rnatype, offset_in_probe = offset)]
}

dir.create("deliverables/probes", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------
# probe_summary.tsv: full sample x probe x orientation grid, 0 where no reads
# ------------------------------------------------------------------------
if (nrow(probes) > 0) {
  grid <- CJ(sample_id = samples, probe_id = probes$probe_id,
             orientation = c("as_given", "reverse_complement"), sorted = FALSE)
  if (length(summary_rows) > 0) {
    grid <- merge(grid, rbindlist(summary_rows),
                  by = c("sample_id", "probe_id", "orientation"), all.x = TRUE)
  } else {
    grid[, n_reads := NA_integer_]
  }
  grid[is.na(n_reads), n_reads := 0]
  grid[, frac_total   := n_reads / pmax(total_by_sample[sample_id], 1)]
  grid[, frac_unknown := ifelse(unknown_by_sample[sample_id] > 0,
                                n_reads / unknown_by_sample[sample_id], NA_real_)]
  setcolorder(grid, c("sample_id", "probe_id", "orientation", "n_reads", "frac_total", "frac_unknown"))
  setorder(grid, sample_id, probe_id, orientation)
  fwrite(grid, "deliverables/probes/probe_summary.tsv", sep = "\t", quote = FALSE, na = "")
}

if (length(hits_rows) > 0) {
  fwrite(rbindlist(hits_rows), "deliverables/probes/probe_hits.tsv",
         sep = "\t", quote = FALSE, na = "")
}

if (opt$top_unknown_n > 0 && length(unknown_rows) > 0) {
  fwrite(rbindlist(unknown_rows), "deliverables/probes/top_unknown_by_sample.tsv",
         sep = "\t", quote = FALSE, na = "")
}

message(sprintf("read_probe ok: %d probes, %d samples, top_unknown_n=%d.",
                nrow(probes), length(samples), opt$top_unknown_n))
