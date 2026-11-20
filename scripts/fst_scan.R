#!/usr/bin/env Rscript

# Genome-wide pairwise Fst between two groups defined in the metadata.
#
# Two estimators are available:
#
#   weir-cockerham (default)
#     The Weir & Cockerham (1984) theta estimator, using the a/b/c variance
#     components. For haploid data b and c vanish, so theta reduces to a / (a + ...)
#     with the sample-size correction retained. This is the estimator to cite
#     when Methods say "Weir and Cockerham's Fst".
#
#   wright
#     A sample-size-weighted Wright's Fst:
#       (n1(p1 - pbar)^2 + n2(p2 - pbar)^2) / ((n1 + n2) pbar (1 - pbar))
#     Simpler, uncorrected for finite sampling, and biased upward when group
#     sizes are small or uneven. Provided for continuity with earlier analyses.
#
# The two give similar rankings but different absolute values, so do not mix
# them within a manuscript, and name the one used.
#
# Input is the binary SNP matrix (chr, pos, ref, <samples>) with 0 = reference,
# 1 = alternate, 0.5 = mixed and N / . = missing.

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Output directory [default %default]", metavar = "character"),
  make_option(c("-b", "--matrix_binary"), type = "character", default = NULL,
              help = "Binary SNP matrix: chr, pos, ref, <samples> [required]",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Metadata TSV [required]", metavar = "character"),
  make_option(c("--label_id"), type = "character", default = "sample_id",
              help = "Metadata column holding sample IDs [default %default]",
              metavar = "character"),
  make_option(c("--label_category"), type = "character", default = "country",
              help = "Metadata column defining the groups [default %default]",
              metavar = "character"),
  make_option(c("--pop1"), type = "character", default = NULL,
              help = "Comma-separated category value(s) forming population 1 [required]",
              metavar = "character"),
  make_option(c("--pop2"), type = "character", default = NULL,
              help = paste("Comma-separated category value(s) forming population 2.",
                           "If omitted, every category except pop1 is used."),
              metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = NULL,
              help = "Optional Fws column; if given, samples below --fws_th are dropped",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold, applied only with --label_fws [default %default]",
              metavar = "numeric"),
  make_option(c("--estimator"), type = "character", default = "weir-cockerham",
              help = "weir-cockerham or wright [default %default]",
              metavar = "character"),
  make_option(c("--max_missing"), type = "numeric", default = 0.2,
              help = "Drop a SNP if either population exceeds this missing fraction [default %default]",
              metavar = "numeric"),
  make_option(c("--maf"), type = "numeric", default = 0.01,
              help = "Drop SNPs with pooled MAF below this [default %default]",
              metavar = "numeric"),
  make_option(c("--mixed_as"), type = "character", default = "alt",
              help = "Treat 0.5 calls as 'alt' or 'missing' [default %default]",
              metavar = "character"),
  make_option(c("--na_char"), type = "character", default = "N",
              help = "Missing-genotype character in the matrix [default %default]",
              metavar = "character"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Comma-separated chromosomes to drop (e.g. Pf3D7_API_v3,Pf3D7_MIT_v3)",
              metavar = "character"),
  make_option(c("--gene-file"), type = "character", default = NULL,
              help = "Optional genome product TSV; annotates SNPs above --fst_thresh",
              metavar = "character"),
  make_option(c("--fst_thresh"), type = "numeric", default = 0.95,
              help = "Fst above which a SNP is reported as extreme [default %default]",
              metavar = "numeric"),
  make_option(c("--out_prefix"), type = "character", default = NULL,
              help = "Output file prefix (default: <pop1>_vs_<pop2>)",
              metavar = "character"),
  make_option(c("--threads"), type = "integer", default = 4,
              help = "data.table threads [default %default]", metavar = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$matrix_binary) || is.null(opt$metadata) || is.null(opt$pop1)) {
  stop("--matrix_binary, --metadata and --pop1 are required.")
}
if (!opt$estimator %in% c("weir-cockerham", "wright")) {
  stop("--estimator must be 'weir-cockerham' or 'wright'.")
}
if (!opt$mixed_as %in% c("alt", "missing")) {
  stop("--mixed_as must be 'alt' or 'missing'.")
}

setDTthreads(opt$threads)
dir.create(opt$workdir, showWarnings = FALSE, recursive = TRUE)
split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

# ---- metadata and group membership ----------------------------------
meta <- read.csv(opt$metadata, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
for (col in c(opt$label_id, opt$label_category)) {
  if (!col %in% colnames(meta)) stop("Metadata has no column '", col, "'.")
}
if (!is.null(opt$label_fws)) {
  if (!opt$label_fws %in% colnames(meta)) stop("Metadata has no column '", opt$label_fws, "'.")
  meta <- meta[!is.na(meta[[opt$label_fws]]) & meta[[opt$label_fws]] >= opt$fws_th, ]
  message("Fws >= ", opt$fws_th, ": ", nrow(meta), " samples retained")
}

p1_vals <- split_csv(opt$pop1)
p2_vals <- if (is.null(opt$pop2)) setdiff(unique(meta[[opt$label_category]]), p1_vals) else split_csv(opt$pop2)
if (!length(p2_vals)) stop("Population 2 is empty.")

s1 <- meta[[opt$label_id]][meta[[opt$label_category]] %in% p1_vals]
s2 <- meta[[opt$label_id]][meta[[opt$label_category]] %in% p2_vals]
if (!length(s1)) stop("No samples match --pop1.")
if (!length(s2)) stop("No samples match --pop2.")

# ---- genotypes -------------------------------------------------------
message("Reading matrix...")
snp <- fread(opt$matrix_binary, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("chr", "pos") %in% colnames(snp)[1:3])) {
  stop("Matrix must begin with chr, pos, ref columns.")
}
if (!is.null(opt$remove_chr)) {
  snp <- snp[!snp$chr %in% split_csv(opt$remove_chr), , drop = FALSE]
}

s1 <- intersect(s1, colnames(snp))
s2 <- intersect(s2, colnames(snp))
s2 <- setdiff(s2, s1)                       # a sample can only sit in one group
if (!length(s1) || !length(s2)) stop("No overlap between metadata samples and matrix columns.")
message("pop1: ", length(s1), " samples | pop2: ", length(s2), " samples")

legend <- snp[, 1:3, drop = FALSE]
geno <- snp[, c(s1, s2), drop = FALSE]
rm(snp)

geno[geno == opt$na_char] <- NA
geno[geno == "."] <- NA
geno[] <- lapply(geno, function(x) suppressWarnings(as.numeric(x)))
if (opt$mixed_as == "alt") geno[geno == 0.5] <- 1 else geno[geno == 0.5] <- NA
G <- as.matrix(geno); rm(geno)

i1 <- seq_along(s1)
i2 <- seq_along(s2) + length(s1)

miss1 <- rowMeans(is.na(G[, i1, drop = FALSE]))
miss2 <- rowMeans(is.na(G[, i2, drop = FALSE]))
keep <- miss1 <= opt$max_missing & miss2 <= opt$max_missing
message(sum(keep), " / ", length(keep), " SNPs pass missingness <= ", opt$max_missing)

G <- G[keep, , drop = FALSE]; legend <- legend[keep, , drop = FALSE]

n1 <- rowSums(!is.na(G[, i1, drop = FALSE]))
n2 <- rowSums(!is.na(G[, i2, drop = FALSE]))
p1 <- rowMeans(G[, i1, drop = FALSE], na.rm = TRUE)
p2 <- rowMeans(G[, i2, drop = FALSE], na.rm = TRUE)

ok <- n1 > 0 & n2 > 0
pbar_pool <- (n1 * p1 + n2 * p2) / pmax(n1 + n2, 1)
ok <- ok & pmin(pbar_pool, 1 - pbar_pool) >= opt$maf
message(sum(ok), " SNPs pass pooled MAF >= ", opt$maf)
if (!any(ok)) stop("No SNPs survive filtering.")

legend <- legend[ok, , drop = FALSE]
n1 <- n1[ok]; n2 <- n2[ok]; p1 <- p1[ok]; p2 <- p2[ok]; pbar_pool <- pbar_pool[ok]

# ---- Fst -------------------------------------------------------------
if (opt$estimator == "wright") {
  fst <- (n1 * (p1 - pbar_pool)^2 + n2 * (p2 - pbar_pool)^2) /
         ((n1 + n2) * pbar_pool * (1 - pbar_pool))
} else {
  # Weir & Cockerham (1984) theta, haploid form: b and c are zero, and the
  # within-population term carries the sample-size correction.
  r <- 2
  nbar <- (n1 + n2) / r
  nc <- (r * nbar - (n1^2 + n2^2) / (r * nbar)) / (r - 1)
  pbar <- (n1 * p1 + n2 * p2) / (r * nbar)
  s2_ <- ((n1 * (p1 - pbar)^2) + (n2 * (p2 - pbar)^2)) / ((r - 1) * nbar)
  a <- (nbar / nc) * (s2_ - (1 / (nbar - 1)) * (pbar * (1 - pbar) - ((r - 1) / r) * s2_))
  b <- (nbar / (nbar - 1)) * (pbar * (1 - pbar) - ((r - 1) / r) * s2_)
  fst <- a / (a + b)
}
fst[!is.finite(fst)] <- NA_real_
fst <- pmin(pmax(fst, 0), 1)                # clamp small negative estimates

res <- data.frame(legend,
                  n_pop1 = n1, n_pop2 = n2,
                  freq_pop1 = signif(p1, 4), freq_pop2 = signif(p2, 4),
                  fst = signif(fst, 5),
                  stringsAsFactors = FALSE)

prefix <- opt$out_prefix
if (is.null(prefix)) {
  lab <- function(v) if (length(v) > 3) paste0(length(v), "pops") else paste(v, collapse = "-")
  prefix <- paste0(lab(p1_vals), "_vs_", lab(p2_vals))
}
out_all <- file.path(opt$workdir, paste0(prefix, ".fst.", opt$estimator, ".tsv"))
fwrite(res, out_all, sep = "\t", quote = FALSE)

message("\nFst summary (", opt$estimator, "):")
print(summary(res$fst))
message("mean Fst: ", signif(mean(res$fst, na.rm = TRUE), 5))

extreme <- res[!is.na(res$fst) & res$fst >= opt$fst_thresh, , drop = FALSE]
message(nrow(extreme), " SNPs with Fst >= ", opt$fst_thresh)

# ---- optional gene annotation ---------------------------------------
if (!is.null(opt$`gene-file`) && nrow(extreme)) {
  gp <- tryCatch(
    read.csv(opt$`gene-file`, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  need <- c("chr", "pos_start", "pos_end")
  if (!is.null(gp) && all(need %in% colnames(gp))) {
    name_col <- intersect(c("gene_name", "gene_id", "Gene_name"), colnames(gp))[1]
    extreme$gene <- vapply(seq_len(nrow(extreme)), function(i) {
      hit <- gp[gp$chr == extreme$chr[i] &
                gp$pos_start <= extreme$pos[i] &
                gp$pos_end >= extreme$pos[i], , drop = FALSE]
      if (!nrow(hit) || is.na(name_col)) "" else paste(unique(hit[[name_col]]), collapse = ";")
    }, character(1))
  } else {
    message("Gene file missing chr/pos_start/pos_end columns; skipping annotation.")
  }
}
if (nrow(extreme)) {
  fwrite(extreme, file.path(opt$workdir, paste0(prefix, ".fst.extreme.tsv")),
         sep = "\t", quote = FALSE)
}

message("\nWritten to ", opt$workdir, ":")
message("  ", basename(out_all), "   all SNPs with per-population frequencies and Fst")
if (nrow(extreme)) message("  ", prefix, ".fst.extreme.tsv   SNPs with Fst >= ", opt$fst_thresh)
message("\nNOTE: estimator = ", opt$estimator,
        ". Name this in the Methods; 'Weir and Cockerham' applies only to weir-cockerham.")
