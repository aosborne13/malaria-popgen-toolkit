#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ─────────────────────────────────────────────────────────────
# Helper: MAF calculation for hmmIBD-coded matrix
#   0 / 1 / -1 (missing)
# ─────────────────────────────────────────────────────────────
calculate_maf <- function(mat) {
  m <- as.matrix(mat)
  m[m < 0] <- NA
  p <- rowMeans(m == 1, na.rm = TRUE)
  p[is.nan(p)] <- NA_real_
  pmin(p, 1 - p)
}

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────
option_list <- list(

  make_option("--outdir", type = "character", default = ".",
              help = "Output directory [default %default]"),

  make_option("--binary_matrix", type = "character",
              help = "Binary SNP matrix (chr,pos,ref,<samples>)",
              metavar = "file"),

  make_option("--metadata", type = "character",
              help = "Metadata TSV file",
              metavar = "file"),

  make_option("--category", type = "character", default = NULL,
              help = "Category name(s) to run (comma-separated) or omit for ALL"),

  make_option("--label_category", type = "character", default = "country",
              help = "Metadata column defining category [default %default]"),

  make_option("--label_id", type = "character", default = "sample_id",
              help = "Metadata sample ID column [default %default]"),

  make_option("--label_fws", type = "character", default = "fws",
              help = "Metadata Fws column [default %default]"),

  make_option("--fws_th", type = "numeric", default = 0.95,
              help = "Minimum Fws threshold [default %default]"),

  make_option("--maf", type = "numeric", default = 0.01,
              help = "Minor allele frequency threshold [default %default]"),

  make_option("--threads", type = "integer", default = 4,
              help = "data.table threads [default %default]"),

  make_option("--remove_chr", type = "character", default = NULL,
              help = "Comma-separated chromosomes to remove"),

  make_option("--regex_chr", type = "character",
              default = "(.*?)_(.+)_(.*)",
              help = "Regex to parse chromosome names"),

  make_option("--regex_groupid", type = "integer", default = 3,
              help = "Regex capture group for numeric chromosome"),

  make_option("--hmmibd_bin", type = "character", default = "hmmIBD",
              help = "Path to hmmIBD binary"),

  make_option("--skip_hmmibd", action = "store_true", default = FALSE,
              help = "Prepare inputs but do NOT run hmmIBD")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ─────────────────────────────────────────────────────────────
# Validate input
# ─────────────────────────────────────────────────────────────
if (is.null(opt$binary_matrix) || is.null(opt$metadata)) {
  stop("ERROR: --binary_matrix and --metadata are required", call. = FALSE)
}

setDTthreads(opt$threads)

outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("Writing outputs to: ", normalizePath(outdir))

# ─────────────────────────────────────────────────────────────
# Load metadata
# ─────────────────────────────────────────────────────────────
meta <- fread(opt$metadata, sep = "\t", data.table = FALSE)

needed <- c(opt$label_id, opt$label_category, opt$label_fws)
missing_cols <- setdiff(needed, colnames(meta))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing metadata columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

meta[[opt$label_id]] <- as.character(meta[[opt$label_id]])
meta[[opt$label_category]] <- as.character(meta[[opt$label_category]])

all_categories <- sort(unique(meta[[opt$label_category]]))

# Decide which categories to run
if (is.null(opt$category) || opt$category == "" || tolower(opt$category) == "all") {
  categories <- all_categories
  message("Running ALL categories: ", paste(categories, collapse = ", "))
} else {
  categories <- strsplit(opt$category, ",")[[1]] |> trimws()
  bad <- setdiff(categories, all_categories)
  if (length(bad) > 0) {
    stop("ERROR: Unknown category: ", paste(bad, collapse = ", "), call. = FALSE)
  }
}

# ─────────────────────────────────────────────────────────────
# Load SNP matrix
# ─────────────────────────────────────────────────────────────
message("Loading SNP matrix...")
snp_all <- fread(opt$binary_matrix, sep = "\t", data.table = FALSE)

stopifnot(all(c("chr", "pos", "ref") %in% names(snp_all)[1:3]))

# Remove unwanted chromosomes
if (!is.null(opt$remove_chr)) {
  rm_chr <- strsplit(opt$remove_chr, ",")[[1]] |> trimws()
  snp_all <- snp_all[!snp_all$chr %in% rm_chr, , drop = FALSE]
  message("Removed chromosomes: ", paste(rm_chr, collapse = ", "))
}

# Convert chromosome to numeric for hmmIBD
chr_match <- str_match(snp_all$chr, opt$regex_chr)
snp_all$chr <- as.numeric(chr_match[, opt$regex_groupid])

# Write global haplotype legend
hap_leg <- snp_all[, c("chr", "pos")]
fwrite(hap_leg,
       file.path(outdir, "ibd_matrix_hap_leg.tsv"),
       sep = "\t")

# Genotype columns
sample_cols <- setdiff(colnames(snp_all), c("chr", "pos", "ref"))

# ─────────────────────────────────────────────────────────────
# Process each category
# ─────────────────────────────────────────────────────────────
for (cat in categories) {

  message("\n==============================")
  message("Category: ", cat)
  message("==============================")

  meta_cat <- meta %>%
    filter(.data[[opt$label_category]] == cat &
           !is.na(.data[[opt$label_fws]]) &
           .data[[opt$label_fws]] >= opt$fws_th)

  if (nrow(meta_cat) < 2) {
    warning("Skipping ", cat, ": <2 samples after Fws filter")
    next
  }

  samples <- intersect(meta_cat[[opt$label_id]], sample_cols)

  if (length(samples) < 2) {
    warning("Skipping ", cat, ": insufficient overlapping samples")
    next
  }

  snp_cat <- snp_all[, c("chr", "pos", "ref", samples), drop = FALSE]
  geno <- snp_cat[, -(1:3), drop = FALSE]

  geno[geno %in% c("NA", "N", ".")] <- NA
  geno[] <- lapply(geno, as.numeric)

  hmm <- geno
  hmm[hmm == 0.5] <- 1
  hmm[is.na(hmm)] <- -1

  maf_vals <- calculate_maf(hmm)
  keep <- which(!is.na(maf_vals) & maf_vals >= opt$maf)

  if (length(keep) == 0) {
    warning("No SNPs pass MAF in ", cat)
    next
  }

  message("Keeping ", length(keep), " SNPs for ", cat,
          " (", length(samples), " samples)")

  cat_safe <- gsub("[^A-Za-z0-9_.-]", "_", cat)
  maf_lbl <- format(opt$maf, scientific = FALSE)

  out_mat <- cbind(
    chrom = snp_cat$chr[keep],
    pos   = snp_cat$pos[keep],
    hmm[keep, , drop = FALSE]
  )

  input_file <- file.path(
    outdir,
    sprintf("hmmIBD_%s_maf%s.txt", cat_safe, maf_lbl)
  )

  fwrite(out_mat, input_file, sep = "\t")
  message("Wrote hmmIBD input: ", input_file)

  if (opt$skip_hmmibd) {
    message("skip_hmmibd = TRUE, not running hmmIBD")
    next
  }

  out_prefix <- sprintf("hmmIBD_%s_maf%s_out", cat_safe, maf_lbl)
  cmd <- sprintf(
    "cd %s && %s -i %s -o %s",
    shQuote(outdir),
    shQuote(opt$hmmibd_bin),
    shQuote(basename(input_file)),
    shQuote(out_prefix)
  )

  message("Running hmmIBD...")
  system(cmd)
}

message("\nDone.")

