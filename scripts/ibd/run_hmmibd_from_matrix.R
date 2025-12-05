#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ─────────────────────────────────────────────────────────────
# Helper: MAF calculation (0/1/-1 matrix: -1 = missing)
# ─────────────────────────────────────────────────────────────
calculate_maf <- function(mat) {
  m <- as.matrix(mat)
  m[m < 0] <- NA       # treat -1 as missing
  p <- rowMeans(m == 1, na.rm = TRUE)  # alt allele frequency
  p[is.nan(p)] <- NA_real_
  maf <- pmin(p, 1 - p)
  maf
}

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Output / working directory [default %default]",
              metavar = "character"),
  make_option(c("-b", "--binary_matrix"), type = "character", default = NULL,
              help = "Input filename of filtered binary matrix [required]",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Full path to metadata TSV file [required]",
              metavar = "character"),
  make_option(c("-c", "--category"), type = "character", default = NULL,
              help = "Category name (e.g. Ethiopia), comma-separated list, or omit/empty for ALL",
              metavar = "character"),
  make_option(c("--label_category"), type = "character", default = "country",
              help = "Column name in metadata for category [default %default]",
              metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = "fws",
              help = "Column name in metadata for Fws values [default %default]",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold (keep samples with Fws >= this) [default %default]",
              metavar = "number"),
  make_option(c("--label_id"), type = "character", default = "sra_run",
              help = "Column name in metadata for sample ID [default %default]",
              metavar = "character"),
  make_option(c("--maf"), type = "numeric", default = 0.01,
              help = "MAF threshold [default %default]",
              metavar = "number"),
  make_option(c("--na_char"), type = "character", default = "NA",
              help = "Character used for missing genotypes (also treats '.' and 'N' as missing)",
              metavar = "character"),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Number of threads for data.table [default %default]",
              metavar = "numeric"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Comma-separated chromosomes to remove, e.g. Pf3D7_API_v3,Pf3D7_MIT_v3",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection (default matches Pf3D7_01_v3)",
              metavar = "character"),
  make_option("--regex_groupid", type = "integer", default = 3,
              help = "Capture group index in regex_chr giving numeric chromosome [default %default]",
              metavar = "numeric"),
  make_option("--hmmibd_bin", type = "character", default = "hmmIBD",
              help = "Path to hmmIBD binary (must be executable) [default %default]",
              metavar = "character"),
  make_option("--skip_hmmibd", action = "store_true", default = FALSE,
              help = "Only prepare hmmIBD input files; do NOT run hmmIBD")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ─────────────────────────────────────────────────────────────
# Extract options
# ─────────────────────────────────────────────────────────────
workdir        <- opt$workdir
bin_mat_file   <- opt$binary_matrix
met_file       <- opt$metadata
category_arg   <- opt$category
label_category <- opt$label_category
label_id       <- opt$label_id
label_fws      <- opt$label_fws
threshold_fws  <- opt$fws_th
th_maf         <- opt$maf
threads        <- opt$threads
na_char        <- opt$na_char
rm_chr_str     <- opt$remove_chr
pattern        <- opt$regex_chr
groupid        <- opt$regex_groupid
hmmibd_bin     <- path.expand(opt$hmmibd_bin)
skip_hmmibd    <- isTRUE(opt$skip_hmmibd)

if (is.null(bin_mat_file)) stop("ERROR: --binary_matrix is required\n", call. = FALSE)
if (is.null(met_file))     stop("ERROR: --metadata is required\n", call. = FALSE)

dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
setDTthreads(threads)

message("Working directory: ", workdir)
message("Binary matrix: ", bin_mat_file)
message("Metadata: ", met_file)

# ─────────────────────────────────────────────────────────────
# Load metadata
# ─────────────────────────────────────────────────────────────
metadata <- read.csv(met_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

needed_cols <- c(label_category, label_id, label_fws)
missing_cols <- setdiff(needed_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing column(s) in metadata: ",
       paste(missing_cols, collapse = ", "), "\n", call. = FALSE)
}

metadata[[label_category]] <- as.character(metadata[[label_category]])
metadata[[label_id]]       <- as.character(metadata[[label_id]])

available_categories <- sort(unique(metadata[[label_category]]))

# Decide which categories to run
if (is.null(category_arg) || category_arg == "" || tolower(category_arg) == "all") {
  categories <- available_categories
  message("No specific category provided; running ALL categories: ",
          paste(categories, collapse = ", "))
} else {
  categories <- strsplit(category_arg, ",")[[1]]
  categories <- trimws(categories)
  missing_cats <- setdiff(categories, available_categories)
  if (length(missing_cats) > 0) {
    stop("ERROR: Requested category(ies) not found in metadata: ",
         paste(missing_cats, collapse = ", "), "\n",
         "Available: ", paste(available_categories, collapse = ", "), "\n",
         call. = FALSE)
  }
  message("Running category(ies): ", paste(categories, collapse = ", "))
}

# ─────────────────────────────────────────────────────────────
# Load full binary matrix
# ─────────────────────────────────────────────────────────────
message("Loading binary matrix (this may take a while)...")
snp_all <- fread(bin_mat_file, sep = "\t", header = TRUE, data.table = FALSE)

# Expect first three columns: chr, pos, ref
if (!all(c("chr", "pos", "ref") %in% colnames(snp_all)[1:3])) {
  stop("ERROR: Matrix must start with columns: chr, pos, ref\n", call. = FALSE)
}

# Remove specified chromosomes by original name
if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
  rm_chr <- strsplit(rm_chr_str, ",")[[1]] |> trimws()
  present_rm <- intersect(rm_chr, unique(snp_all$chr))
  if (length(present_rm) > 0) {
    message("Removing chromosomes: ", paste(present_rm, collapse = ", "))
    snp_all <- snp_all[!snp_all$chr %in% present_rm, , drop = FALSE]
  } else {
    warning("None of the chromosomes in --remove_chr were found in matrix; continuing.\n")
  }
} else {
  message("No chromosomes removed (API/mito are still present).")
}

# Convert chr string -> numeric according to regex
chr_match <- stringr::str_match(snp_all$chr, pattern)
if (ncol(chr_match) < groupid) {
  stop("ERROR: regex_chr did not produce enough capture groups; check pattern and --regex_groupid\n",
       call. = FALSE)
}
snp_all$chr <- as.numeric(chr_match[, groupid])

# Write shared haplotype legend (chrom,pos) BEFORE per-category filters
hap_leg <- snp_all[, c("chr", "pos")]
colnames(hap_leg) <- c("chrom", "pos")
hap_leg_file <- file.path(workdir, "ibd_matrix_hap_leg.tsv")
fwrite(hap_leg, hap_leg_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote haplotype legend: ", hap_leg_file)

# Genotype columns
all_sample_cols <- setdiff(colnames(snp_all), c("chr", "pos", "ref"))

# ─────────────────────────────────────────────────────────────
# Process categories
# ─────────────────────────────────────────────────────────────
for (cat_val in categories) {
  message("\n==============================")
  message("Category: ", cat_val)
  message("==============================")

  meta_cat <- metadata %>%
    filter(.data[[label_category]] == cat_val & !is.na(.data[[label_fws]]) & .data[[label_fws]] >= threshold_fws)

  if (nrow(meta_cat) < 2) {
    warning("Category ", cat_val, ": fewer than 2 samples pass Fws threshold; skipping.\n")
    next
  }

  samples_meta <- unique(meta_cat[[label_id]])
  samples_in_matrix <- intersect(samples_meta, all_sample_cols)
  missing_samples <- setdiff(samples_meta, all_sample_cols)

  if (length(missing_samples) > 0) {
    warning("Category ", cat_val, ": some metadata samples are missing from matrix: ",
            paste(missing_samples, collapse = ", "))
  }

  if (length(samples_in_matrix) < 2) {
    warning("Category ", cat_val, ": fewer than 2 overlapping samples between metadata and matrix; skipping.\n")
    next
  }

  # Subset matrix for this category: chr,pos,ref + category samples
  snp_cat <- snp_all[, c("chr", "pos", "ref", samples_in_matrix), drop = FALSE]

  # Recode missing
  geno <- snp_cat[, -(1:3), drop = FALSE]
  geno[geno == na_char] <- NA
  geno[geno == "."]     <- NA
  geno[geno == "N"]     <- NA

  # Convert to numeric (0, 0.5, 1, NA)
  geno[] <- lapply(geno, function(x) suppressWarnings(as.numeric(x)))

  # Recode for hmmIBD: 0 -> 0, 0.5 -> 1, 1 -> 1, NA -> -1
  hmm_geno <- geno
  hmm_geno[hmm_geno == 0.5] <- 1
  hmm_geno[is.na(hmm_geno)] <- -1

  # Write "raw" (pre-MAF/Fws) hap matrix for this category, if you still want it:
  cat_safe <- gsub("[^A-Za-z0-9_.-]", "_", cat_val)
  hap_cat_full <- cbind(
    chrom = snp_cat$chr,
    pos   = snp_cat$pos,
    hmm_geno
  )
  hap_cat_file <- file.path(workdir, sprintf("ibd_matrix_hap_%s.tsv", cat_safe))
  fwrite(hap_cat_full, hap_cat_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote category hap matrix (pre-MAF/Fws): ", hap_cat_file)

  # MAF filter using hmm-coded matrix (0/1/-1)
  maf_vals <- calculate_maf(hmm_geno)
  keep_idx <- which(!is.na(maf_vals) & maf_vals >= th_maf)

  if (length(keep_idx) == 0) {
    warning("Category ", cat_val, ": no SNPs pass MAF ≥ ", th_maf, "; skipping.\n")
    next
  }

  message("Category ", cat_val, ": ",
          length(keep_idx), " SNPs pass MAF ≥ ", th_maf,
          " with ", length(samples_in_matrix), " samples (Fws ≥ ", threshold_fws, ").")

  hap_leg_cat <- hap_leg[keep_idx, , drop = FALSE]
  hmm_geno_cat <- hmm_geno[keep_idx, , drop = FALSE]

  maf_label <- format(th_maf, trim = TRUE, scientific = FALSE)
  input_base <- sprintf("hmmIBD_%s_maf%s.txt", cat_safe, maf_label)
  input_file <- file.path(workdir, input_base)

  hmm_input <- cbind(hap_leg_cat, hmm_geno_cat)
  fwrite(hmm_input, input_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote hmmIBD input for ", cat_val, ": ", input_file)

  # Optionally run hmmIBD
  if (skip_hmmibd) {
    message("skip_hmmibd = TRUE: not running hmmIBD for ", cat_val)
    next
  }

  out_prefix_base <- sprintf("hmmIBD_%s_maf%s_out", cat_safe, maf_label)

  # Run from inside workdir with relative paths to avoid AV issues with full paths.
  cmd_inner <- sprintf("%s -i %s -o %s",
                       shQuote(hmmibd_bin),
                       shQuote(input_base),
                       shQuote(out_prefix_base))
  full_cmd <- sprintf("cd %s && %s", shQuote(workdir), cmd_inner)

  message("Running hmmIBD for ", cat_val, ":")
  message("  ", full_cmd)

  res <- tryCatch(
    system(full_cmd, intern = TRUE),
    warning = function(w) w,
    error   = function(e) e
  )

  log_file <- file.path(workdir, sprintf("hmmIBD_run_%s.log", cat_safe))

  if (inherits(res, "error")) {
    warning("hmmIBD failed for category ", cat_val,
            "; see log: ", log_file, "\n  Error: ", conditionMessage(res))
    writeLines(c("COMMAND:", full_cmd,
                 "STATUS: ERROR",
                 conditionMessage(res)),
               con = log_file)
  } else {
    status <- attr(res, "status")
    status_msg <- if (is.null(status) || status == 0) "OK" else paste("EXIT", status)
    writeLines(c("COMMAND:", full_cmd,
                 paste("STATUS:", status_msg),
                 "",
                 res),
               con = log_file)
    message("hmmIBD completed for ", cat_val,
            "; log written to ", log_file)
  }
}

message("\nDone.")
