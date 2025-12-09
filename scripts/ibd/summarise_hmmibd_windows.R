#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ─────────────────────────────────────────────────────────────
# Helper: annotate_candidate_regions
#   Replaces malaria-hub helpers.R dependency.
#   regions: data.frame(chr, start, end)
#   annotation: data.frame(chr, pos_start, pos_end, gene_id, gene_name, gene_product)
# ─────────────────────────────────────────────────────────────
annotate_candidate_regions <- function(regions, annotation) {
  reg_dt <- as.data.table(regions)
  setnames(reg_dt, c("chr", "start", "end"))

  ann_dt <- as.data.table(annotation)
  setnames(ann_dt, c("chr", "pos_start", "pos_end"),
           c("chr", "start", "end"))

  setkey(reg_dt, chr, start, end)
  setkey(ann_dt, chr, start, end)

  ov <- foverlaps(ann_dt, reg_dt, type = "any", nomatch = 0L)

  # One row per gene per IBD window
  res <- ov[, .(
    chr        = i.chr,
    start      = i.start,
    end        = i.end,
    pos_start  = start,
    pos_end    = end,
    gene_id    = gene_id,
    gene_name  = gene_name,
    gene_product = gene_product
  )]

  as.data.frame(res)
}

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = NULL,
              help = "Main hmmIBD output directory (same as --outdir used in hmmibd-matrix)",
              metavar = "character"),
  make_option(c("--list_category"), type = "character", default = NULL,
              help = "File with one category label per line (e.g. Ethiopia_2013).\nIf omitted, auto-detect subdirectories in workdir.",
              metavar = "character"),
  make_option(c("-l", "--legend"), type = "character",
              default = "ibd_matrix_hap_leg.tsv",
              help = "SNP legend file (default: ibd_matrix_hap_leg.tsv in workdir)",
              metavar = "character"),
  make_option("--gene_product", type = "character", default = NULL,
              help = "Gene product (annotation) TSV with columns: chr, pos_start, pos_end, gene_id, gene_name, gene_product",
              metavar = "character"),
  make_option(c("-r", "--ref_index"), type = "character", default = NULL,
              help = "Reference FASTA index (.fai)",
              metavar = "character"),
  make_option("--maf", type = "numeric", default = 0.01,
              help = "MAF threshold used in hmmIBD inputs (default: 0.01)",
              metavar = "numeric"),
  make_option(c("-w", "--window_size"), type = "numeric", default = 50000,
              help = "Window size in bp for summarising IBD (default: 50000)",
              metavar = "numeric"),
  make_option("--quantile_cutoff", type = "numeric", default = 0.95,
              help = "Quantile cut-off for annotated IBD segments (default: 0.95)",
              metavar = "numeric"),
  make_option("--suffix", type = "character",
              default = format(Sys.time(), "%d_%m_%Y"),
              help = "Prefix/suffix for summary files (default: current date dd_mm_YYYY)",
              metavar = "character"),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Number of threads for data.table (default: 4)",
              metavar = "numeric"),
  make_option("--remove_chr", type = "character", default = NULL,
              help = "Comma-separated chromosomes to remove, e.g. Pf3D7_API_v3,Pf3D7_MIT_v3",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection (default matches Pf3D7_01_v3)",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Capture group index in regex_chr for numeric chromosome (default: 3)",
              metavar = "numeric"),
  make_option("--NSNP", type = "numeric", default = 0,
              help = "Minimum SNPs per segment (default: 0 = no filter)",
              metavar = "numeric"),
  make_option("--LSEGMENT", type = "numeric", default = 0,
              help = "Minimum segment length in bp (default: 0 = no filter)",
              metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ─────────────────────────────────────────────────────────────
# Extract options
# ─────────────────────────────────────────────────────────────
workdir        <- opt$workdir
category_list  <- opt$list_category
legend         <- opt$legend
gene_product_file <- opt$gene_product
ref_index      <- opt$ref_index
th_maf         <- opt$maf
window_size    <- opt$window_size
th_quantile    <- opt$quantile_cutoff
suffix         <- opt$suffix
threads        <- opt$threads
rm_chr         <- opt$remove_chr
pattern        <- opt$regex_chr
groupid        <- opt$regex_groupid
min_n_snp      <- opt$NSNP
min_l_seg      <- opt$LSEGMENT

if (is.null(workdir)) {
  stop("ERROR: --workdir is required\n", call. = FALSE)
}

setDTthreads(threads)
message("Workdir: ", workdir)

# ─────────────────────────────────────────────────────────────
# Load category list
# - If --list_category is a file → read it
# - If NULL → auto-detect subdirectories in workdir
# ─────────────────────────────────────────────────────────────
if (!is.null(category_list)) {
  if (file.exists(category_list)) {
    category_list <- read.table(category_list, header = FALSE,
                                stringsAsFactors = FALSE)$V1
  } else {
    stop("Cannot locate file with category list: ", category_list, "\n", call. = FALSE)
  }
} else {
  # auto-detect: immediate subdirectories in workdir (e.g. Ethiopia_2013, Ethiopia_2017,…)
  subdirs <- list.dirs(workdir, full.names = FALSE, recursive = FALSE)
  category_list <- subdirs[nzchar(subdirs)]
  if (length(category_list) == 0) {
    stop("No subdirectories found in workdir and no --list_category provided.\n",
         "Expected per-category folders created by hmmibd-matrix (e.g. Ethiopia_2013/).",
         call. = FALSE)
  }
  message("Auto-detected categories: ", paste(category_list, collapse = ", "))
}

# ─────────────────────────────────────────────────────────────
# Load SNP legend
# ─────────────────────────────────────────────────────────────
legend_path <- legend
if (!file.exists(legend_path)) {
  # try relative to workdir
  candidate <- file.path(workdir, legend)
  if (file.exists(candidate)) {
    legend_path <- candidate
  }
}

if (!file.exists(legend_path)) {
  stop("Cannot locate SNP legend file: ", legend, "\n", call. = FALSE)
}

snp_hmmibd_02_1 <- read_tsv(legend_path, col_types = cols()) %>% as.data.frame()
message("Loaded legend: ", legend_path)

# ─────────────────────────────────────────────────────────────
# Load reference index (.fai)
# ─────────────────────────────────────────────────────────────
if (is.null(ref_index) || !file.exists(ref_index)) {
  stop("ERROR: --ref_index must point to an existing FASTA .fai file\n", call. = FALSE)
}

fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  rename(chr = V1, end_chr = V2) %>%
  select(chr, end_chr) %>%
  mutate(start_chr = 1)

# Remove chromosomes if requested
if (!is.null(rm_chr)) {
  rm_chr_vec <- strsplit(rm_chr, ",")[[1]] |> trimws()
  if (any(rm_chr_vec %in% unique(fai$chr))) {
    fai <- fai %>% filter(!chr %in% rm_chr_vec)
  } else {
    stop("None of the chromosomes in --remove_chr are present in the FAI; check names.\n",
         call. = FALSE)
  }
} else {
  message("No chromosomes removed (API/mito will be included unless filtered earlier).")
}

fai$chr <- as.numeric(stringr::str_match(fai$chr, pattern)[, groupid])

# Define filename extension based on window / filters
if (min_n_snp == 0 & min_l_seg == 0) {
  filename_ext <- sprintf("win%dkb", window_size / 1000)
} else {
  filename_ext <- sprintf("win%dkb_minsnp%d_minlen%d",
                          window_size / 1000, min_n_snp, min_l_seg)
}

# ─────────────────────────────────────────────────────────────
# Summarise per-category hmmIBD results
# ─────────────────────────────────────────────────────────────
combined_ibd_r      <- data.frame()
combined_fraction_r <- data.frame()

maf_label <- as.character(th_maf)

for (category_n in category_list) {
  message("\nCategory: ", category_n)

  cat_dir <- file.path(workdir, category_n)

  frac_file <- file.path(cat_dir,
                         sprintf("hmmIBD_%s_maf%s_out.hmm_fract.txt",
                                 category_n, maf_label))
  ibd_file  <- file.path(cat_dir,
                         sprintf("hmmIBD_%s_maf%s_out.hmm.txt",
                                 category_n, maf_label))

  message("  Fraction file: ", frac_file)
  message("  Segments file: ", ibd_file)

  if (!file.exists(frac_file) || !file.exists(ibd_file)) {
    message("  -> hmmIBD results not found for ", category_n, "; skipping.")
    next
  }

  ibd_frac <- read_tsv(frac_file, col_types = cols()) %>% as.data.frame()
  ibd_data <- read_tsv(ibd_file,  col_types = cols()) %>% as.data.frame()

  if (nrow(ibd_frac) == 0 || nrow(ibd_data) == 0) {
    message("  -> hmmIBD results empty for ", category_n, "; skipping.")
    next
  }

  # Parse IBD segments
  ibd_data <- ibd_data %>%
    unite(id, sample1, sample2, sep = "_", remove = FALSE) %>%
    mutate(total = n_distinct(id),
           length = end - start + 1)

  ibd_conf <- ibd_data %>%
    filter(different == 0,
           Nsnp  > min_n_snp,
           length > min_l_seg) %>%
    arrange(chr, start, end) %>%
    rename(pos_start = start, pos_end = end)

  if (nrow(ibd_conf) == 0) {
    message("  -> No segments pass NSNP/LSEGMENT filters for ", category_n, "; skipping.")
    next
  }

  # Loop chromosomes
  chr_vec <- sort(unique(fai$chr))
  results <- data.frame()
  total_pairs <- unique(ibd_conf$total)

  for (k in seq_along(chr_vec)) {
    chr_k <- chr_vec[k]

    data_chr <- ibd_conf %>% filter(chr == chr_k)
    if (nrow(data_chr) == 0) next

    length_chr <- fai %>%
      filter(chr == chr_k) %>%
      pull(end_chr) %>%
      unique()

    # ensure a single scalar
    if (length(length_chr) != 1) {
      length_chr <- length_chr[1]
      warning("FAI has multiple entries for chr ", chr_k,
              "; using the first length = ", length_chr)
    }

    # Windows for this chromosome:
    # start at 1, step by window_size; end = start+window_size-1 capped at chr length
    win_start <- seq(1, length_chr, by = window_size)
    win_end   <- pmin(win_start + window_size - 1, length_chr)

    windows <- data.frame(
      chr       = chr_k,
      win_start = win_start,
      win_end   = win_end
    )

    # Overlap IBD segments with windows
    x <- data.table(
      start = as.numeric(as.character(data_chr$pos_start)),
      end   = as.numeric(as.character(data_chr$pos_end))
    )
    y <- data.table(
      start = as.numeric(as.character(windows$win_start)),
      end   = as.numeric(as.character(windows$win_end))
    )

    setkey(y, start, end)
    overlaps <- foverlaps(x, y, type = "any", which = TRUE)

    matched <- cbind(data_chr[overlaps$xid, ],
                     as.data.frame(y[overlaps$yid])) %>%
      mutate(rel_start = ifelse(pos_start > start, pos_start, start),
             rel_end   = ifelse(pos_end   < end,   pos_end,   end),
             rel_length = rel_end - rel_start + 1)

    matched_w <- windows %>%
      full_join(matched,
                by = c("chr", "win_start" = "start", "win_end" = "end"))

    matched_wm <- matched_w %>%
      group_by(win_start, win_end) %>%
      mutate(sum_seg_len = sum(rel_length, na.rm = TRUE),
             av_seg_len  = ifelse(!is.na(sum_seg_len),
                                  sum_seg_len / window_size, 0),
             fraction    = av_seg_len / total_pairs) %>%
      select(chr, win_start, win_end, sum_seg_len, av_seg_len, fraction) %>%
      distinct()

    results <- bind_rows(results, matched_wm)
  }

  if (nrow(results) == 0) {
    message("  -> No windowed results for ", category_n, "; skipping.")
    next
  }

  results_it <- results %>%
    as.data.frame() %>%
    mutate(category = category_n)

  combined_ibd_r      <- bind_rows(combined_ibd_r, results_it)

  frac_it <- ibd_frac %>%
    unite(id, sample1, sample2, sep = "_") %>%
    mutate(category = category_n) %>%
    select(id, category, fract_sites_IBD)

  combined_fraction_r <- bind_rows(combined_fraction_r, frac_it)
}

# ─────────────────────────────────────────────────────────────
# Save combined summaries
# ─────────────────────────────────────────────────────────────
if (nrow(combined_ibd_r) == 0 || nrow(combined_fraction_r) == 0) {
  stop("No summarised IBD or fraction data to save. Stopping.\n", call. = FALSE)
}

message("\nSaving combined summaries...")

ibd_out_file <- file.path(
  workdir,
  sprintf("%s_hmmIBD_ibd_%s.tsv", suffix, filename_ext)
)
write.table(combined_ibd_r, ibd_out_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

colnames(combined_fraction_r) <- c("id", "category", "fraction")
frac_out_file <- file.path(
  workdir,
  sprintf("%s_hmmIBD_fraction.tsv", suffix)
)
write.table(combined_fraction_r, frac_out_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

message("  -> ", ibd_out_file)
message("  -> ", frac_out_file)

# ─────────────────────────────────────────────────────────────
# Annotation step (genes / products)
# ─────────────────────────────────────────────────────────────
message("\nAnnotating high-IBD windows...")

if (is.null(gene_product_file) || !file.exists(gene_product_file)) {
  stop("Specify a valid --gene_product annotation file. Stopping.\n", call. = FALSE)
}

annotation <- readr::read_tsv(gene_product_file, col_types = cols())

# Remove chromosomes if requested
if (!is.null(rm_chr)) {
  rm_chr_vec <- strsplit(rm_chr, ",")[[1]] |> trimws()
  if (any(rm_chr_vec %in% unique(annotation$chr))) {
    annotation <- annotation %>% filter(!chr %in% rm_chr_vec)
  } else {
    stop("None of the chromosomes in --remove_chr were found in annotation; stopping.\n",
         call. = FALSE)
  }
} else {
  message("Annotation: no chromosomes removed (API/mito included if present).")
}

annotation$chr <- as.numeric(stringr::str_match(annotation$chr, pattern)[, groupid])

ibd_regions <- combined_ibd_r %>%
  rename(start = win_start, end = win_end) %>%
  select(chr, start, end) %>%
  distinct()

res_annot <- annotate_candidate_regions(ibd_regions, annotation)

# Per-category  quantiles
quantile_tbl <- combined_ibd_r %>%
  group_by(category) %>%
  mutate(qfrac = quantile(fraction, th_quantile, na.rm = TRUE)) %>%
  filter(fraction >= qfrac) %>%
  ungroup() %>%
  distinct()

quantile_annot <- quantile_tbl %>%
  rename(start = win_start, end = win_end) %>%
  inner_join(res_annot,
             by = c("chr", "start", "end")) %>%
  mutate(gene_product = gsub("\\t", "", gene_product)) %>%
  ungroup()

annot_file <- file.path(
  workdir,
  sprintf("%s_hmmIBD_ibd_%s_annotated_q%s.tsv",
          suffix, filename_ext, as.character(th_quantile))
)
write.table(quantile_annot, annot_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
message("  -> ", annot_file)

# Wide format
quantile_regs <- quantile_tbl %>%
  select(chr, win_start, win_end) %>%
  distinct()

quantile_annot_wide <- combined_ibd_r %>%
  right_join(quantile_regs,
             by = c("chr", "win_start", "win_end")) %>%
  rename(start = win_start, end = win_end) %>%
  left_join(res_annot,
            by = c("chr", "start", "end")) %>%
  select(-c(pos_start, pos_end)) %>%
  group_by(chr, start, end, category, fraction) %>%
  summarise(
    genes    = paste0(gene_id, "(", gene_name, ")", collapse = "; "),
    products = paste0(gene_product, collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    genes    = gsub("\\(NA\\)", "", genes),
    products = gsub("\\t", "", products),
    fraction = round(fraction, 3)
  ) %>%
  tidyr::pivot_wider(names_from = "category", values_from = "fraction") %>%
  ungroup()

annot_wide_file <- file.path(
  workdir,
  sprintf("%s_hmmIBD_ibd_%s_annotated_q%s_wide.tsv",
          suffix, filename_ext, as.character(th_quantile))
)
write.table(quantile_annot_wide, annot_wide_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

message("  -> ", annot_wide_file)
message("\nDone.")
