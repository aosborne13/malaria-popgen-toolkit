#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(data.table)
})

# -------------------------------------------------------------------
# Local helper: annotate_candidate_regions using data.table::foverlaps
# -------------------------------------------------------------------
annotate_candidate_regions <- function(regions, annotation) {
  # regions: data.frame with chr, start, end
  # annotation: data.frame with chr, pos_start, pos_end, gene_id, gene_name, gene_product
  reg <- as.data.table(regions)
  ann <- as.data.table(annotation)

  setkey(ann, chr, pos_start, pos_end)
  setkey(reg, chr, start, end)

  ov <- foverlaps(
    reg,
    ann,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "pos_start", "pos_end"),
    type = "any",
    nomatch = 0L
  )

  as.data.frame(ov)
}

# -------------------------------------------------------------------
# Options
# -------------------------------------------------------------------
option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = NULL,
              help = "Main directory (where per-category hmmIBD outputs live)",
              metavar = "character"),
  make_option("--list_category", type = "character", default = NULL,
              help = "Optional file with category list (one per line). "
                   "If omitted, categories are auto-detected from hmmIBD_* files.",
              metavar = "character"),
  make_option(c("-l", "--legend"), type = "character",
              default = "ibd_matrix_hap_leg.tsv",
              help = "SNP legend file [default %default]",
              metavar = "character"),
  make_option("--gene_product", type = "character", default = NULL,
              help = "Gene product annotation TSV",
              metavar = "character"),
  make_option(c("-r", "--ref_index"), type = "character", default = NULL,
              help = "Reference FASTA index (.fai)",
              metavar = "character"),
  make_option("--maf", type = "numeric", default = 0.01,
              help = "MAF used when generating hmmIBD inputs [default %default]",
              metavar = "numeric"),
  make_option(c("-w", "--window_size"), type = "numeric", default = 10000,
              help = "Window size in bp [default %default]",
              metavar = "numeric"),
  make_option("--quantile_cutoff", type = "numeric", default = 0.95,
              help = "Quantile cutoff for high-IBD windows [default %default]",
              metavar = "numeric"),
  make_option("--suffix", type = "character",
              default = format(Sys.time(), "%d_%m_%Y"),
              help = "Suffix for output files [default current date]",
              metavar = "character"),
  make_option(c("-t", "--threads"), type = "integer", default = 4,
              help = "Threads for data.table [default %default]",
              metavar = "numeric"),
  make_option("--remove_chr", type = "character", default = NULL,
              help = "Comma-separated chromosomes to remove, e.g. Pf3D7_API_v3,Pf3D7_MIT_v3",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex to parse chromosome names from FAI [default %default]",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Capture group index for numeric chromosome [default %default]",
              metavar = "numeric"),
  make_option("--NSNP", type = "numeric", default = 0,
              help = "Min number of SNPs per IBD segment [default %default]",
              metavar = "numeric"),
  make_option("--LSEGMENT", type = "numeric", default = 0,
              help = "Min length (bp) per IBD segment [default %default]",
              metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

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

if (is.null(workdir)) stop("ERROR: --workdir is required", call. = FALSE)
if (is.null(ref_index)) stop("ERROR: --ref_index is required", call. = FALSE)
if (is.null(gene_product_file)) stop("ERROR: --gene_product is required", call. = FALSE)

cat("Workdir:", workdir, "\n")

setDTthreads(threads)

# -------------------------------------------------------------------
# Category list (auto-detect if not provided)
# -------------------------------------------------------------------
if (!is.null(category_list) && file.exists(category_list)) {
  categories <- read.table(category_list, header = FALSE, stringsAsFactors = FALSE)$V1
} else {
  # Auto-detect from per-category hmmIBD outputs in subdirs:
  # e.g. workdir/<cat>/hmmIBD_<cat>_maf0.01_out.hmm.txt
  subdirs <- list.dirs(workdir, full.names = FALSE, recursive = FALSE)
  candidates <- character(0)
  for (cat_name in subdirs) {
    seg_file <- file.path(
      workdir,
      cat_name,
      sprintf("hmmIBD_%s_maf%s_out.hmm.txt", cat_name, as.character(th_maf))
    )
    if (file.exists(seg_file)) {
      candidates <- c(candidates, cat_name)
    }
  }
  if (length(candidates) == 0) {
    stop("ERROR: could not auto-detect any categories in ", workdir,
         ". Provide --list_category or ensure per-category hmmIBD outputs exist.",
         call. = FALSE)
  }
  categories <- sort(unique(candidates))
}
cat("Auto-detected categories:", paste(categories, collapse = ", "), "\n")

# -------------------------------------------------------------------
# Legend
# -------------------------------------------------------------------
legend_path <- file.path(workdir, legend)
if (!file.exists(legend_path)) {
  stop("Cannot locate SNP legend file: ", legend_path, call. = FALSE)
}
cat("Loaded legend:", legend_path, "\n")

snp_hmmibd_02_1 <- read_tsv(legend_path, col_types = cols()) %>% as.data.frame()

# -------------------------------------------------------------------
# Reference index / FAI
# -------------------------------------------------------------------
if (!file.exists(ref_index)) {
  stop("Cannot locate reference index: ", ref_index, call. = FALSE)
}

fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  rename(chr = V1, end_chr = V2) %>%
  select(chr, end_chr) %>%
  mutate(start_chr = 1)

# Remove unwanted chromosomes
if (!is.null(rm_chr)) {
  rm_chr_ <- strsplit(rm_chr, ",")[[1]] |> trimws()
  if (any(rm_chr_ %in% unique(fai$chr))) {
    rm_chr_ <- rm_chr_[rm_chr_ %in% unique(fai$chr)]
    fai <- fai %>% filter(!chr %in% rm_chr_)
  } else {
    stop("Wrong name for chromosomes to remove. Stopping...")
  }
} else {
  message("No chromosomes removed. API/mito should usually be removed!")
}

# Convert chr string -> numeric (Pf3D7_01_v3 -> 1, etc.)
fai$chr <- as.numeric(stringr::str_match(fai$chr, pattern)[, groupid])

# -------------------------------------------------------------------
# Output filename suffix for window size / filters
# -------------------------------------------------------------------
if (min_n_snp == 0 & min_l_seg == 0) {
  filename_ext <- sprintf("win%dkb", window_size / 1000)
} else {
  filename_ext <- sprintf("win%dkb_minsnp%d_minlen%d",
                          window_size / 1000, min_n_snp, min_l_seg)
}

combined_ibd_r <- data.frame()
combined_fraction_r <- data.frame()

# -------------------------------------------------------------------
# Per-category processing
# -------------------------------------------------------------------
for (category_n in categories) {

  message("\nCategory: ", category_n)

  frac_file <- file.path(
    workdir,
    category_n,
    sprintf("hmmIBD_%s_maf%s_out.hmm_fract.txt", category_n, as.character(th_maf))
  )
  seg_file <- file.path(
    workdir,
    category_n,
    sprintf("hmmIBD_%s_maf%s_out.hmm.txt", category_n, as.character(th_maf))
  )

  message("  Fraction file: ", frac_file)
  message("  Segments file: ", seg_file)

  if (!file.exists(frac_file) || !file.exists(seg_file)) {
    message("  -> hmmIBD results not found for ", category_n, "; skipping.")
    next
  }

  ibd_frac <- read_tsv(frac_file, col_types = cols()) %>% as.data.frame()
  ibd_data <- read_tsv(seg_file, col_types = cols()) %>% as.data.frame()

  if (nrow(ibd_frac) == 0 || nrow(ibd_data) == 0) {
    message("  -> hmmIBD results empty for ", category_n, "; skipping.")
    next
  }

  # Parse IBD segments
  ibd_data <- ibd_data %>%
    unite(id, sample1, sample2, sep = "_", remove = FALSE) %>%
    mutate(
      total  = n_distinct(id),
      length = end - start + 1
    )

  ibd_conf <- ibd_data %>%
    filter(different == 0, Nsnp > min_n_snp, length > min_l_seg) %>%
    arrange(chr, start, end) %>%
    rename(pos_start = start, pos_end = end)

  if (nrow(ibd_conf) == 0) {
    message("  -> no segments pass filters for ", category_n, "; skipping.")
    next
  }

  # Build windows per chr, compute fraction IBD per window
  chr_vec <- sort(unique(ibd_conf$chr))
  results <- data.frame()

  for (chr_k in chr_vec) {
    data_chr <- ibd_conf %>% filter(chr == chr_k)
    length_chr <- fai %>%
      filter(chr == chr_k) %>%
      pull(end_chr) %>%
      unique()

    if (length(length_chr) != 1) {
      length_chr <- length_chr[1]
      warning("FAI has multiple entries for chr ", chr_k,
              "; using first length = ", length_chr)
    }

    # Window boundaries
    win_start <- seq(1, length_chr, by = window_size)
    win_end   <- pmin(win_start + window_size - 1, length_chr)

    windows <- data.frame(
      chr       = chr_k,
      win_start = win_start,
      win_end   = win_end
    )

    # Overlap IBD segments with windows
    x <- data.table(
      start = as.numeric(data_chr$pos_start),
      end   = as.numeric(data_chr$pos_end)
    )
    y <- data.table(
      start = as.numeric(windows$win_start),
      end   = as.numeric(windows$win_end)
    )

    setkey(y, start, end)
    overlaps <- foverlaps(x, y, type = "any", which = TRUE)

    matched <- cbind(data_chr[overlaps$xid, ], as.data.frame(y[overlaps$yid])) %>%
      mutate(
        rel_start  = ifelse(pos_start > start, pos_start, start),
        rel_end    = ifelse(pos_end   < end,   pos_end,   end),
        rel_length = rel_end - rel_start + 1
      )

    matched_w <- windows %>%
      full_join(matched, by = c("chr", "win_start" = "start", "win_end" = "end"))

    total_pairs <- unique(ibd_conf$total)

    matched_wm <- matched_w %>%
      group_by(win_start, win_end) %>%
      mutate(
        sum_seg_len = sum(rel_length, na.rm = TRUE),
        av_seg_len  = ifelse(!is.na(sum_seg_len), sum_seg_len / window_size, 0),
        fraction    = av_seg_len / total_pairs
      ) %>%
      select(chr, win_start, win_end, sum_seg_len, av_seg_len, fraction) %>%
      distinct()

    results <- bind_rows(results, matched_wm)
  }

  results_it <- results %>%
    mutate(category = category_n)

  combined_ibd_r <- bind_rows(combined_ibd_r, results_it)

  frac_it <- ibd_frac %>%
    unite(id, sample1, sample2, sep = "_") %>%
    mutate(category = category_n) %>%
    select(id, category, fract_sites_IBD)

  combined_fraction_r <- bind_rows(combined_fraction_r, frac_it)
}

# -------------------------------------------------------------------
# Save combined results
# -------------------------------------------------------------------
message("\nSaving combined summaries...")

if (nrow(combined_ibd_r) == 0 || nrow(combined_fraction_r) == 0) {
  stop("None results to save. Stopping...", call. = FALSE)
}

ibd_out <- file.path(workdir,
                     sprintf("%s_hmmIBD_ibd_%s.tsv", suffix, filename_ext))
frac_out <- file.path(workdir,
                      sprintf("%s_hmmIBD_fraction.tsv", suffix))

write.table(
  combined_ibd_r,
  ibd_out,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

colnames(combined_fraction_r) <- c("id", "category", "fraction")
write.table(
  combined_fraction_r,
  frac_out,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

cat("  -> ", ibd_out, "\n")
cat("  -> ", frac_out, "\n")

# -------------------------------------------------------------------
# Annotation of high-IBD windows
# -------------------------------------------------------------------
message("\nAnnotating high-IBD windows...")

if (!file.exists(gene_product_file)) {
  stop("Annotation file not found: ", gene_product_file, call. = FALSE)
}

annotation <- readr::read_tsv(gene_product_file, col_types = cols())

# Remove chromosomes
if (!is.null(rm_chr)) {
  rm_chr_ <- strsplit(rm_chr, ",")[[1]] |> trimws()
  if (any(rm_chr_ %in% unique(annotation$chr))) {
    rm_chr_ <- rm_chr_[rm_chr_ %in% unique(annotation$chr)]
    annotation <- annotation %>% filter(!chr %in% rm_chr_)
  } else {
    stop("Wrong name for chromosomes to remove (annotation). Stopping...")
  }
} else {
  message("No chromosomes removed from annotation (API/mito still present).")
}

annotation$chr <- as.numeric(stringr::str_match(annotation$chr, pattern)[, groupid])

ibd_regions <- combined_ibd_r %>%
  rename(start = win_start, end = win_end) %>%
  select(chr, start, end) %>%
  distinct()

res_annot <- annotate_candidate_regions(ibd_regions, annotation)

quantile <- combined_ibd_r %>%
  group_by(category) %>%
  mutate(qfrac = quantile(fraction, th_quantile, na.rm = TRUE)) %>%
  filter(fraction >= qfrac) %>%
  ungroup() %>%
  distinct()

quantile_annot <- quantile %>%
  rename(start = win_start, end = win_end) %>%
  inner_join(res_annot, by = c("chr", "start", "end")) %>%
  mutate(gene_product = gsub("\\t", "", gene_product))

annot_out <- file.path(
  workdir,
  sprintf("%s_hmmIBD_ibd_%s_annotated_q%s.tsv",
          suffix, filename_ext, as.character(th_quantile))
)
write.table(
  quantile_annot,
  annot_out,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

# Wide format
quantile_regs <- quantile %>%
  select(chr, win_start, win_end) %>%
  distinct()

quantile_annot_wide <- combined_ibd_r %>%
  right_join(quantile_regs,
             by = c("chr", "win_start", "win_end")) %>%
  rename(start = win_start, end = win_end) %>%
  left_join(res_annot, by = c("chr", "start", "end")) %>%
  select(-c(pos_start, pos_end)) %>%
  group_by(chr, start, end, category, fraction) %>%
  summarize(
    genes = paste0(gene_id, "(", gene_name, ")", collapse = "; "),
    products = paste0(gene_product, collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    genes = gsub("\\(NA\\)", "", genes),
    products = gsub("\\t", "", products),
    fraction = round(fraction, 3)
  ) %>%
  tidyr::pivot_wider(
    names_from  = "category",
    values_from = "fraction"
  )

annot_wide_out <- file.path(
  workdir,
  sprintf("%s_hmmIBD_ibd_%s_annotated_q%s_wide.tsv",
          suffix, filename_ext, as.character(th_quantile))
)

write.table(
  quantile_annot_wide,
  annot_wide_out,
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

cat("  -> ", annot_out, "\n")
cat("  -> ", annot_wide_out, "\n")

message("\nDone.")

