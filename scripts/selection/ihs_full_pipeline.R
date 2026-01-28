#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(rehh)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
})

# ─────────────────────────────────────────────────────────────
# Helper: MAF calculation (0/1/NA matrix)
# ─────────────────────────────────────────────────────────────
calculate_maf <- function(mat) {
  m <- as.matrix(mat)
  m[m < 0] <- NA
  p <- rowMeans(m == 1, na.rm = TRUE)  # alt allele freq (assuming 1 encodes ALT)
  p[is.nan(p)] <- NA_real_
  maf <- pmin(p, 1 - p)
  maf
}

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Working/output directory [default %default]",
              metavar = "character"),
  make_option(c("-b", "--matrix_binary"), type = "character", default = NULL,
              help = "Input filename of filtered binary matrix [required]",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Full path to metadata TSV file [required]",
              metavar = "character"),

  # Primary grouping (e.g. country) and optional subgroup (e.g. year)
  make_option(c("-c", "--category"), type = "character", default = NULL,
              help = "Primary category value(s) in label_category column, comma-separated, or 'ALL' for all.",
              metavar = "character"),
  make_option(c("--label_category"), type = "character", default = "country",
              help = "Column name in metadata for primary category (e.g. country) [default %default]",
              metavar = "character"),
  make_option(c("--subgroup_col"), type = "character", default = NULL,
              help = "Optional subgroup column in metadata (e.g. year); if set, iHS is run per (category, subgroup) combo.",
              metavar = "character"),

  make_option(c("--label_id"), type = "character", default = "sample_id",
              help = "Column name in metadata for sample ID [default %default]",
              metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = "fws",
              help = "Column name in metadata for Fws values [default %default]",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold (keep samples with Fws >= this) [default %default]",
              metavar = "number"),
  make_option(c("--maf"), type = "numeric", default = 0.01,
              help = "MAF threshold for SNPs [default %default]",
              metavar = "number"),
  make_option(c("--rehh_min_perc_hap"), type = "numeric", default = 80,
              help = "rehh::data2haplohh min_perc_geno.hap [default %default]",
              metavar = "number"),
  make_option(c("--rehh_min_perc_mrk"), type = "numeric", default = 70,
              help = "rehh::data2haplohh min_perc_geno.mrk [default %default]",
              metavar = "number"),
  make_option(c("--na_char"), type = "character", default = "NA",
              help = "Missing genotype character (also treats '.' as missing) [default %default]",
              metavar = "character"),
  make_option("--forced_recode", action = "store_true", default = FALSE,
              help = "Recode NA -> REF (0) and mixed (0.5) -> ALT (1)"),
  make_option("--forced_mixed", action = "store_true", default = FALSE,
              help = "Recode mixed (0.5) -> ALT (1), keep NA as NA"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Chromosomes to remove, e.g. Pf3D7_API_v3,Pf3D7_MIT_v3",
              metavar = "character"),
  make_option("--regex_chr", type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection (default matches Pf3D7_01_v3)",
              metavar = "character"),
  make_option("--regex_groupid", type = "numeric", default = 3,
              help = "Capture group index in regex_chr giving numeric chromosome [default %default]",
              metavar = "numeric"),
  make_option(c("--threads"), type = "integer", default = 4,
              help = "Number of threads for data.table [default %default]",
              metavar = "number"),

  # Gene-level annotation for iHS (only genome product TSV now)
  make_option(c("--genome-file"), type = "character", default = NULL,
              help = "Genome product TSV with columns: chr,pos_start,pos_end,gene_id,gene_product,gene_name [required]",
              metavar = "character"),
  make_option(c("--focus-pop"), type = "character", default = NULL,
              help = "Focus population label for detailed plots (must match one of the resulting categories, e.g. Ethiopia_2017)",
              metavar = "character"),
  make_option(c("--min-maf-ihs"), type = "numeric", default = 0.0,
              help = "min_maf argument to ihh2ihs [default %default]",
              metavar = "numeric"),
  make_option(c("--freqbin"), type = "numeric", default = 0.05,
              help = "freqbin argument to ihh2ihs [default %default]",
              metavar = "numeric"),
  make_option(c("--ihs-thresh"), type = "numeric", default = 2.5,
              help = "Absolute iHS threshold for labeling (|IHS| > this) [default %default]",
              metavar = "numeric"),
  make_option(c("--logp-thresh"), type = "numeric", default = 5.0,
              help = "-log10(p) threshold for labeling [default %default]",
              metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ─────────────────────────────────────────────────────────────
# Extract options
# ─────────────────────────────────────────────────────────────
workdir         <- opt$workdir
bin_mat_file    <- opt$matrix_binary
metadata_file   <- opt$metadata
category_arg    <- opt$category
label_category  <- opt$label_category
subgroup_col    <- opt$subgroup_col
label_id        <- opt$label_id
label_fws       <- opt$label_fws
threshold_fws   <- opt$fws_th
th_maf          <- opt$maf
th_min_perc_sam <- opt$rehh_min_perc_hap
th_min_perc_snp <- opt$rehh_min_perc_mrk
na_char         <- opt$na_char
forced_recode   <- isTRUE(opt$forced_recode)
forced_mixed    <- isTRUE(opt$forced_mixed)
rm_chr_str      <- opt$remove_chr
pattern         <- opt$regex_chr
groupid         <- opt$regex_groupid
threads         <- opt$threads

genome_file  <- opt$`genome-file`
focus_pop    <- opt$`focus-pop`
min_maf_ihs  <- opt$`min-maf-ihs`
freqbin      <- opt$freqbin
ihs_thresh   <- opt$`ihs-thresh`
logp_thresh  <- opt$`logp-thresh`

# ─────────────────────────────────────────────────────────────
# Validate inputs
# ─────────────────────────────────────────────────────────────
if (is.null(bin_mat_file)) stop("ERROR: --matrix_binary is required\n", call. = FALSE)
if (!file.exists(bin_mat_file)) stop("ERROR: matrix_binary file not found: ", bin_mat_file, "\n", call. = FALSE)

if (is.null(metadata_file)) stop("ERROR: --metadata is required\n", call. = FALSE)
if (!file.exists(metadata_file)) stop("ERROR: metadata file not found: ", metadata_file, "\n", call. = FALSE)

if (is.null(genome_file) || !file.exists(genome_file)) {
  stop("ERROR: --genome-file must be provided and must exist.\n", call. = FALSE)
}

dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
setDTthreads(threads)

message("Workdir:        ", workdir)
message("Binary matrix:  ", bin_mat_file)
message("Metadata:       ", metadata_file)
message("Genome file:    ", genome_file)

# ─────────────────────────────────────────────────────────────
# Load metadata
# ─────────────────────────────────────────────────────────────
metadata_all <- read_tsv(metadata_file, show_col_types = FALSE)

needed_cols <- c(label_id, label_category, label_fws)
if (!all(needed_cols %in% colnames(metadata_all))) {
  stop("ERROR: Missing columns in metadata: ",
       paste(setdiff(needed_cols, colnames(metadata_all)), collapse = ", "),
       "\n", call. = FALSE)
}

metadata_all[[label_category]] <- as.character(metadata_all[[label_category]])
metadata_all[[label_id]]       <- as.character(metadata_all[[label_id]])

if (!is.null(subgroup_col) && nzchar(subgroup_col)) {
  if (!subgroup_col %in% colnames(metadata_all)) {
    stop("ERROR: subgroup_col '", subgroup_col,
         "' not found in metadata.\n", call. = FALSE)
  }
}

# ─────────────────────────────────────────────────────────────
# Build grouping definition table (cat_def)
# cat_def has columns: parent_cat, subgroup (may be NA), cat_label
# ─────────────────────────────────────────────────────────────
if (is.null(subgroup_col) || !nzchar(subgroup_col)) {
  available_categories <- metadata_all %>%
    dplyr::select(all_of(label_category)) %>%
    distinct() %>%
    pull() %>%
    sort()

  if (length(available_categories) == 0) {
    stop("ERROR: No categories found in metadata column '", label_category, "'.\n",
         call. = FALSE)
  }

  message("Available categories in ", label_category, ": ",
          paste(available_categories, collapse = ", "))

  if (is.null(category_arg) || category_arg == "" || tolower(category_arg) == "all") {
    categories_vec <- available_categories
    message("No specific category provided; using ALL categories.")
  } else {
    categories_vec <- strsplit(category_arg, ",")[[1]] |> trimws()
    missing_cats <- setdiff(categories_vec, available_categories)
    if (length(missing_cats) > 0) {
      stop("ERROR: Requested category(ies) not found in metadata: ",
           paste(missing_cats, collapse = ", "), "\n",
           "Available: ", paste(available_categories, collapse = ", "), "\n",
           call. = FALSE)
    }
  }

  cat_def <- tibble(
    parent_cat = categories_vec,
    subgroup   = NA_character_,
    cat_label  = categories_vec
  )
  message("Will run iHS for categories: ", paste(cat_def$cat_label, collapse = ", "))

} else {
  combo_tbl <- metadata_all %>%
    filter(!is.na(.data[[label_category]]), !is.na(.data[[subgroup_col]])) %>%
    distinct(parent_cat = .data[[label_category]],
             subgroup   = .data[[subgroup_col]])

  if (nrow(combo_tbl) == 0) {
    stop("ERROR: No (", label_category, ", ", subgroup_col,
         ") combinations found in metadata.\n", call. = FALSE)
  }

  message("Available primary categories in ", label_category, ": ",
          paste(sort(unique(combo_tbl$parent_cat)), collapse = ", "))

  if (is.null(category_arg) || category_arg == "" || tolower(category_arg) == "all") {
    requested_parents <- sort(unique(combo_tbl$parent_cat))
    message("No specific primary category provided; using ALL: ",
            paste(requested_parents, collapse = ", "))
  } else {
    requested_parents <- strsplit(category_arg, ",")[[1]] |> trimws()
    missing_parents <- setdiff(requested_parents, combo_tbl$parent_cat)
    if (length(missing_parents) > 0) {
      stop("ERROR: Requested primary category(ies) not found in metadata: ",
           paste(missing_parents, collapse = ", "), "\n",
           "Available: ", paste(sort(unique(combo_tbl$parent_cat)), collapse = ", "), "\n",
           call. = FALSE)
    }
  }

  cat_def <- combo_tbl %>%
    filter(parent_cat %in% requested_parents) %>%
    mutate(
      subgroup  = as.character(subgroup),
      cat_label = paste(parent_cat, subgroup, sep = "_")
    ) %>%
    arrange(parent_cat, subgroup)

  if (nrow(cat_def) == 0) {
    stop("ERROR: No (category, subgroup) combinations left after filtering.\n", call. = FALSE)
  }

  message("Will run iHS for groups: ", paste(cat_def$cat_label, collapse = ", "))
}

# Decide focus population label (must match cat_label)
if (is.null(focus_pop) || focus_pop == "" || !(focus_pop %in% cat_def$cat_label)) {
  if (!is.null(focus_pop) && focus_pop != "" && !(focus_pop %in% cat_def$cat_label)) {
    warning("Requested focus-pop '", focus_pop, "' not among resulting groups; defaulting to first.\n")
  }
  focus_pop <- cat_def$cat_label[1]
}
message("Focus population: ", focus_pop)

# ─────────────────────────────────────────────────────────────
# Define drug-resistance gene regions (coordinates are numeric CHR/POS)
# ─────────────────────────────────────────────────────────────
drug_genes <- tibble::tribble(
  ~gene,   ~CHR, ~START,   ~END,
  "CRT",        7,  403000,  406000,
  "K13",       13, 1724817, 1728000,
  "MDR1",       5,  957500,  962000,
  "DHFR",       4,  746000,  751000,
  "DHPS",       8,  548000,  552000,
  "PX1",        7,  889800,  899213,
  "UBP1",       1,  188400,  201400,
  "AP2-MU",    12,  716700,  720700,
  "AAT1",       6, 1213100, 1217350
)

# ─────────────────────────────────────────────────────────────
# Load genome product annotation (intervals)
# Expected headers: chr pos_start pos_end gene_id gene_product gene_name
# ─────────────────────────────────────────────────────────────
genes <- read_tsv(genome_file, show_col_types = FALSE)

genes_clean <- genes %>%
  mutate(CHR = as.numeric(str_extract(chr, "\\d{2}"))) %>%
  rename(
    START     = pos_start,
    END       = pos_end,
    GENE_ID   = gene_id,
    GENE_NAME = gene_name,
    PRODUCT   = gene_product
  ) %>%
  filter(!is.na(CHR)) %>%
  select(CHR, START, END, GENE_ID, GENE_NAME, PRODUCT)

# ─────────────────────────────────────────────────────────────
# Main loop: for each (parent_cat, subgroup, cat_label)
# ─────────────────────────────────────────────────────────────
all_ihs <- list()

for (i_grp in seq_len(nrow(cat_def))) {
  parent_cat   <- cat_def$parent_cat[i_grp]
  subgroup_val <- cat_def$subgroup[i_grp]
  cat_label    <- cat_def$cat_label[i_grp]

  message("\n==============================")
  message("Category group: ", cat_label)
  message("==============================")

  if (is.na(subgroup_val)) {
    meta_cat <- metadata_all %>%
      filter(.data[[label_category]] == parent_cat)
  } else {
    meta_cat <- metadata_all %>%
      filter(.data[[label_category]] == parent_cat,
             .data[[subgroup_col]] == subgroup_val)
  }

  if (nrow(meta_cat) == 0) {
    warning("  Group ", cat_label, " has no metadata rows; skipping.")
    next
  }

  category_str <- gsub(" ", "_", cat_label)

  # Columns to select from binary matrix: chr,pos,ref + samples in this group
  sample_ids <- meta_cat %>%
    dplyr::select(all_of(label_id)) %>%
    pull() %>%
    unique()

  if (length(sample_ids) == 0) {
    warning("  Group ", cat_label, ": no sample IDs in metadata; skipping.")
    next
  }

  samples_select <- c("chr", "pos", "ref", sample_ids)

  message("  Loading SNP matrix subset for group (", length(sample_ids), " samples)...")

  snp <- fread(bin_mat_file, sep = "\t",
               select = samples_select,
               header = TRUE, data.table = FALSE)

  if (!all(c("chr", "pos", "ref") %in% colnames(snp)[1:3])) {
    stop("ERROR: Matrix must start with columns: chr, pos, ref\n", call. = FALSE)
  }

  # Filter chromosomes if requested
  if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
    rm_chr <- strsplit(rm_chr_str, ",")[[1]] |> trimws()
    if (any(rm_chr %in% unique(snp$chr))) {
      message("  Removing chromosomes: ", paste(rm_chr, collapse = ", "))
      snp <- snp %>% filter(!chr %in% rm_chr)
    } else {
      warning("  None of the chromosomes in --remove_chr found in matrix for this group; continuing.")
    }
  } else {
    message("  No chromosomes removed (API/mito may still be present).")
  }

  # Transform chromosome from string to numeric via regex
  snp$chr <- as.numeric(stringr::str_match(snp$chr, pattern)[, groupid])

  # Keep metadata rows whose IDs are present in matrix columns
  meta_cat <- meta_cat %>%
    filter(.data[[label_id]] %in% colnames(snp)[-(1:3)])
  if (nrow(meta_cat) < 2) {
    warning("  Group ", cat_label, ": fewer than 2 samples in matrix; skipping.")
    next
  }

  # Recode missing data
  if (!is.na(na_char)) {
    snp[snp == na_char] <- NA
    snp[snp == "."]     <- NA
  }

  # Split genotype and SNP description
  snp_c <- as.data.frame(snp[, -(1:3)])
  snp_c[] <- lapply(snp_c, function(x) suppressWarnings(as.numeric(x)))
  snp_d <- as.data.frame(snp[, 1:3])

  rm(snp)

  # Recode mixed/missing according to flags
  maj3 <- snp_c
  if (forced_recode) {
    maj3[is.na(maj3)] <- 0
    maj3[maj3 == 0.5] <- 1
  } else if (forced_mixed) {
    maj3[maj3 == 0.5] <- 1
  } else {
    maj3[maj3 == 0.5] <- NA
  }
  rm(snp_c)

  # MAF filter
  maf_vals <- calculate_maf(maj3)
  to_keep_sti <- which(!is.na(maf_vals) & maf_vals >= th_maf)
  rm(maf_vals)

  if (length(to_keep_sti) == 0) {
    warning("  Group ", cat_label, ": no SNPs pass MAF ≥ ", th_maf, "; skipping.")
    next
  }

  snp_d_filt <- snp_d[to_keep_sti, ]
  maj4       <- maj3[to_keep_sti, , drop = FALSE]
  rm(maj3, snp_d)

  # Apply Fws filter to samples
  samples_keep <- meta_cat %>%
    filter(.data[[label_fws]] >= threshold_fws) %>%
    pull(all_of(label_id)) %>%
    unique()

  if (length(samples_keep) < 2) {
    warning("  Group ", cat_label, ": fewer than 2 samples pass Fws ≥ ", threshold_fws, "; skipping.")
    next
  }

  maj4 <- maj4[, samples_keep, drop = FALSE]

  # ─────────────────────────────────────────────────────────────
  # Build SNP map WITHOUT external SNP annotation file
  # We only need a consistent marker ID + chr + pos.
  # ref comes from the matrix. alt is a placeholder (not used for iHS stats).
  # ─────────────────────────────────────────────────────────────
  snp_d_filt$chr <- as.numeric(snp_d_filt$chr)
  snp_d_filt$pos <- suppressWarnings(as.numeric(snp_d_filt$pos))

  map <- snp_d_filt %>%
    mutate(
      info = paste(chr, pos, sep = "_"),
      alt  = "ALT"
    ) %>%
    select(info, chr, pos, ref, alt) %>%
    distinct()

  # Make haplotypes (whole-genome matrix; we’ll subset per chr later)
  hap <- t(maj4)
  hap_c <- hap
  hap_c[hap == 1]   <- 2
  hap_c[hap == 0]   <- 1
  hap_c[is.na(hap)] <- 0
  colnames(hap_c) <- NULL
  rownames(hap_c) <- NULL

  rm(maj4, hap)

  hap_c <- as.data.frame(hap_c, stringsAsFactors = FALSE)
  is_char <- sapply(hap_c, is.character)
  hap_c[is_char] <- lapply(hap_c[is_char], function(x) suppressWarnings(as.numeric(x)))

  # Create haplotypes per chromosome and run scan_hh
  u_chr <- sort(unique(map$chr))
  results_hh <- NULL

  log_file <- file.path(workdir, paste0(category_str, "_rehh.log"))
  sink(log_file, append = FALSE)
  cat("## Create haplotypes: data2haplohh() & scan_hh() for group ", category_str, " ##\n", sep = "")

  for (uchr in u_chr) {
    cat("Chr ", uchr, " ...\n", sep = "")

    idx_chr <- which(map$chr == uchr)
    if (length(idx_chr) == 0) {
      cat("  No markers for chr ", uchr, " after filtering; skipping.\n", sep = "")
      next
    }

    hap_chr_s <- hap_c[, idx_chr, drop = FALSE]
    map_chr   <- map[idx_chr, , drop = FALSE]

    if (ncol(hap_chr_s) != nrow(map_chr)) {
      stop("For group ", category_str, " chr ", uchr,
           ": number of markers in hap (", ncol(hap_chr_s),
           ") != markers in map (", nrow(map_chr), ").", call. = FALSE)
    }

    hap_ids <- paste0("h", seq_len(nrow(hap_chr_s)))
    hap_chr_out <- cbind(hap_ids, hap_chr_s)

    hap_file <- file.path(workdir, sprintf("hap_chr%d_%s", uchr, category_str))
    write.table(hap_chr_out, hap_file,
                sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

    map_file_chr <- file.path(workdir, sprintf("snp.info.chr%d.%s", uchr, category_str))
    write.table(map_chr, map_file_chr,
                quote = FALSE, row.names = FALSE, col.names = FALSE)

    hap_chr_pop <- data2haplohh(
      hap_file          = hap_file,
      map_file          = map_file_chr,
      recode.allele     = FALSE,
      chr.name          = uchr,
      min_perc_geno.hap = th_min_perc_sam,
      min_perc_geno.mrk = th_min_perc_snp,
      min_maf           = 0
    )

    res_chr_s <- scan_hh(hap_chr_pop)
    results_hh <- rbind(results_hh, res_chr_s)
    cat("\n")
  }
  sink()

  if (is.null(results_hh) || nrow(results_hh) == 0) {
    warning("  scan_hh produced no results for group ", cat_label, "; skipping.")
    next
  }

  scanned_file <- file.path(workdir, sprintf("scanned_haplotypes_%s.tsv", category_str))
  results_hh$category <- category_str
  write.table(results_hh, scanned_file,
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  message("  Wrote scan_hh results: ", scanned_file)

  message("  Computing iHS via ihh2ihs()...")
  ihs_obj <- ihh2ihs(results_hh, min_maf = min_maf_ihs, freqbin = freqbin)

  if (!is.null(ihs_obj$ihs)) {
    ihs_tbl <- as_tibble(ihs_obj$ihs) %>%
      mutate(category_name = cat_label)
    all_ihs[[cat_label]] <- ihs_tbl
    message("  iHS markers: ", nrow(ihs_tbl))
  } else {
    warning("  No valid iHS markers for group: ", cat_label)
  }
}

if (length(all_ihs) == 0) {
  stop("No iHS data could be computed for any group.\n", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────
# Combine iHS for all groups, coerce CHR/POSITION to numeric, annotate extremes
# ─────────────────────────────────────────────────────────────
ihs_all_countries <- bind_rows(all_ihs) %>%
  mutate(
    CHR      = as.numeric(CHR),
    POSITION = as.numeric(POSITION)
  )

ihs_file <- file.path(workdir, "iHS_all_countries.tsv")
write_tsv(ihs_all_countries, ihs_file)
message("Wrote combined iHS file: ", ihs_file)

ihs_all <- ihs_all_countries %>%
  filter(CHR %in% 1:14)

ihs_extreme <- ihs_all %>%
  filter(abs(IHS) > ihs_thresh)

ihs_annotated <- ihs_extreme %>%
  inner_join(genes_clean, by = "CHR") %>%
  filter(POSITION >= START, POSITION <= END) %>%
  select(category_name, CHR, POSITION, IHS, GENE_ID, GENE_NAME, PRODUCT)

annot_file <- file.path(workdir, "iHS_extreme_sites_annotated.tsv")
write_tsv(ihs_annotated, annot_file)
message("Wrote extreme iHS annotation: ", annot_file)

# ─────────────────────────────────────────────────────────────
# Per-group Manhattan plots (all groups together)
# ─────────────────────────────────────────────────────────────
plot_dir <- file.path(workdir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

ihs_all <- ihs_all %>%
  rowwise() %>%
  mutate(
    drug_gene = drug_genes$gene[which(
      CHR == drug_genes$CHR &
        POSITION >= drug_genes$START &
        POSITION <= drug_genes$END
    )[1]]
  ) %>%
  ungroup() %>%
  mutate(
    drug_gene   = ifelse(is.na(drug_gene), "None", drug_gene),
    is_drug_snp = drug_gene != "None"
  ) %>%
  arrange(CHR, POSITION)

chr_offsets <- ihs_all %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POSITION), .groups = "drop") %>%
  mutate(chr_offset = dplyr::lag(cumsum(chr_len), default = 0))

ihs_all <- ihs_all %>%
  left_join(chr_offsets, by = "CHR") %>%
  mutate(
    POS_cum = POSITION + chr_offset,
    logp    = LOGPVALUE
  )

axis_df <- ihs_all %>%
  group_by(CHR) %>%
  summarise(center = (min(POS_cum) + max(POS_cum)) / 2, .groups = "drop")

p_ihs <- ggplot(ihs_all,
                aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
  geom_point(alpha = 0.4, size = 1.0) +
  geom_point(data = ihs_all %>% filter(is_drug_snp),
             color = "red", size = 2.0) +
  facet_wrap(~category_name, ncol = 1) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center / 1e6) +
  scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
  geom_hline(yintercept = c(-ihs_thresh, ihs_thresh),
             linetype = "dashed", color = "black") +
  labs(x = "Chromosome", y = "iHS") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 8),
    strip.text      = element_text(size = 9)
  )

ggsave(file.path(plot_dir, "manhattan_ihs_per_group.png"),
       p_ihs, width = 10, height = 12, dpi = 300)

p_logp <- ggplot(ihs_all,
                 aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
  geom_point(alpha = 0.4, size = 1.0) +
  geom_point(data = ihs_all %>% filter(is_drug_snp),
             color = "red", size = 2.0) +
  facet_wrap(~category_name, ncol = 1) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center / 1e6) +
  scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
  geom_hline(yintercept = logp_thresh, linetype = "dashed", color = "black") +
  labs(x = "Chromosome", y = "-log10(p-value)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(size = 8),
    strip.text      = element_text(size = 9)
  )

ggsave(file.path(plot_dir, "manhattan_logp_ihs_per_group.png"),
       p_logp, width = 10, height = 12, dpi = 300)

# ─────────────────────────────────────────────────────────────
# Per-category standalone plots + TSVs
# ─────────────────────────────────────────────────────────────
all_groups <- sort(unique(ihs_all$category_name))

for (cat_nm in all_groups) {

  ihs_cat <- ihs_all %>%
    filter(category_name == cat_nm, CHR %in% 1:14)

  if (nrow(ihs_cat) == 0) next

  axis_df_cat <- ihs_cat %>%
    group_by(CHR) %>%
    summarise(center = (min(POS_cum) + max(POS_cum)) / 2, .groups = "drop")

  ihs_cat_annot <- ihs_cat %>%
    inner_join(genes_clean, by = "CHR") %>%
    filter(POSITION >= START, POSITION <= END)

  high_logp_tbl <- ihs_cat_annot %>%
    filter(logp > logp_thresh) %>%
    select(category_name, CHR, POSITION, IHS, logp,
           GENE_ID, GENE_NAME, PRODUCT,
           drug_gene, is_drug_snp)

  high_ihs_tbl <- ihs_cat_annot %>%
    filter(abs(IHS) > ihs_thresh) %>%
    select(category_name, CHR, POSITION, IHS, logp,
           GENE_ID, GENE_NAME, PRODUCT,
           drug_gene, is_drug_snp)

  write_tsv(high_logp_tbl, file.path(plot_dir, sprintf("iHS_%s_logp_high_snps.tsv", cat_nm)))
  write_tsv(high_ihs_tbl,  file.path(plot_dir, sprintf("iHS_%s_extreme_ihs_snps.tsv", cat_nm)))

  p_ihs_cat <- ggplot(ihs_cat,
                      aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_point(data = ihs_cat %>% filter(is_drug_snp),
               color = "red", size = 2.2) +
    scale_x_continuous(label = axis_df_cat$CHR, breaks = axis_df_cat$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    geom_hline(yintercept = c(-ihs_thresh, ihs_thresh),
               linetype = "dashed", color = "black") +
    labs(x = "Chromosome", y = "iHS", title = paste("iHS:", cat_nm)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir, sprintf("manhattan_ihs_%s.png", cat_nm)),
         p_ihs_cat, width = 10, height = 5, dpi = 300)

  p_logp_cat <- ggplot(ihs_cat,
                       aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_point(data = ihs_cat %>% filter(is_drug_snp),
               color = "red", size = 2.2) +
    scale_x_continuous(label = axis_df_cat$CHR, breaks = axis_df_cat$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    geom_hline(yintercept = logp_thresh,
               linetype = "dashed", color = "black") +
    labs(x = "Chromosome", y = "-log10(p-value)", title = paste("-log10(p):", cat_nm)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir, sprintf("manhattan_logp_ihs_%s.png", cat_nm)),
         p_logp_cat, width = 10, height = 5, dpi = 300)

  # DR-only labels: one peak SNP per drug gene
  label_logp_raw <- ihs_cat_annot %>% filter(is_drug_snp, logp > logp_thresh)
  label_ihs_raw  <- ihs_cat_annot %>% filter(is_drug_snp, abs(IHS) > ihs_thresh)

  label_logp <- label_logp_raw %>%
    group_by(drug_gene) %>%
    slice_max(logp, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(POSITION, CHR, logp, drug_gene, POS_cum)

  label_ihs <- label_ihs_raw %>%
    group_by(drug_gene) %>%
    slice_max(abs(IHS), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(POSITION, CHR, IHS, drug_gene, POS_cum)

  p_logp_lab <- ggplot(ihs_cat,
                       aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
    geom_point(alpha = 0.35, size = 1.2) +
    geom_hline(yintercept = logp_thresh,
               linetype = "dashed", color = "black") +
    geom_point(data = ihs_cat %>% filter(POSITION %in% label_logp$POSITION),
               color = "red", size = 2.4) +
    geom_text_repel(
      data = label_logp,
      aes(x = POS_cum / 1e6, y = logp, label = drug_gene),
      size = 2.8, max.overlaps = 10,
      segment.size = 0.2, box.padding = 0.3,
      min.segment.length = 0,
      color = "black"
    ) +
    scale_x_continuous(label = axis_df_cat$CHR, breaks = axis_df_cat$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    labs(x = "Chromosome", y = "-log10(p-value)",
         title = paste("-log10(p):", cat_nm, "(DR genes)")) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir, sprintf("manhattan_logp_ihs_%s_annotated.png", cat_nm)),
         p_logp_lab, width = 10, height = 5, dpi = 300)

  p_ihs_lab <- ggplot(ihs_cat,
                      aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
    geom_point(alpha = 0.35, size = 1.2) +
    geom_hline(yintercept = c(-ihs_thresh, ihs_thresh),
               linetype = "dashed", color = "black") +
    geom_point(data = ihs_cat %>% filter(POSITION %in% label_ihs$POSITION),
               color = "red", size = 2.4) +
    geom_text_repel(
      data = label_ihs,
      aes(x = POS_cum / 1e6, y = IHS, label = drug_gene),
      size = 2.8, max.overlaps = 10,
      segment.size = 0.2, box.padding = 0.3,
      min.segment.length = 0,
      color = "black"
    ) +
    scale_x_continuous(label = axis_df_cat$CHR, breaks = axis_df_cat$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    labs(x = "Chromosome", y = "iHS",
         title = paste("iHS:", cat_nm, "(DR genes)")) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir, sprintf("manhattan_ihs_%s_annotated.png", cat_nm)),
         p_ihs_lab, width = 10, height = 5, dpi = 300)
}

message("\nDone with iHS scanning + plotting.")
