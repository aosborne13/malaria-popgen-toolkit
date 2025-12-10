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
  m[m < 0] <- NA          # just in case
  p <- rowMeans(m == 1, na.rm = TRUE)  # alt allele freq
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
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "SNP-level annotation file (chr,pos,ref,alt,...) [required]",
              metavar = "character"),
  make_option(c("-c", "--category"), type = "character", default = NULL,
              help = "Category name(s), comma-separated, or 'ALL' for all categories",
              metavar = "character"),
  make_option(c("--label_category"), type = "character", default = "country",
              help = "Column name in metadata for category [default %default]",
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

  # Gene-level annotation for iHS
  make_option(c("--genome-file"), type = "character", default = NULL,
              help = "Pf genome product annotation TSV (pf_genome_product_v3.tsv) [required]",
              metavar = "character"),
  make_option(c("--focus-pop"), type = "character", default = NULL,
              help = "Focus population for detailed plots [default: first category]",
              metavar = "character"),
  make_option(c("--min-maf-ihs"), type = "numeric", default = 0.0,
              help = "min_maf argument to ihh2ihs [default %default]",
              metavar = "numeric"),
  make_option(c("--freqbin"), type = "numeric", default = 0.05,
              help = "freqbin argument to ihh2ihs [default %default]",
              metavar = "numeric"),
  make_option(c("--ihs-thresh"), type = "numeric", default = 2.0,
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
annotation_file <- opt$annotation
category_arg    <- opt$category
label_category  <- opt$label_category
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

if (is.null(bin_mat_file)) {
  stop("ERROR: --matrix_binary is required\n", call. = FALSE)
}
if (!file.exists(bin_mat_file)) {
  stop("ERROR: matrix_binary file not found: ", bin_mat_file, "\n", call. = FALSE)
}
if (is.null(metadata_file)) {
  stop("ERROR: --metadata is required\n", call. = FALSE)
}
if (!file.exists(metadata_file)) {
  stop("ERROR: metadata file not found: ", metadata_file, "\n", call. = FALSE)
}
if (is.null(annotation_file)) {
  stop("ERROR: --annotation is required\n", call. = FALSE)
}
if (!file.exists(annotation_file)) {
  stop("ERROR: annotation file not found: ", annotation_file, "\n", call. = FALSE)
}
if (is.null(genome_file) || !file.exists(genome_file)) {
  stop("ERROR: --genome-file must be provided and must exist.\n", call. = FALSE)
}

dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
setDTthreads(threads)

message("Workdir:        ", workdir)
message("Binary matrix:  ", bin_mat_file)
message("Metadata:       ", metadata_file)
message("Annotation:     ", annotation_file)
message("Genome file:    ", genome_file)

# ─────────────────────────────────────────────────────────────
# Load annotation & metadata
# ─────────────────────────────────────────────────────────────
annotation_raw <- read_tsv(annotation_file, show_col_types = FALSE)
metadata_all   <- read_tsv(metadata_file, show_col_types = FALSE)

needed_cols <- c(label_id, label_category, label_fws)
if (!all(needed_cols %in% colnames(metadata_all))) {
  stop("ERROR: Missing columns in metadata: ",
       paste(setdiff(needed_cols, colnames(metadata_all)), collapse = ", "),
       "\n", call. = FALSE)
}

metadata_all[[label_category]] <- as.character(metadata_all[[label_category]])
metadata_all[[label_id]]       <- as.character(metadata_all[[label_id]])

available_categories <- metadata_all %>%
  dplyr::select(all_of(label_category)) %>%
  distinct() %>%
  pull() %>%
  sort()

if (length(available_categories) == 0) {
  stop("ERROR: No categories found in metadata column '", label_category, "'.\n",
       call. = FALSE)
}

message("Available categories: ", paste(available_categories, collapse = ", "))

# Decide which categories to run
if (is.null(category_arg) || category_arg == "" ||
    tolower(category_arg) == "all") {
  categories <- available_categories
  message("No specific category provided; using ALL categories.")
} else {
  categories <- strsplit(category_arg, ",")[[1]] |> trimws()
  missing_cats <- setdiff(categories, available_categories)
  if (length(missing_cats) > 0) {
    stop("ERROR: Requested category(ies) not found in metadata: ",
         paste(missing_cats, collapse = ", "), "\n",
         "Available: ", paste(available_categories, collapse = ", "), "\n",
         call. = FALSE)
  }
  message("Using selected category(ies): ", paste(categories, collapse = ", "))
}

# Decide focus population (default: first category after filtering)
if (is.null(focus_pop) || focus_pop == "" || !(focus_pop %in% categories)) {
  if (!is.null(focus_pop) && focus_pop != "" && !(focus_pop %in% categories)) {
    warning("Requested focus-pop '", focus_pop,
            "' not in selected categories; defaulting to first category.\n")
  }
  focus_pop <- categories[1]
}
message("Focus population: ", focus_pop)

# ─────────────────────────────────────────────────────────────
# Define drug-resistance gene regions
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
# Load genome product annotation
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
# Main loop: for each category
# ─────────────────────────────────────────────────────────────
all_ihs <- list()

for (category in categories) {
  message("\n==============================")
  message("Category: ", category)
  message("==============================")

  meta_cat <- metadata_all %>%
    filter(.data[[label_category]] == category)

  if (nrow(meta_cat) == 0) {
    warning("  Category ", category, " has no metadata rows; skipping.")
    next
  }

  category_str <- gsub(" ", "_", category)

  # Columns to select from binary matrix: chr,pos,ref + samples in this category
  sample_ids <- meta_cat %>%
    dplyr::select(all_of(label_id)) %>%
    pull() %>%
    unique()

  if (length(sample_ids) == 0) {
    warning("  Category ", category, ": no sample IDs in metadata; skipping.")
    next
  }

  samples_select <- c("chr", "pos", "ref", sample_ids)

  message("  Loading SNP matrix subset for category (",
          length(sample_ids), " samples)...")

  snp <- fread(bin_mat_file, sep = "\t",
               select = samples_select,
               header = TRUE, data.table = FALSE)

  if (!all(c("chr", "pos", "ref") %in% colnames(snp)[1:3])) {
    stop("ERROR: Matrix must start with columns: chr, pos, ref\n", call. = FALSE)
  }

  # Copy annotation so we don't mutate the global version
  annotation <- annotation_raw

  # Filter chromosomes if requested
  if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
    rm_chr <- strsplit(rm_chr_str, ",")[[1]] |> trimws()
    if (any(rm_chr %in% unique(snp$chr))) {
      message("  Removing chromosomes: ", paste(rm_chr, collapse = ", "))
      snp        <- snp        %>% filter(!chr %in% rm_chr)
      annotation <- annotation %>% filter(!chr %in% rm_chr)
    } else {
      warning("  None of the chromosomes in --remove_chr found in matrix for this category; continuing.")
    }
  } else {
    message("  No chromosomes removed (API/mito may still be present).")
  }

  # Transform chromosome from string to numeric via regex
  snp$chr        <- as.numeric(stringr::str_match(snp$chr, pattern)[, groupid])
  annotation$chr <- as.numeric(stringr::str_match(annotation$chr, pattern)[, groupid])

  # Keep metadata rows whose IDs are present in matrix columns
  meta_cat <- meta_cat %>%
    filter(.data[[label_id]] %in% colnames(snp)[-(1:3)])
  if (nrow(meta_cat) < 2) {
    warning("  Category ", category, ": fewer than 2 samples in matrix; skipping.")
    next
  }

  # Re-check ordering / overlap
  if (!all(meta_cat[[label_id]] %in% colnames(snp)[-(1:3)])) {
    warning("  Some metadata samples are missing from matrix for category ", category)
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
    # NA -> 0 (REF), 0.5 -> 1 (ALT)
    maj3[is.na(maj3)] <- 0
    maj3[maj3 == 0.5] <- 1
  } else if (forced_mixed) {
    # 0.5 -> 1, NA remain NA
    maj3[maj3 == 0.5] <- 1
  } else {
    # default: mixed -> NA
    maj3[maj3 == 0.5] <- NA
  }
  rm(snp_c)

  # MAF filter
  maf_vals <- calculate_maf(maj3)
  to_keep_sti <- which(!is.na(maf_vals) & maf_vals >= th_maf)
  rm(maf_vals)

  if (length(to_keep_sti) == 0) {
    warning("  Category ", category, ": no SNPs pass MAF ≥ ", th_maf, "; skipping.")
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
    warning("  Category ", category, ": fewer than 2 samples pass Fws ≥ ",
            threshold_fws, "; skipping.")
    next
  }

  maj4 <- maj4[, samples_keep, drop = FALSE]

  # Build SNP map and join with annotation
  snp_annot <- snp_d_filt %>%
    left_join(annotation, by = c("chr", "pos", "ref")) %>%
    tidyr::unite("info", c(chr, pos), sep = "_", remove = FALSE)

  if (!all(c("alt") %in% colnames(snp_annot))) {
    stop("ERROR: annotation must contain at least columns chr,pos,ref,alt.\n",
         call. = FALSE)
  }

  map <- snp_annot %>%
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
  i <- sapply(hap_c, is.character)
  hap_c[i] <- lapply(hap_c[i], function(x) suppressWarnings(as.numeric(x)))

  # Create haplotypes per chromosome and run scan_hh
  u_chr <- sort(unique(map$chr))
  results_hh <- NULL

  log_file <- file.path(workdir, paste0(category_str, "_rehh.log"))
  sink(log_file, append = FALSE)
  cat("## Create haplotypes: data2haplohh() & scan_hh() for category ",
      category_str, " ##\n", sep = "")

  for (uchr in u_chr) {
    cat("Chr ", uchr, " ...\n", sep = "")

    # indices of markers for this chromosome in the map
    idx_chr <- which(map$chr == uchr)

    if (length(idx_chr) == 0) {
      cat("  No markers for chr ", uchr, " after filtering; skipping.\n", sep = "")
      next
    }

    # subset haplotypes and map for this chromosome
    hap_chr_s <- hap_c[, idx_chr, drop = FALSE]
    map_chr   <- map[idx_chr, , drop = FALSE]

    # sanity check: genotype columns = markers in map
    if (ncol(hap_chr_s) != nrow(map_chr)) {
      stop("For category ", category_str, " chr ", uchr,
           ": number of markers in hap (", ncol(hap_chr_s),
           ") != markers in map (", nrow(map_chr), ").", call. = FALSE)
    }

    # add explicit unique hap IDs as first column
    hap_ids <- paste0("h", seq_len(nrow(hap_chr_s)))
    hap_chr_out <- cbind(hap_ids, hap_chr_s)

    hap_file <- file.path(workdir, sprintf("hap_chr%d_%s", uchr, category_str))
    write.table(hap_chr_out, hap_file,
                sep = "\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

    map_file_chr <- file.path(workdir,
                              sprintf("snp.info.chr%d.%s", uchr, category_str))
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
    warning("  scan_hh produced no results for category ", category, "; skipping.")
    next
  }

  # Store scan_hh results per category
  scanned_file <- file.path(workdir,
                            sprintf("scanned_haplotypes_%s.tsv", category_str))
  results_hh$category <- category_str
  write.table(results_hh, scanned_file,
              quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  message("  Wrote scan_hh results: ", scanned_file)

  # Compute iHS for this category
  message("  Computing iHS via ihh2ihs()...")
  ihs_obj <- ihh2ihs(results_hh, min_maf = min_maf_ihs, freqbin = freqbin)

  if (!is.null(ihs_obj$ihs)) {
    ihs_tbl <- as_tibble(ihs_obj$ihs) %>%
      mutate(category_name = category)
    all_ihs[[category]] <- ihs_tbl
    message("  iHS markers: ", nrow(ihs_tbl))
  } else {
    warning("  No valid iHS markers for category: ", category)
  }
}

if (length(all_ihs) == 0) {
  stop("No iHS data could be computed for any category.\n", call. = FALSE)
}

# ─────────────────────────────────────────────────────────────
# Combine iHS for all categories, annotate extreme sites
# ─────────────────────────────────────────────────────────────
ihs_all_countries <- bind_rows(all_ihs)
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
# Per-category Manhattan plots
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
  mutate(chr_offset = lag(cumsum(chr_len), default = 0))

ihs_all <- ihs_all %>%
  left_join(chr_offsets, by = "CHR") %>%
  mutate(POS_cum = POSITION + chr_offset)

axis_df <- ihs_all %>%
  group_by(CHR) %>%
  summarise(center = (min(POS_cum) + max(POS_cum)) / 2, .groups = "drop")

# iHS plot
p_ihs <- ggplot(ihs_all, aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_point(data = ihs_all %>% filter(is_drug_snp),
             color = "red", size = 1) +
  facet_wrap(~category_name, ncol = 4) +
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

ggsave(file.path(plot_dir, "manhattan_ihs_per_country.png"),
       p_ihs, width = 16, height = 10)

# -log10(p) – LOGPVALUE from ihh2ihs is already -log10(p)
ihs_all <- ihs_all %>%
  mutate(logp = LOGPVALUE)

p_logp <- ggplot(ihs_all, aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_point(data = ihs_all %>% filter(is_drug_snp),
             color = "red", size = 1) +
  facet_wrap(~category_name, ncol = 4) +
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

ggsave(file.path(plot_dir, "manhattan_logp_ihs_per_country.png"),
       p_logp, width = 16, height = 10)

# ─────────────────────────────────────────────────────────────
# Focus population – detailed plots + TSVs
# ─────────────────────────────────────────────────────────────
ihs_focus <- ihs_all_countries %>%
  filter(category_name == focus_pop, CHR %in% 1:14)

if (nrow(ihs_focus) > 0) {

  ihs_focus <- ihs_focus %>%
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

  chr_offsets_f <- ihs_focus %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POSITION), .groups = "drop") %>%
    mutate(chr_offset = lag(cumsum(chr_len), default = 0))

  ihs_focus <- ihs_focus %>%
    left_join(chr_offsets_f, by = "CHR") %>%
    mutate(POS_cum = POSITION + chr_offset,
           logp    = LOGPVALUE)

  axis_df_f <- ihs_focus %>%
    group_by(CHR) %>%
    summarize(center = (min(POS_cum) + max(POS_cum)) / 2, .groups = "drop")

  # iHS
  p_ihs_f <- ggplot(ihs_focus,
                    aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_point(data = ihs_focus %>% filter(is_drug_snp),
               color = "red", size = 1.5) +
    scale_x_continuous(label = axis_df_f$CHR,
                       breaks = axis_df_f$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    geom_hline(yintercept = c(-ihs_thresh, ihs_thresh),
               linetype = "dashed", color = "black") +
    labs(x = "Chromosome", y = "iHS",
         title = paste("iHS:", focus_pop)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir,
                   paste0("manhattan_ihs_", focus_pop, ".png")),
         p_ihs_f, width = 12, height = 6)

  # -log10(p)
  p_logp_f <- ggplot(ihs_focus,
                     aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_point(data = ihs_focus %>% filter(is_drug_snp),
               color = "red", size = 1.5) +
    scale_x_continuous(label = axis_df_f$CHR,
                       breaks = axis_df_f$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    geom_hline(yintercept = logp_thresh,
               linetype = "dashed", color = "black") +
    labs(x = "Chromosome", y = "-log10(p-value)",
         title = paste("-log10(p):", focus_pop)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir,
                   paste0("manhattan_logp_ihs_", focus_pop, ".png")),
         p_logp_f, width = 12, height = 6)

  # Annotated / labeled versions
  ihs_focus_ord <- ihs_focus %>%
    arrange(CHR, POSITION)

  ihs_focus_full <- ihs_focus_ord %>%
    inner_join(genes_clean, by = "CHR") %>%
    filter(POSITION >= START, POSITION <= END)

  label_logp <- ihs_focus_full %>%
    filter(logp > logp_thresh) %>%
    distinct(POSITION, CHR, logp, GENE_ID, GENE_NAME, PRODUCT)

  label_ihs <- ihs_focus_full %>%
    filter(abs(IHS) > ihs_thresh) %>%
    distinct(POSITION, CHR, IHS, GENE_ID, GENE_NAME, PRODUCT)

  ihs_focus_pos <- ihs_focus_ord %>%
    select(POSITION, CHR, POS_cum)

  label_logp <- label_logp %>%
    left_join(ihs_focus_pos, by = c("POSITION", "CHR"))

  label_ihs <- label_ihs %>%
    left_join(ihs_focus_pos, by = c("POSITION", "CHR"))

  p_logp_lab <- ggplot(ihs_focus_ord,
                       aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_hline(yintercept = logp_thresh,
               linetype = "dashed", color = "black") +
    geom_point(data = ihs_focus_ord %>%
                 filter(POSITION %in% label_logp$POSITION),
               color = "red", size = 1.5) +
    geom_text_repel(
      data = label_logp,
      aes(x = POS_cum / 1e6, y = logp, label = GENE_ID),
      size = 2.7, max.overlaps = 30,
      segment.size = 0.2, box.padding = 0.4
    ) +
    scale_x_continuous(label = axis_df_f$CHR,
                       breaks = axis_df_f$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    labs(x = "Chromosome", y = "-log10(p-value)",
         title = paste("-log10(p):", focus_pop, "(annotated)")) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir,
                   paste0("manhattan_logp_ihs_", focus_pop, "_annotated.png")),
         p_logp_lab, width = 12, height = 6)

  p_ihs_lab <- ggplot(ihs_focus_ord,
                      aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_hline(yintercept = c(-ihs_thresh, ihs_thresh),
               linetype = "dashed", color = "black") +
    geom_point(data = ihs_focus_ord %>%
                 filter(POSITION %in% label_ihs$POSITION),
               color = "red", size = 1.5) +
    geom_text_repel(
      data = label_ihs,
      aes(x = POS_cum / 1e6, y = IHS, label = GENE_ID),
      size = 2.7, max.overlaps = 25,
      segment.size = 0.2, box.padding = 0.5
    ) +
    scale_x_continuous(label = axis_df_f$CHR,
                       breaks = axis_df_f$center / 1e6) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    labs(x = "Chromosome", y = "iHS",
         title = paste("iHS:", focus_pop, "(annotated)")) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir,
                   paste0("manhattan_ihs_", focus_pop, "_annotated.png")),
         p_ihs_lab, width = 12, height = 6)

  # Save annotated tables
  write_tsv(label_logp,
            file.path(plot_dir,
                      paste0("iHS_", focus_pop, "_logp_high_snps.tsv")))
  write_tsv(label_ihs,
            file.path(plot_dir,
                      paste0("iHS_", focus_pop, "_extreme_ihs_snps.tsv")))
} else {
  warning("No iHS data found for focus population: ", focus_pop)
}

message("\nDone with iHS scanning + plotting.")

