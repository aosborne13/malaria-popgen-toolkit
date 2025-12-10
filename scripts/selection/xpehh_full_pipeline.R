#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(rehh)
})

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Directory containing scanned_haplotypes_<pop>.tsv and where outputs will be written [default %default]",
              metavar = "character"),
  make_option(c("--scan_prefix"), type = "character", default = "scanned_haplotypes_",
              help = "Prefix of scan_hh result files (e.g. scanned_haplotypes_) [default %default]",
              metavar = "character"),
  make_option(c("--ref_pop"), type = "character", default = NULL,
              help = "Reference population label (e.g. Ethiopia) [required]",
              metavar = "character"),
  make_option(c("--targets"), type = "character", default = NULL,
              help = "Comma-separated list of target populations (e.g. Cameroon,DRC,Ghana). If omitted, all populations except ref_pop are used.",
              metavar = "character"),
  make_option(c("--genome-file"), type = "character", default = NULL,
              help = "Pf genome product annotation TSV (pf_genome_product_v3.tsv) [required]",
              metavar = "character"),
  make_option(c("--xpehh-thresh"), type = "numeric", default = 2.0,
              help = "Absolute XP-EHH threshold for 'extreme' sites [default %default]",
              metavar = "numeric"),
  make_option(c("--logp-thresh"), type = "numeric", default = 1.3,
              help = "-log10(p) threshold for 'extreme' sites [default %default]",
              metavar = "numeric"),
  make_option(c("--out_prefix"), type = "character", default = NULL,
              help = "Prefix for combined output files (default: <ref_pop>_vs_all)",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

workdir      <- opt$workdir
scan_prefix  <- opt$scan_prefix
ref_pop      <- opt$ref_pop
targets_arg  <- opt$targets
genome_file  <- opt$`genome-file`
xpehh_thresh <- opt$`xpehh-thresh`
logp_thresh  <- opt$`logp-thresh`
out_prefix   <- opt$out_prefix

if (is.null(ref_pop) || ref_pop == "") {
  stop("ERROR: --ref_pop is required.\n", call. = FALSE)
}
if (is.null(genome_file) || !file.exists(genome_file)) {
  stop("ERROR: --genome-file must be provided and must exist.\n", call. = FALSE)
}
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)

if (is.null(out_prefix) || out_prefix == "") {
  out_prefix <- paste0(ref_pop, "_vs_all")
}

message("Workdir:      ", workdir)
message("Scan prefix:  ", scan_prefix)
message("Ref pop:      ", ref_pop)
message("Genome file:  ", genome_file)
message("Out prefix:   ", out_prefix)

# ─────────────────────────────────────────────────────────────
# Determine target populations (auto-detect like iHS)
# ─────────────────────────────────────────────────────────────

scan_files <- list.files(workdir, pattern = paste0("^", scan_prefix), full.names = FALSE)
if (length(scan_files) == 0) {
  stop("ERROR: No files matching '", scan_prefix, "*' found in workdir.\n", call. = FALSE)
}

all_pops <- scan_files |>
  gsub(paste0("^", scan_prefix), "", x = _) |>
  gsub("\\.tsv$", "", x = _) |>
  unique()

if (!(ref_pop %in% all_pops)) {
  stop("ERROR: ref_pop '", ref_pop, "' has no file ", scan_prefix, ref_pop, ".tsv in workdir.\n",
       "Available pops: ", paste(all_pops, collapse = ", "), "\n",
       call. = FALSE)
}

if (!is.null(targets_arg) && targets_arg != "") {
  requested <- strsplit(targets_arg, ",")[[1]] |> trimws()
  missing   <- setdiff(requested, all_pops)
  if (length(missing) > 0) {
    warning("Requested targets not found in scanned files: ",
            paste(missing, collapse = ", "), "\n",
            "Available pops: ", paste(all_pops, collapse = ", "))
  }
  targets <- intersect(requested, all_pops)
} else {
  # default: all pops except ref
  targets <- setdiff(all_pops, ref_pop)
}

if (length(targets) == 0) {
  stop("ERROR: No target populations after filtering (ref_pop removed).\n",
       "Available pops: ", paste(all_pops, collapse = ", "), "\n",
       call. = FALSE)
}

message("Target populations: ", paste(targets, collapse = ", "))

# ─────────────────────────────────────────────────────────────
# Load reference scan_hh
# ─────────────────────────────────────────────────────────────
ref_file <- file.path(workdir, sprintf("%s%s.tsv", scan_prefix, ref_pop))
message("Loading reference scan_hh: ", ref_file)
scan_ref <- read.table(ref_file, header = TRUE)

# ─────────────────────────────────────────────────────────────
# Load genome annotation, define DR genes
# ─────────────────────────────────────────────────────────────
genes <- readr::read_tsv(genome_file, show_col_types = FALSE)

genes_clean <- genes %>%
  mutate(CHR = as.numeric(stringr::str_extract(chr, "\\d{2}"))) %>%
  dplyr::rename(
    START     = pos_start,
    END       = pos_end,
    GENE_ID   = gene_id,
    GENE_NAME = gene_name,
    PRODUCT   = gene_product
  ) %>%
  dplyr::filter(!is.na(CHR)) %>%
  dplyr::select(CHR, START, END, GENE_ID, GENE_NAME, PRODUCT)

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
# XP-EHH per comparison
# ─────────────────────────────────────────────────────────────
xpehh_list <- list()

for (target in targets) {
  targ_file <- file.path(workdir, sprintf("%s%s.tsv", scan_prefix, target))
  if (!file.exists(targ_file)) {
    warning("Target scan_hh file not found for ", target, ": ", targ_file)
    next
  }

  message("\n▶ ", ref_pop, " vs ", target)
  scan_target <- read.table(targ_file, header = TRUE)

  xpehh_out <- ies2xpehh(scan_pop1 = scan_ref, scan_pop2 = scan_target)

  df <- as_tibble(xpehh_out)
  # drop rows that are all NA (can happen for some markers)
  df <- df[!apply(df, 1, function(x) all(is.na(x))), , drop = FALSE]

  if (nrow(df) == 0) {
    warning("  No XP-EHH results for comparison ", ref_pop, " vs ", target)
    next
  }

  df <- df %>%
    mutate(
      ref_pop    = ref_pop,
      target_pop = target,
      comparison = paste(ref_pop, "vs", target, sep = "_"),
      CHR        = as.numeric(CHR),
      POSITION   = as.numeric(POSITION)
    )

  # save per-comparison table
  out_comp <- file.path(workdir,
                        sprintf("xpehh_%s_vs_%s.tsv", ref_pop, target))
  readr::write_tsv(df, out_comp)
  message("  Wrote: ", out_comp)

  xpehh_list[[target]] <- df
}

if (length(xpehh_list) == 0) {
  stop("No XP-EHH results produced for any comparison.\n", call. = FALSE)
}

xpehh_all <- bind_rows(xpehh_list)

# ─────────────────────────────────────────────────────────────
# Combined output & annotation
# ─────────────────────────────────────────────────────────────
combined_file <- file.path(workdir,
                           sprintf("%s_xpehh_all.tsv", out_prefix))
readr::write_tsv(xpehh_all, combined_file)
message("\nWrote combined XP-EHH file: ", combined_file)

# Only nuclear chromosomes, with gene annotation
xpehh_all_nuc <- xpehh_all %>%
  filter(CHR %in% 1:14)

xpehh_annotated <- xpehh_all_nuc %>%
  inner_join(genes_clean, by = "CHR") %>%
  filter(POSITION >= START, POSITION <= END) %>%
  select(comparison, ref_pop, target_pop,
         CHR, POSITION, XPEHH, LOGPVALUE,
         GENE_ID, GENE_NAME, PRODUCT)

# Mark DR regions on annotated table as well
xpehh_annotated <- xpehh_annotated %>%
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
    is_drug_snp = drug_gene != "None",
    logp        = LOGPVALUE  # LOGPVALUE already -log10(p)
  )

xpehh_extreme <- xpehh_annotated %>%
  filter((XPEHH > xpehh_thresh | XPEHH < -xpehh_thresh) &
           logp > logp_thresh)

extreme_file <- file.path(workdir,
                          sprintf("%s_xpehh_extreme_sites_annotated.tsv",
                                  out_prefix))
readr::write_tsv(xpehh_extreme, extreme_file)
message("Wrote extreme XP-EHH annotation: ", extreme_file)

# ─────────────────────────────────────────────────────────────
# Manhattan plots per comparison
# ─────────────────────────────────────────────────────────────
plot_dir <- file.path(workdir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (cmp in unique(xpehh_all$comparison)) {

  df <- xpehh_all %>% filter(comparison == cmp, CHR %in% 1:14)

  if (nrow(df) == 0) next

  df <- df %>%
    filter(!is.na(CHR), !is.na(POSITION)) %>%
    arrange(CHR, POSITION)

  # Build cumulative positions
  chr_lengths <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POSITION), .groups = "drop") %>%
    mutate(
      chr_start  = cumsum(dplyr::lag(chr_len, default = 0)),
      chr_center = chr_start + chr_len / 2
    )

  df <- df %>%
    left_join(chr_lengths, by = "CHR") %>%
    mutate(
      pos_cum = POSITION + chr_start,
      logp    = LOGPVALUE
    )

  # flag DR SNPs
  df <- df %>%
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
    )

  # pretty label for title
  title_str <- gsub("_", " ", cmp)

  # XP-EHH Manhattan
  p_xpehh <- ggplot(df, aes(x = pos_cum / 1e6, y = XPEHH, color = as.factor(CHR))) +
    geom_point(size = 0.9, alpha = 0.7) +
    geom_point(data = df %>% filter(is_drug_snp),
               color = "red", size = 1.6) +
    geom_hline(yintercept = c(-xpehh_thresh, xpehh_thresh),
               linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(
      breaks = chr_lengths$chr_center / 1e6,
      labels = chr_lengths$CHR
    ) +
    scale_color_manual(values = rep(c("grey30", "grey60"), length(unique(df$CHR)))) +
    labs(
      title = paste("XP-EHH:", title_str),
      x = "Chromosome",
      y = "XP-EHH"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(filename = file.path(plot_dir,
                              sprintf("XP-EHH_%s.png", cmp)),
         plot = p_xpehh, width = 10, height = 4, dpi = 300)

  # -log10(p) Manhattan
  p_logp <- ggplot(df, aes(x = pos_cum / 1e6, y = logp, color = as.factor(CHR))) +
    geom_point(size = 0.9, alpha = 0.7) +
    geom_point(data = df %>% filter(is_drug_snp),
               color = "red", size = 1.6) +
    geom_hline(yintercept = logp_thresh,
               linetype = "dashed", color = "black", linewidth = 0.4) +
    scale_x_continuous(
      breaks = chr_lengths$chr_center / 1e6,
      labels = chr_lengths$CHR
    ) +
    scale_color_manual(values = rep(c("grey30", "grey60"), length(unique(df$CHR)))) +
    labs(
      title = paste("-log10(p):", title_str),
      x = "Chromosome",
      y = "-log10(p)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(filename = file.path(plot_dir,
                              sprintf("XP-EHH_logp_%s.png", cmp)),
         plot = p_logp, width = 10, height = 4, dpi = 300)
}

message("\nDone with XP-EHH comparisons and plotting.")
