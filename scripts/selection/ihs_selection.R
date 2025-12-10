#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(gt)
  library(tidyverse)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(rehh)
  library(tibble)
  library(ggrepel)
  library(data.table)
})

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-d", "--data-dir"), type = "character", default = ".",
              help = "Directory with scanned_haplotypes_*.tsv and where outputs will be written [default %default]",
              metavar = "character"),
  make_option(c("--country-list"), type = "character", default = NULL,
              help = "Text file with one country name per line (for scanned_haplotypes_<country>.tsv)",
              metavar = "character"),
  make_option(c("--genome-file"), type = "character", default = NULL,
              help = "Pf genome product annotation TSV (e.g. pf_genome_product_v3.tsv) [required]",
              metavar = "character"),
  make_option(c("--focus-pop"), type = "character", default = "Ethiopia",
              help = "Name of focal population for detailed plots (default %default)",
              metavar = "character"),
  make_option(c("--years"), type = "character", default = "2013,2017,2021",
              help = "Comma-separated list of years for focus-pop by-year iHS (scanned_haplotypes_<year>.tsv) [default %default]",
              metavar = "character"),
  make_option(c("--min-maf"), type = "numeric", default = 0.0,
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

opt <- parse_args(OptionParser(option_list = option_list))

data_dir    <- opt$`data-dir`
country_lst <- opt$`country-list`
genome_file <- opt$`genome-file`
focus_pop   <- opt$`focus-pop`
years_str   <- opt$years
min_maf     <- opt$`min-maf`
freqbin     <- opt$freqbin
ihs_thresh  <- opt$`ihs-thresh`
logp_thresh <- opt$`logp-thresh`

if (is.null(genome_file) || !file.exists(genome_file)) {
  stop("ERROR: --genome-file must be provided and must exist.\n", call. = FALSE)
}
if (is.null(country_lst) || !file.exists(country_lst)) {
  stop("ERROR: --country-list must be provided and must exist.\n", call. = FALSE)
}

years <- strsplit(years_str, ",")[[1]] |> trimws()
dir.create(file.path(data_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
plot_dir <- file.path(data_dir, "plots")

message("Data dir:     ", data_dir)
message("Genome file:  ", genome_file)
message("Country list: ", country_lst)
message("Focus pop:    ", focus_pop)
message("Years:        ", paste(years, collapse = ", "))

# ─────────────────────────────────────────────────────────────
# DR gene regions
# ─────────────────────────────────────────────────────────────

drug_genes <- tibble::tribble(
  ~gene,     ~CHR, ~START,   ~END,
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
# 1. Compute iHS per country
# ─────────────────────────────────────────────────────────────

countries <- readLines(country_lst) |> trimws()
countries <- countries[countries != ""]

if (length(countries) == 0) {
  stop("Country list is empty.\n", call. = FALSE)
}

all_ihs <- list()

for (country in countries) {
  file_path <- file.path(data_dir, sprintf("scanned_haplotypes_%s.tsv", country))
  if (!file.exists(file_path)) {
    warning("File not found for country ", country, ": ", file_path)
    next
  }

  cat("Computing iHS for country:", country, "\n")
  scan_data <- read.table(file_path, header = TRUE)
  ihs <- ihh2ihs(scan_data, min_maf = min_maf, freqbin = freqbin)

  if (!is.null(ihs$ihs)) {
    ihs_tbl <- as_tibble(ihs$ihs) %>%
      mutate(category_name = country)
    all_ihs[[country]] <- ihs_tbl
  } else {
    warning("No valid iHS markers for country: ", country)
  }
}

if (length(all_ihs) == 0) {
  stop("No iHS data could be computed for any country.\n", call. = FALSE)
}

ihs_all_countries <- bind_rows(all_ihs)
ihs_file <- file.path(data_dir, "iHS_all_countries.tsv")
write_tsv(ihs_all_countries, ihs_file)
message("Wrote combined iHS file: ", ihs_file)

# ─────────────────────────────────────────────────────────────
# 2. Global annotation of extreme |iHS|
# ─────────────────────────────────────────────────────────────

genes <- read_tsv(genome_file, show_col_types = FALSE)

genes_clean <- genes %>%
  mutate(CHR = as.numeric(str_extract(chr, "\\d{2}"))) %>%
  rename(
    START    = pos_start,
    END      = pos_end,
    GENE_ID  = gene_id,
    GENE_NAME= gene_name,
    PRODUCT  = gene_product
  ) %>%
  filter(!is.na(CHR)) %>%
  select(CHR, START, END, GENE_ID, GENE_NAME, PRODUCT)

ihs_all <- ihs_all_countries %>%
  filter(CHR %in% 1:14)

ihs_extreme <- ihs_all %>%
  filter(abs(IHS) > ihs_thresh)

ihs_annotated <- ihs_extreme %>%
  inner_join(genes_clean, by = "CHR") %>%
  filter(POSITION >= START, POSITION <= END) %>%
  select(category_name, CHR, POSITION, IHS, GENE_ID, GENE_NAME, PRODUCT)

annot_file <- file.path(data_dir, "iHS_extreme_sites_annotated.tsv")
write_tsv(ihs_annotated, annot_file)
message("Wrote extreme iHS annotation: ", annot_file)

# ─────────────────────────────────────────────────────────────
# 3. Per-country Manhattan plots
# ─────────────────────────────────────────────────────────────

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

# iHS
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

# -log10(p)
ihs_all <- ihs_all %>%
  mutate(logp = -log10(10^(-LOGPVALUE)))

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
# 4. Focus population (e.g. Ethiopia) – detailed plots + TSVs
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
    mutate(POS_cum = POSITION + chr_offset)

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
         title = paste("iHS", focus_pop)) +
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
  ihs_focus <- ihs_focus %>%
    mutate(logp = -log10(10^(-LOGPVALUE)))

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
    mutate(logp = -log10(10^(-LOGPVALUE))) %>%
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
         title = paste("-log10(p):", focus_pop)) +
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
         title = paste("iHS:", focus_pop)) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(size = 8),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir,
                   paste0("manhattan_ihs_", focus_pop, "_annotated.png")),
         p_ihs_lab, width = 12, height = 6)

  write_tsv(label_logp,
            file.path(plot_dir,
                      paste0("iHS_", focus_pop, "_logp_high_snps.tsv")))
  write_tsv(label_ihs,
            file.path(plot_dir,
                      paste0("iHS_", focus_pop, "_extreme_ihs_snps.tsv")))
} else {
  warning("No iHS data found for focus population: ", focus_pop)
}

# ─────────────────────────────────────────────────────────────
# 5. By-year iHS for focus population
# ─────────────────────────────────────────────────────────────

all_ihs_years <- list()

for (year in years) {
  scan_file <- file.path(data_dir, sprintf("scanned_haplotypes_%s.tsv", year))
  if (!file.exists(scan_file)) {
    warning("Scan file not found for year ", year, ": ", scan_file)
    next
  }

  cat("Computing iHS for year:", year, "\n")
  scan_data <- read.table(scan_file, header = TRUE)
  ihs_result <- ihh2ihs(scan_data, min_maf = min_maf, freqbin = freqbin)

  if (!is.null(ihs_result$ihs)) {
    ihs_df <- as_tibble(ihs_result$ihs) %>%
      mutate(year = year)
    all_ihs_years[[year]] <- ihs_df
  } else {
    warning("No valid iHS markers for year: ", year)
  }
}

if (length(all_ihs_years) > 0) {

  ihs_all_years <- bind_rows(all_ihs_years) %>%
    filter(CHR %in% 1:14) %>%
    mutate(logp = -log10(10^(-LOGPVALUE))) %>%
    arrange(year, CHR, POSITION)

  chr_offsets_y <- ihs_all_years %>%
    group_by(CHR) %>%
    summarise(chr_len = max(POSITION), .groups = "drop") %>%
    mutate(chr_offset = lag(cumsum(chr_len), default = 0))

  ihs_all_years <- ihs_all_years %>%
    left_join(chr_offsets_y, by = "CHR") %>%
    mutate(POS_cum = POSITION + chr_offset)

  ihs_all_years <- ihs_all_years %>%
    left_join(genes_clean, by = c("CHR")) %>%
    mutate(overlap_gene = POSITION >= START & POSITION <= END) %>%
    mutate(
      label_ihs  = ifelse(abs(IHS) > ihs_thresh & overlap_gene, GENE_ID, NA),
      label_logp = ifelse(logp > logp_thresh & overlap_gene, GENE_ID, NA)
    ) %>%
    left_join(
      drug_genes %>% rename(DRUG_GENE = gene,
                            DRUG_START = START,
                            DRUG_END = END),
      by = "CHR"
    ) %>%
    mutate(
      in_drug_region = !is.na(DRUG_GENE) &
        POSITION >= DRUG_START & POSITION <= DRUG_END
    )

  axis_df_y <- ihs_all_years %>%
    group_by(CHR) %>%
    summarise(center = (min(POS_cum) + max(POS_cum)) / 2, .groups = "drop")

  # iHS by year
  p_ihs_years <- ggplot(ihs_all_years,
                        aes(x = POS_cum / 1e6, y = IHS, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 0.9) +
    geom_point(data = ihs_all_years %>% filter(in_drug_region),
               color = "red", size = 1) +
    geom_hline(yintercept = c(-ihs_thresh, ihs_thresh),
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = ihs_all_years %>% filter(!is.na(label_ihs)),
      aes(label = label_ihs),
      size = 2.3, max.overlaps = 20,
      segment.size = 0.2, box.padding = 0.4
    ) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    facet_wrap(~year, ncol = 1) +
    scale_x_continuous(label = axis_df_y$CHR,
                       breaks = axis_df_y$center / 1e6) +
    labs(x = "Chromosome", y = "iHS", title = "iHS by Year") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text      = element_text(size = 10),
      axis.text.x     = element_text(size = 7),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir, "manhattan_iHS_by_year.png"),
         p_ihs_years, width = 12, height = 10)

  # -log10(p) by year
  p_logp_years <- ggplot(ihs_all_years,
                         aes(x = POS_cum / 1e6, y = logp, color = factor(CHR))) +
    geom_point(alpha = 0.4, size = 0.9) +
    geom_point(data = ihs_all_years %>% filter(in_drug_region),
               color = "red", size = 1) +
    geom_hline(yintercept = logp_thresh,
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = ihs_all_years %>% filter(!is.na(label_logp)),
      aes(label = label_logp),
      size = 2.3, max.overlaps = 25,
      segment.size = 0.2, box.padding = 0.4
    ) +
    scale_color_manual(values = rep(c("grey30", "grey60"), 7)) +
    facet_wrap(~year, ncol = 1) +
    scale_x_continuous(label = axis_df_y$CHR,
                       breaks = axis_df_y$center / 1e6) +
    labs(x = "Chromosome", y = "-log10(p-value)", title = "-log10(p) by Year") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.text      = element_text(size = 10),
      axis.text.x     = element_text(size = 7),
      plot.title      = element_text(size = 14, face = "bold")
    )

  ggsave(file.path(plot_dir, "manhattan_logp_by_year.png"),
         p_logp_years, width = 12, height = 10)

  # Matched tables
  label_df_logp_y <- ihs_all_years %>%
    filter(!is.na(label_logp)) %>%
    arrange(desc(logp)) %>%
    select(year, CHR, POSITION, logp, GENE_ID, GENE_NAME, PRODUCT)

  label_df_ihs_y <- ihs_all_years %>%
    filter(!is.na(label_ihs)) %>%
    arrange(desc(abs(IHS))) %>%
    select(year, CHR, POSITION, IHS, GENE_ID, GENE_NAME, PRODUCT)

  write_tsv(label_df_logp_y,
            file.path(plot_dir, "iHS_high_logp_by_year.tsv"))
  write_tsv(label_df_ihs_y,
            file.path(plot_dir, "iHS_high_abs_iHS_by_year.tsv"))
} else {
  warning("No by-year iHS results computed (check scanned_haplotypes_<year>.tsv files).")
}

message("Done with iHS plotting/selection.")

