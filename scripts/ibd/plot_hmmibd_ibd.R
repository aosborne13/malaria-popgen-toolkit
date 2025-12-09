#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(data.table)
  library(ggrepel)
})

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Directory containing <suffix>_hmmIBD_ibd_win<kb>.tsv and <suffix>_hmmIBD_fraction.tsv [default %default]",
              metavar = "character"),
  make_option(c("--ref_index"), type = "character", default = NULL,
              help = "FAI index for reference genome (required)",
              metavar = "character"),
  make_option(c("--gene_product"), type = "character", default = NULL,
              help = "Gene product annotation TSV (pf_genome_product_v3.tsv) [required for gene/drug highlighting]",
              metavar = "character"),
  make_option(c("--suffix"), type = "character", default = NULL,
              help = "Prefix used by summarise_hmmibd_windows.R (e.g. 13_08_2025) [required]",
              metavar = "character"),
  make_option(c("-w", "--window_size"), type = "numeric", default = 50000,
              help = "Sliding window size (bp) used in summarise_hmmibd_windows.R [default %default]",
              metavar = "numeric"),
  make_option(c("--quantile_cutoff"), type = "numeric", default = 0.95,
              help = "Quantile cutoff used for high-IBD windows [default %default]",
              metavar = "numeric"),
  make_option(c("--remove_chr"), type = "character",
              default = "Pf3D7_API_v3,Pf3D7_MIT_v3",
              help = "Comma-separated chromosomes to drop (default %default)",
              metavar = "character"),
  make_option(c("--regex_chr"), type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex to derive numeric chromosome from chr name [default %default]",
              metavar = "character"),
  make_option(c("--regex_groupid"), type = "integer", default = 3,
              help = "Capture group index in regex_chr giving numeric chromosome [default %default]",
              metavar = "integer"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for plots (default: <workdir>/win_<kb>)",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

workdir        <- opt$workdir
ref_index      <- opt$ref_index
gene_file      <- opt$gene_product
suffix         <- opt$suffix
win_bp         <- opt$window_size
win_kb         <- win_bp / 1000
th_quantile    <- opt$quantile_cutoff
rm_chr_str     <- opt$remove_chr
pattern        <- opt$regex_chr
groupid        <- opt$regex_groupid
outdir         <- opt$outdir

if (is.null(ref_index) || !file.exists(ref_index)) {
  stop("ERROR: --ref_index must be provided and must exist.\n", call. = FALSE)
}
if (is.null(suffix)) {
  stop("ERROR: --suffix is required (e.g. 13_08_2025).\n", call. = FALSE)
}
if (!dir.exists(workdir)) {
  stop("ERROR: workdir does not exist: ", workdir, "\n", call. = FALSE)
}

if (is.null(outdir)) {
  outdir <- file.path(workdir, sprintf("win_%dkb", win_kb))
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Workdir: ", workdir)
message("Outdir:  ", outdir)
message("Suffix:  ", suffix)
message("Window:  ", win_bp, " bp")

# ─────────────────────────────────────────────────────────────
# Load combined IBD & fraction files
# ─────────────────────────────────────────────────────────────

ibd_file <- file.path(workdir,
                      sprintf("%s_hmmIBD_ibd_win%dkb.tsv", suffix, win_kb))
frac_file <- file.path(workdir,
                       sprintf("%s_hmmIBD_fraction.tsv", suffix))

if (!file.exists(ibd_file)) {
  stop("ERROR: IBD window file not found: ", ibd_file, "\n", call. = FALSE)
}
if (!file.exists(frac_file)) {
  stop("ERROR: fraction file not found: ", frac_file, "\n", call. = FALSE)
}

combined_ibd <- read_tsv(ibd_file, col_types = cols())
fraction_ibd <- read_tsv(frac_file, col_types = cols())

if (!"category" %in% colnames(combined_ibd)) {
  stop("combined IBD file must contain a 'category' column", call. = FALSE)
}
category_order <- combined_ibd$category |> unique() |> sort()

# Save a clean fraction table (for downstream)
write_tsv(fraction_ibd, file.path(outdir, "fraction_ibd.tsv"))

# ─────────────────────────────────────────────────────────────
# Process FAI (chromosome coordinates)
# ─────────────────────────────────────────────────────────────

fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  rename(chr = V1, end_chr = V2) %>%
  filter(!is.na(chr))

if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
  rm_chr <- strsplit(rm_chr_str, ",")[[1]] |> trimws()
  fai <- fai %>% filter(!chr %in% rm_chr)
}

fai <- fai %>%
  mutate(chr_num = as.numeric(stringr::str_match(chr, pattern)[, groupid])) %>%
  arrange(chr_num) %>%
  mutate(
    tr_chr    = lag(cumsum(end_chr), default = 0),
    fill_color = rep(c("grey95", "white"), length.out = n())
  )

transpose_chr <- fai %>% select(chr = chr_num, tr_chr)
chr_labels    <- fai %>% mutate(chr = chr_num, mid = tr_chr + end_chr / 2)

# ─────────────────────────────────────────────────────────────
# Boxplot of IBD fractions
# ─────────────────────────────────────────────────────────────

fraction_ibd_plot <- fraction_ibd %>%
  filter(!is.na(fraction), !is.na(category)) %>%
  mutate(category = factor(category, levels = category_order))

g_box <- ggplot(fraction_ibd_plot,
                aes(x = category, y = fraction, fill = category)) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 1,
               color = "yellow") +
  theme_classic(base_size = 12) +
  ylim(0, 0.5) +
  labs(x = "", y = "Pairwise fraction IBD") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title  = element_text(size = 30),
    legend.position = "none"
  )

ggsave(file.path(outdir, "ibd_fraction_boxplot.png"),
       g_box, width = 8, height = 6)

# ─────────────────────────────────────────────────────────────
# Combine IBD + genomic position
# ─────────────────────────────────────────────────────────────

combined_ibd_tr <- combined_ibd %>%
  filter(!is.na(fraction),
         !is.na(category),
         !is.na(chr),
         !is.na(win_start)) %>%
  mutate(chr = as.numeric(chr)) %>%
  left_join(transpose_chr, by = "chr") %>%
  mutate(
    pos_bp_ed = as.numeric(win_start) + tr_chr,
    category  = factor(category, levels = category_order)
  )

signif_cutoff <- quantile(combined_ibd_tr$fraction,
                          probs = th_quantile, na.rm = TRUE)
message("95th %-tile (or chosen quantile) of fraction: ", signif_cutoff)

# ─────────────────────────────────────────────────────────────
# Gene annotation (for high-IBD windows and drug genes)
# ─────────────────────────────────────────────────────────────

genes <- NULL
genes_of_interest <- NULL
gene_labels_clean <- NULL

if (!is.null(gene_file) && file.exists(gene_file)) {
  genes <- read_tsv(gene_file, col_types = cols()) %>%
    rename(
      gene_chr_raw = chr,
      gene_start   = pos_start,
      gene_end     = pos_end,
      gene_id      = gene_id,
      product      = gene_product
    ) %>%
    mutate(
      gene_chr = as.numeric(str_extract(gene_chr_raw, "(?<=Pf3D7_)\\d{2}"))
    )

  # Identify high-IBD windows overlapping genes
  ibd_outliers <- combined_ibd_tr %>%
    filter(fraction > signif_cutoff) %>%
    select(chr, win_start, win_end, fraction, category)

  genes_of_interest <- inner_join(
    ibd_outliers, genes,
    by = c("chr" = "gene_chr")
  ) %>%
    filter(win_start <= gene_end & win_end >= gene_start)

  write_tsv(genes_of_interest,
            file.path(outdir, "ibd_selection_candidate_genes.tsv"))

  # Only label drug-resistance genes
  res_genes <- c("Kelch13", "CRT", "DHFR", "PPPK-DHPS", "MDR1",
                 "UBP1", "coronin", "AP2-MU", "PX1")

  gene_labels_clean <- genes_of_interest %>%
    mutate(chr = as.numeric(chr)) %>%
    left_join(transpose_chr, by = "chr") %>%
    mutate(
      label_pos = (gene_start + gene_end) / 2 + tr_chr,
      resistance_flag = if_else(
        str_detect(product,
                   regex(paste(res_genes, collapse = "|"), ignore_case = TRUE)),
        "resistance", "other"
      )
    ) %>%
    filter(resistance_flag == "resistance") %>%
    distinct(gene_id, label_pos, category)
} else {
  message("No gene_product file given or file not found; skipping gene-based annotation.")
}

# ─────────────────────────────────────────────────────────────
# Genome-wide IBD fraction plot (per category)
# ─────────────────────────────────────────────────────────────

p_clean <- ggplot(combined_ibd_tr) +
  geom_rect(
    data = fai,
    aes(xmin = tr_chr,
        xmax = tr_chr + end_chr,
        ymin = -Inf, ymax = Inf,
        fill = fill_color),
    inherit.aes = FALSE,
    alpha = 0.5, color = NA
  ) +
  scale_fill_identity() +
  geom_line(aes(x = pos_bp_ed, y = fraction, color = category),
            size = 0.8) +
  facet_wrap(~ category, ncol = 1, scales = "free_x") +
  scale_x_continuous(breaks = chr_labels$mid, labels = chr_labels$chr) +
  scale_y_continuous(limits = c(0, 0.5)) +
  geom_hline(yintercept = signif_cutoff,
             linetype = "dashed", color = "red") +
  labs(x = "Chromosome", y = "IBD Fraction", color = "Category") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 15) +
  theme(
    strip.text       = element_blank(),
    strip.background = element_blank(),
    axis.text.x      = element_text(size = 12),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(size = 17),
    legend.title     = element_text(size = 15),
    legend.text      = element_text(size = 14),
    legend.position  = "bottom"
  )

ggsave(file.path(outdir, "ibd_genomewide_fraction_cleaned.png"),
       p_clean, width = 12, height = 15)

# ─────────────────────────────────────────────────────────────
# Chromosome painting of high-IBD segments
# (with drug-resistance genes highlighted in red)
# ─────────────────────────────────────────────────────────────

# Collapse top-IBD windows into contiguous segments per category & chr
ibd_top <- combined_ibd_tr %>%
  filter(fraction > signif_cutoff) %>%
  transmute(
    category,
    chr   = as.integer(chr),
    start = as.numeric(win_start),
    end   = as.numeric(win_end)
  )

collapse_intervals <- function(df) {
  if (nrow(df) == 0) return(df)
  setDT(df)[order(category, chr, start, end)]
  df[, lag_end := shift(end, type = "lag"), by = .(category, chr)]
  df[, grp := cumsum(ifelse(is.na(lag_end) | start > lag_end, 1L, 0L)),
     by = .(category, chr)]
  df[, .(start = min(start), end = max(end)),
     by = .(category, chr, grp)][, grp := NULL][]
}

ibd_segments <- collapse_intervals(ibd_top) %>% as_tibble()

# “Pill” scaffolding
pill_x   <- tibble(chr = 1:14, x = chr)
chr_lens <- fai %>% transmute(chr = chr_num, len_bp = end_chr)
to_mb    <- function(x) x / 1e6

pill_backbone <- expand_grid(
  category = factor(category_order, levels = category_order),
  pill_x
) %>%
  left_join(chr_lens, by = "chr") %>%
  mutate(y0 = 0, y1 = to_mb(len_bp))

paint_segs <- ibd_segments %>%
  left_join(pill_x, by = "chr") %>%
  mutate(y0 = to_mb(start),
         y1 = to_mb(end))

# Drug-resistance gene bars for chromosome painting
dr_paint <- NULL
if (!is.null(genes_of_interest) && nrow(genes_of_interest) > 0) {
  res_genes <- c("Kelch13", "CRT", "DHFR", "PPPK-DHPS", "MDR1",
                 "UBP1", "coronin", "AP2-MU", "PX1")

  dr_paint <- genes_of_interest %>%
    filter(str_detect(product,
                      regex(paste(res_genes, collapse = "|"),
                            ignore_case = TRUE))) %>%
    mutate(chr = as.integer(chr)) %>%
    left_join(pill_x, by = "chr") %>%
    transmute(
      category,
      chr,
      x,
      y0 = gene_start / 1e6,
      y1 = gene_end   / 1e6
    )
}

region_cols <- scales::hue_pal()(length(category_order))
names(region_cols) <- category_order

paint_plot <- ggplot() +
  # chromosome outline
  geom_segment(
    data = pill_backbone,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 8.5, lineend = "round", color = "grey30"
  ) +
  # chromosome body
  geom_segment(
    data = pill_backbone,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 6.0, lineend = "round", color = "white"
  ) +
  # high-IBD segments
  geom_segment(
    data = paint_segs,
    aes(x = x, xend = x, y = y0, yend = y1, color = category),
    linewidth = 5, lineend = "butt", alpha = 0.9
  ) +
  # drug-resistance genes in red
  { if (!is.null(dr_paint) && nrow(dr_paint) > 0)
      geom_segment(
        data = dr_paint,
        aes(x = x, xend = x, y = y0, yend = y1),
        linewidth = 3, lineend = "butt", color = "red"
      )
    else NULL } +
  scale_color_manual(values = region_cols, guide = "none") +
  facet_wrap(~ factor(category, levels = category_order),
             ncol = 1, strip.position = "left") +
  geom_text(
    data = pill_backbone %>% group_by(category, chr, x) %>% slice_tail(),
    aes(x = x, y = y1 + 0.4, label = chr),
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, max(pill_backbone$y1) + 0.6),
    expand = expansion(mult = c(0, 0.005))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 12) +
  theme(
    strip.background   = element_blank(),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(size = 15, face = "bold",
                                      margin = margin(r = 0), angle = 0),
    aspect.ratio       = 0.58,
    panel.spacing      = unit(8, "pt"),
    plot.margin        = margin(5, 10, 5, 10)
  )

ggsave(file.path(outdir, "ibd_chromosome_painting.png"),
       paint_plot, width = 8, height = 13, dpi = 300, bg = "white")

message("Done plotting IBD.")


