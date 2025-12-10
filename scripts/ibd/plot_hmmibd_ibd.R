#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
  library(stringr)
  library(ggrepel)
  library(scales)
})

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Directory containing *_hmmIBD_ibd_winXkb.tsv and *_hmmIBD_fraction.tsv [default %default]",
              metavar = "character"),
  make_option(c("--ref_index"), type = "character", default = NULL,
              help = "FAI index for reference genome (required)",
              metavar = "character"),
  make_option(c("--gene_product"), type = "character", default = NULL,
              help = "Gene annotation TSV (pf_genome_product_v3.tsv). Used for candidate gene table.",
              metavar = "character"),
  make_option(c("--suffix"), type = "character", default = NULL,
              help = "Prefix used by summarise_hmmibd_windows.R (e.g. 13_08_2025). If NULL, first *_hmmIBD_ibd_*.tsv is used to infer.",
              metavar = "character"),
  make_option(c("--window_size"), type = "numeric", default = 50000,
              help = "Sliding window size in bp that was used in summarise step [default %default]",
              metavar = "numeric"),
  make_option(c("--quantile_cutoff"), type = "numeric", default = 0.95,
              help = "Quantile cutoff for high-IBD windows [default %default]",
              metavar = "numeric"),
  make_option(c("--remove_chr"), type = "character",
              default = "Pf3D7_API_v3,Pf3D7_MIT_v3",
              help = "Comma-separated list of chromosomes to drop [default %default]",
              metavar = "character"),
  make_option(c("--regex_chr"), type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex to derive numeric chromosome from chr name [default %default]",
              metavar = "character"),
  make_option(c("--regex_groupid"), type = "integer", default = 3,
              help = "Capture group index in regex_chr giving numeric chromosome [default %default]",
              metavar = "integer"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for plots (default: workdir/win_Xkb)",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ─────────────────────────────────────────────────────────────
# Extract options
# ─────────────────────────────────────────────────────────────

workdir        <- opt$workdir
ref_index      <- opt$ref_index
gene_file      <- opt$gene_product
suffix         <- opt$suffix
window_size    <- opt$window_size
th_quantile    <- opt$quantile_cutoff
rm_chr_str     <- opt$remove_chr
pattern        <- opt$regex_chr
groupid        <- opt$regex_groupid
outdir         <- opt$outdir

if (is.null(ref_index)) {
  stop("ERROR: --ref_index is required\n", call. = FALSE)
}

if (is.null(outdir)) {
  outdir <- file.path(workdir, sprintf("win_%dkb", window_size / 1000))
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Workdir: ", workdir)
message("Outdir: ", outdir)

# ─────────────────────────────────────────────────────────────
# Locate summarised IBD files
# ─────────────────────────────────────────────────────────────

if (is.null(suffix)) {
  # Try to infer suffix from first *_hmmIBD_ibd_winXkb.tsv
  pattern_ibd <- sprintf(".*_hmmIBD_ibd_win%dkb.tsv$", window_size / 1000)
  files_ibd <- list.files(workdir, pattern = pattern_ibd, full.names = FALSE)
  if (length(files_ibd) == 0) {
    stop("Could not find any *_hmmIBD_ibd_win", window_size/1000,
         "kb.tsv in workdir. Please provide --suffix.\n", call. = FALSE)
  }
  suffix <- sub("_hmmIBD_ibd_win.*$", "", files_ibd[1])
  message("Inferred suffix: ", suffix)
}

ibd_file   <- file.path(workdir,
                        sprintf("%s_hmmIBD_ibd_win%dkb.tsv", suffix, window_size / 1000))
fract_file <- file.path(workdir,
                        sprintf("%s_hmmIBD_fraction.tsv", suffix))

if (!file.exists(ibd_file)) {
  stop("IBD window file not found: ", ibd_file, call. = FALSE)
}
if (!file.exists(fract_file)) {
  stop("Fraction file not found: ", fract_file, call. = FALSE)
}

message("IBD windows: ", ibd_file)
message("Fractions:   ", fract_file)

combined_ibd <- readr::read_tsv(ibd_file, col_types = cols())
fraction_ibd <- readr::read_tsv(fract_file, col_types = cols())

# Coerce category to character
combined_ibd$category <- as.character(combined_ibd$category)
fraction_ibd$category <- as.character(fraction_ibd$category)

category_order <- combined_ibd$category %>%
  unique() %>% sort()

# ─────────────────────────────────────────────────────────────
# Reference index (FAI)
# ─────────────────────────────────────────────────────────────

fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  dplyr::rename(chr = V1, end_chr = V2) %>%
  dplyr::select(chr, end_chr)

if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
  rm_chr_vec <- strsplit(rm_chr_str, ",")[[1]] |> trimws()
  present_rm <- intersect(rm_chr_vec, unique(fai$chr))
  if (length(present_rm) > 0) {
    fai <- fai %>% dplyr::filter(!chr %in% present_rm)
  } else {
    warning("None of the --remove_chr names were found in FAI; continuing.\n")
  }
}

fai <- fai %>%
  mutate(chr_num = as.numeric(str_match(chr, pattern)[, groupid])) %>%
  arrange(chr_num) %>%
  mutate(
    tr_chr     = lag(cumsum(end_chr), default = 0),
    fill_color = rep(c("grey95", "white"), length.out = n())
  )

transpose_chr <- fai %>% select(chr = chr_num, tr_chr)
chr_labels    <- fai %>% mutate(chr = chr_num, mid = tr_chr + end_chr / 2)

# ─────────────────────────────────────────────────────────────
# IBD fraction boxplot (with nicer font sizes)
# ─────────────────────────────────────────────────────────────

fraction_ibd_plot <- fraction_ibd %>%
  filter(!is.na(fraction), !is.na(category)) %>%
  mutate(category = factor(category, levels = category_order))

g_box <- ggplot(fraction_ibd_plot,
                aes(x = category, y = fraction, fill = category)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.6) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 2.5, fill = "yellow", color = "black") +
  ylim(0, 0.5) +
  labs(x = NULL, y = "Pairwise fraction IBD") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 13),
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  filename = file.path(outdir, "ibd_fraction_boxplot.png"),
  plot     = g_box,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# ─────────────────────────────────────────────────────────────
# Genome-wide IBD fraction line plot
# ─────────────────────────────────────────────────────────────

combined_ibd_tr <- combined_ibd %>%
  filter(!is.na(fraction),
         !is.na(category),
         !is.na(chr),
         !is.na(win_start)) %>%
  mutate(
    chr      = as.numeric(chr),
    category = factor(category, levels = category_order)
  ) %>%
  left_join(transpose_chr, by = "chr") %>%
  mutate(pos_bp_ed = as.numeric(win_start) + tr_chr)

signif_cutoff <- quantile(combined_ibd_tr$fraction,
                          probs = th_quantile, na.rm = TRUE)

p_clean <- ggplot(combined_ibd_tr) +
  geom_rect(
    data = fai,
    aes(xmin = tr_chr, xmax = tr_chr + end_chr,
        ymin = -Inf, ymax = Inf, fill = fill_color),
    inherit.aes = FALSE, alpha = 0.5, color = NA
  ) +
  scale_fill_identity() +
  geom_line(
    aes(x = pos_bp_ed, y = fraction, color = category),
    linewidth = 0.6
  ) +
  facet_wrap(~ category, ncol = 1, scales = "free_x") +
  scale_x_continuous(
    breaks = chr_labels$mid,
    labels = chr_labels$chr
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  geom_hline(
    yintercept = signif_cutoff,
    linetype   = "dashed",
    color      = "red"
  ) +
  labs(x = "Chromosome", y = "IBD fraction", color = "Category") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 11) +
  theme(
    strip.text       = element_blank(),
    strip.background = element_blank(),
    axis.text.x      = element_text(size = 8),
    axis.text.y      = element_text(size = 9),
    axis.title       = element_text(size = 13),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    legend.position  = "bottom"
  )

ggsave(
  filename = file.path(outdir, "ibd_genomewide_fraction.png"),
  plot     = p_clean,
  width    = 10,
  height   = 10,
  dpi      = 300
)

# ─────────────────────────────────────────────────────────────
# Optional: gene overlap + candidate list
# ─────────────────────────────────────────────────────────────

if (!is.null(gene_file) && file.exists(gene_file)) {
  genes <- readr::read_tsv(gene_file, col_types = cols()) %>%
    rename(
      gene_chr_raw = chr,
      gene_start   = pos_start,
      gene_end     = pos_end,
      gene_id      = gene_id,
      product      = gene_product
    ) %>%
    mutate(gene_chr = as.numeric(str_extract(gene_chr_raw, "(?<=Pf3D7_)\\d{2}")))

  ibd_outliers <- combined_ibd_tr %>%
    filter(fraction > signif_cutoff) %>%
    select(chr, win_start, win_end, fraction, category)

  genes_of_interest <- inner_join(
    ibd_outliers,
    genes,
    by = c("chr" = "gene_chr")
  ) %>%
    filter(win_start <= gene_end & win_end >= gene_start)

  readr::write_tsv(
    genes_of_interest,
    file.path(outdir, "ibd_selection_candidate_genes.tsv")
  )
} else {
  message("Gene annotation file not provided or not found; skipping candidate gene table.")
}

# ─────────────────────────────────────────────────────────────
# Chromosome painting with high-IBD segments + drug genes
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
  df[, .(start = min(start), end = max(end)), by = .(category, chr, grp)][, grp := NULL][]
}

ibd_segments <- collapse_intervals(ibd_top) %>% as_tibble()

# Chromosome “pill” scaffolding for each facet
pill_x <- tibble(chr = 1:14, x = chr)

chr_lens <- fai %>% transmute(chr = chr_num, len_bp = end_chr)
to_mb <- function(x) x / 1e6

pill_backbone <- expand_grid(
  category = factor(category_order, levels = category_order),
  pill_x
) %>%
  left_join(chr_lens, by = "chr") %>%
  mutate(
    y0 = 0,
    y1 = to_mb(len_bp)
  )

paint_segs <- ibd_segments %>%
  left_join(pill_x, by = "chr") %>%
  mutate(
    y0 = to_mb(start),
    y1 = to_mb(end)
  )

# Drug-resistance genes (fixed genomic coordinates)
drug_genes <- tibble::tribble(
  ~gene,     ~chr, ~start,   ~end,
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

drug_df <- drug_genes %>%
  mutate(
    chr = as.integer(chr),
    y0  = to_mb(start),
    y1  = to_mb(end)
  ) %>%
  left_join(pill_x, by = "chr")

# replicate drug bars across all categories
drug_segs <- expand_grid(
  category = factor(category_order, levels = category_order),
  drug_df
)

region_cols <- hue_pal()(length(category_order))
names(region_cols) <- category_order

paint_plot <- ggplot() +
  # chromosome outline
  geom_segment(
    data = pill_backbone,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 8.5,
    lineend   = "round",
    color     = "grey30"
  ) +
  # chromosome body
  geom_segment(
    data = pill_backbone,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 6.0,
    lineend   = "round",
    color     = "white"
  ) +
  # drug-resistance gene bars (red, constant per facet)
  geom_segment(
    data = drug_segs,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 4,
    lineend   = "butt",
    color     = "red",
    alpha     = 0.9
  ) +
  # high-IBD segments (per category)
  geom_segment(
    data = paint_segs,
    aes(x = x, xend = x, y = y0, yend = y1, color = category),
    linewidth = 5,
    lineend   = "butt",
    alpha     = 0.9
  ) +
  scale_color_manual(values = region_cols, guide = "none") +
  facet_wrap(
    ~ factor(category, levels = category_order),
    ncol = 1,
    strip.position = "left"
  ) +
  geom_text(
    data = pill_backbone %>% group_by(category, chr, x) %>% slice_tail(),
    aes(x = x, y = y1 + 0.4, label = chr),
    size = 4
  ) +
  scale_y_continuous(
    limits = c(0, max(pill_backbone$y1) + 0.6),
    expand = expansion(mult = c(0, 0.005))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(
    strip.background   = element_blank(),
    strip.placement    = "outside",
    strip.text.y.left  = element_text(
      size   = 12,
      face   = "bold",
      margin = margin(r = 0),
      angle  = 0
    ),
    aspect.ratio       = 0.6,
    panel.spacing      = unit(6, "pt"),
    plot.margin        = margin(5, 10, 5, 10)
  )

ggsave(
  filename = file.path(outdir, "ibd_chromosome_painting.png"),
  plot     = paint_plot,
  width    = 7,
  height   = 11,
  dpi      = 300,
  bg       = "white"
)

message("Done.")
