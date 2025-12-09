#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)
  library(tidyverse)
  library(stringr)
  library(data.table)
  library(ggrepel)
})

# ─────────────────────────────────────────────────────────────
# CLI options
# ─────────────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Directory containing the *_hmmIBD_ibd_win50kb.tsv and *_hmmIBD_fraction.tsv files [default %default]",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Metadata TSV with grouping variable [required]",
              metavar = "character"),
  make_option(c("-f", "--ref_fai"), type = "character", default = NULL,
              help = "FAI index of Pf3D7 reference genome [required]",
              metavar = "character"),
  make_option(c("-g", "--gene_file"), type = "character", default = NULL,
              help = "Pf3D7 genome product TSV annotation [required]",
              metavar = "character"),
  make_option(c("-s", "--suffix"), type = "character", default = NULL,
              help = "Prefix used for hmmIBD summary files, e.g. '13_08_2025' if files are 13_08_2025_hmmIBD_ibd_win50kb.tsv and 13_08_2025_hmmIBD_fraction.tsv [required]",
              metavar = "character"),
  make_option(c("--group_var"), type = "character", default = "collection",
              help = "Metadata column used as 'category' (e.g. collection, year, region) [default %default]",
              metavar = "character"),
  make_option(c("--remove_chr"), type = "character",
              default = "Pf3D7_API_v3,Pf3D7_MIT_v3",
              help = "Comma-separated chromosomes to remove from FAI [default %default]",
              metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = "win_50kb",
              help = "Subdirectory inside workdir for plots and tables [default %default]",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$metadata) || is.null(opt$ref_fai) ||
    is.null(opt$gene_file) || is.null(opt$suffix)) {
  stop("ERROR: --metadata, --ref_fai, --gene_file, and --suffix are all required.\n",
       call. = FALSE)
}

workdir     <- opt$workdir
metadata_file <- opt$metadata
ref_index   <- opt$ref_fai
gene_file   <- opt$gene_file
suffix      <- opt$suffix
group_var   <- opt$group_var
rm_chr_str  <- opt$remove_chr
output_dir  <- opt$output_dir

setwd(workdir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Workdir:     ", workdir, "\n")
cat("Suffix:      ", suffix, "\n")
cat("Metadata:    ", metadata_file, "\n")
cat("Ref FAI:     ", ref_index, "\n")
cat("Gene file:   ", gene_file, "\n")
cat("Group var:   ", group_var, "\n")
cat("Remove chr:  ", rm_chr_str, "\n")
cat("Output dir:  ", file.path(workdir, output_dir), "\n\n")

# ─────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────
metadata <- read_tsv(metadata_file, col_types = cols())

# Order of categories (e.g. years / collections)
category_order <- metadata %>%
  pull(!!sym(group_var)) %>%
  unique() %>%
  sort()

ibd_file      <- sprintf("%s_hmmIBD_ibd_win50kb.tsv", suffix)
fraction_file <- sprintf("%s_hmmIBD_fraction.tsv", suffix)

if (!file.exists(ibd_file)) {
  stop("ERROR: Could not find combined IBD file: ", ibd_file, call. = FALSE)
}
if (!file.exists(fraction_file)) {
  stop("ERROR: Could not find fraction IBD file: ", fraction_file, call. = FALSE)
}

combined_ibd <- read_tsv(ibd_file, col_types = cols())
fraction_ibd <- read_tsv(fraction_file, col_types = cols())

# Attach category as factor with desired order (if already present)
if ("category" %in% colnames(combined_ibd)) {
  combined_ibd <- combined_ibd %>%
    mutate(category = factor(category, levels = category_order))
}
if ("category" %in% colnames(fraction_ibd)) {
  fraction_ibd <- fraction_ibd %>%
    mutate(category = factor(category, levels = category_order))
}

# For simplicity, assume combined_ibd and fraction_ibd already have "category"
# from your prior summarisation pipeline. We keep the left_join logic minimal.
metadata_sub <- metadata %>%
  select(category = !!sym(group_var)) %>%
  distinct()

combined_ibd <- combined_ibd %>%
  left_join(metadata_sub, by = "category")
fraction_ibd <- fraction_ibd %>%
  left_join(metadata_sub, by = "category")

write_tsv(fraction_ibd, file.path(output_dir, "fraction_ibd.tsv"))

# ─────────────────────────────────────────────────────────────
# Process FAI (chromosome coordinates)
# ─────────────────────────────────────────────────────────────
rm_chr <- if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
  strsplit(rm_chr_str, ",")[[1]] |> trimws()
} else {
  character(0)
}

pattern <- "(.*?)_(.+)_(.*)"
groupid <- 3

fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  dplyr::rename(chr = V1, end_chr = V2) %>%
  dplyr::filter(!chr %in% rm_chr) %>%
  mutate(chr_num = as.numeric(stringr::str_match(chr, pattern)[, groupid])) %>%
  arrange(chr_num) %>%
  mutate(tr_chr = dplyr::lag(cumsum(end_chr), default = 0),
         fill_color = rep(c("grey95", "white"), length.out = dplyr::n()))

transpose_chr <- fai %>%
  select(chr = chr_num, tr_chr)

chr_labels <- fai %>%
  mutate(chr = chr_num,
         mid = tr_chr + end_chr / 2)

# ─────────────────────────────────────────────────────────────
# BOXPLOT OF IBD FRACTIONS
# ─────────────────────────────────────────────────────────────
fraction_ibd_plot <- fraction_ibd %>%
  filter(!is.na(fraction), !is.na(category)) %>%
  mutate(category = factor(category, levels = category_order))

g_box <- ggplot(fraction_ibd_plot, aes(x = category, y = fraction, fill = category)) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 9, size = 1, color = "yellow") +
  theme_classic(base_size = 12) +
  ylim(0, 0.5) +
  labs(x = "", y = "Pairwise fraction IBD") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title  = element_text(size = 30),
    legend.position = "none"
  )

ggsave(file.path(output_dir, "ibd_fraction_boxplot.png"),
       g_box, width = 8, height = 6)

# ─────────────────────────────────────────────────────────────
# GENOME-WIDE IBD FRACTION PLOT
# ─────────────────────────────────────────────────────────────
combined_ibd_tr <- combined_ibd %>%
  filter(!is.na(fraction), !is.na(category),
         !is.na(chr), !is.na(win_start)) %>%
  mutate(
    chr = as.numeric(chr),
    category = factor(category, levels = category_order)
  ) %>%
  left_join(transpose_chr, by = "chr") %>%
  mutate(pos_bp_ed = as.numeric(win_start) + tr_chr)

signif_cutoff <- quantile(combined_ibd_tr$fraction, probs = 0.95, na.rm = TRUE)

# Load gene annotation (for candidate gene table)
genes <- read_tsv(gene_file, col_types = cols()) %>%
  rename(
    gene_chr_raw = chr,
    gene_start   = pos_start,
    gene_end     = pos_end,
    gene_id      = gene_id,
    product      = gene_product
  ) %>%
  mutate(gene_chr = as.numeric(str_extract(gene_chr_raw, "(?<=Pf3D7_)\\d{2}")))

# High-IBD windows overlapping genes
ibd_outliers <- combined_ibd_tr %>%
  filter(fraction > signif_cutoff) %>%
  select(chr, win_start, win_end, fraction, category)

genes_of_interest <- inner_join(
  ibd_outliers,
  genes,
  by = c("chr" = "gene_chr")
) %>%
  filter(win_start <= gene_end & win_end >= gene_start)

write_tsv(genes_of_interest,
          file.path(output_dir, "ibd_selection_candidate_genes.tsv"))

# Drug-resistance genes for labeling
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

# Genome-wide fraction plot
p_clean <- ggplot(combined_ibd_tr) +
  geom_rect(
    data = fai,
    aes(xmin = tr_chr, xmax = tr_chr + end_chr,
        ymin = -Inf, ymax = Inf, fill = fill_color),
    inherit.aes = FALSE, alpha = 0.5, color = NA
  ) +
  scale_fill_identity() +
  geom_line(aes(x = pos_bp_ed, y = fraction, color = category), size = 0.8) +
  facet_wrap(~ category, ncol = 1, scales = "free_x") +
  scale_x_continuous(
    breaks = chr_labels$mid,
    labels = chr_labels$chr
  ) +
  scale_y_continuous(limits = c(0, 0.5)) +
  geom_hline(yintercept = signif_cutoff, linetype = "dashed", color = "red") +
  labs(x = "Chromosome", y = "IBD Fraction", color = "Category") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 15) +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 17),
    legend.title = element_text(size = 15),
    legend.text  = element_text(size = 14),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "ibd_genomewide_fraction_cleaned.png"),
       p_clean, width = 12, height = 15)

# ─────────────────────────────────────────────────────────────
# CHROMOSOME PAINTING PLOT
# ─────────────────────────────────────────────────────────────
# Collapse top-IBD windows into contiguous segments
ibd_top <- combined_ibd_tr %>%
  filter(fraction > signif_cutoff) %>%
  transmute(category,
            chr = as.integer(chr),
            start = as.numeric(win_start),
            end   = as.numeric(win_end))

collapse_intervals <- function(df) {
  if (nrow(df) == 0) return(df)
  dt <- as.data.table(df)
  setorder(dt, category, chr, start, end)
  dt[, lag_end := shift(end, type = "lag"), by = .(category, chr)]
  dt[, grp := cumsum(ifelse(is.na(lag_end) | start > lag_end, 1L, 0L)),
     by = .(category, chr)]
  dt[, .(start = min(start), end = max(end)), by = .(category, chr, grp)][
    , grp := NULL][]
}

ibd_segments <- collapse_intervals(ibd_top) %>% as_tibble()

# Chromosome “pill” scaffold
pill_x <- tibble(chr = 1:14, x = chr)

chr_lens <- fai %>%
  transmute(chr = chr_num, len_bp = end_chr)
to_mb <- function(x) x / 1e6

pill_backbone <- tidyr::expand_grid(
  category = factor(category_order, levels = category_order),
  pill_x
) %>%
  left_join(chr_lens, by = "chr") %>%
  mutate(y0 = 0, y1 = to_mb(len_bp))

paint_segs <- ibd_segments %>%
  left_join(pill_x, by = "chr") %>%
  mutate(
    y0 = to_mb(start),
    y1 = to_mb(end)
  )

region_cols <- scales::hue_pal()(length(category_order))
names(region_cols) <- category_order

paint_plot <- ggplot() +
  # chromosome outline (outer)
  geom_segment(
    data = pill_backbone,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 8.5, lineend = "round", color = "grey30"
  ) +
  # chromosome body (inner)
  geom_segment(
    data = pill_backbone,
    aes(x = x, xend = x, y = y0, yend = y1),
    linewidth = 6.0, lineend = "round", color = "white"
  ) +
  # painted high-IBD segments
  geom_segment(
    data = paint_segs,
    aes(x = x, xend = x, y = y0, yend = y1, color = category),
    linewidth = 5, lineend = "butt", alpha = 0.9
  ) +
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
    strip.text.y.left  = element_text(
      size = 15, face = "bold",
      margin = margin(r = 0), angle = 0
    ),
    aspect.ratio = 0.58,
    panel.spacing = unit(8, "pt"),
    plot.margin = margin(5, 10, 5, 10)
  )

ggsave(file.path(output_dir, "ibd_chromosome_painting.png"),
       paint_plot, width = 8, height = 13, dpi = 300, bg = "white")

cat("\nAll IBD plots written to: ", file.path(workdir, output_dir), "\n")
