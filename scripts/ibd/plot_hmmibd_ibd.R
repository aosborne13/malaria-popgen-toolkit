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
              help = "Directory containing *_hmmIBD_ibd_winXkb.tsv / *_hmmIBD_fraction.tsv",
              metavar = "character"),
  make_option(c("--ref_index"), type = "character", default = NULL,
              help = "FAI index for reference genome (required)",
              metavar = "character"),
  make_option(c("--gene_product"), type = "character", default = NULL,
              help = "Gene product annotation TSV (pf_genome_product_v3.tsv)",
              metavar = "character"),
  make_option(c("--suffix"), type = "character", default = NULL,
              help = "Prefix used by summarise_hmmibd_windows.R (e.g. 13_08_2025)",
              metavar = "character"),
  make_option(c("--window_size"), type = "numeric", default = 50000,
              help = "Sliding window size in bp used in summary [default %default]",
              metavar = "numeric"),
  make_option(c("--quantile_cutoff"), type = "numeric", default = 0.95,
              help = "Quantile cutoff for high-IBD windows [default %default]",
              metavar = "numeric"),
  make_option(c("--remove_chr"), type = "character",
              default = "Pf3D7_API_v3,Pf3D7_MIT_v3",
              help = "Comma-separated chromosomes to drop [default %default]",
              metavar = "character"),
  make_option(c("--regex_chr"), type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex for numeric chromosome from chr name [default %default]",
              metavar = "character"),
  make_option(c("--regex_groupid"), type = "integer", default = 3,
              help = "Capture group index in regex_chr giving numeric chromosome [default %default]",
              metavar = "integer"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory for plots (default: workdir/win_Xkb)",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

workdir         <- opt$workdir
ref_index       <- opt$ref_index
gene_file       <- opt$gene_product
suffix          <- opt$suffix
window_size     <- opt$window_size
quantile_cutoff <- opt$quantile_cutoff
rm_chr_str      <- opt$remove_chr
pattern         <- opt$regex_chr
groupid         <- opt$regex_groupid

if (is.null(suffix)) {
  stop("ERROR: --suffix is required (must match summarise_hmmibd_windows.R outputs).",
       call. = FALSE)
}
if (is.null(ref_index) || !file.exists(ref_index)) {
  stop("ERROR: --ref_index not provided or file does not exist.", call. = FALSE)
}
if (is.null(gene_file) || !file.exists(gene_file)) {
  stop("ERROR: --gene_product not provided or file does not exist.", call. = FALSE)
}

window_kb <- window_size / 1000

if (is.null(opt$outdir)) {
  outdir <- file.path(workdir, sprintf("win_%dkb", window_kb))
} else {
  outdir <- opt$outdir
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Workdir: ", workdir)
message("Output dir: ", outdir)

# ─────────────────────────────────────────────────────────────
# Load summary + annotation files
# ─────────────────────────────────────────────────────────────

ibd_file <- file.path(
  workdir,
  sprintf("%s_hmmIBD_ibd_win%dkb.tsv", suffix, window_kb)
)
frac_file <- file.path(
  workdir,
  sprintf("%s_hmmIBD_fraction.tsv", suffix)
)

# Robust handling of annotated filename: support both q0.90 and q0.9
q_str_full  <- sprintf("%.2f", quantile_cutoff)          # e.g. "0.90"
q_str_short <- sub("0+$", "", q_str_full)                # e.g. "0.9"

annot_candidates <- c(
  file.path(
    workdir,
    sprintf("%s_hmmIBD_ibd_win%dkb_annotated_q%s.tsv",
            suffix, window_kb, q_str_full)
  ),
  file.path(
    workdir,
    sprintf("%s_hmmIBD_ibd_win%dkb_annotated_q%s.tsv",
            suffix, window_kb, q_str_short)
  )
)

annot_file_existing <- annot_candidates[file.exists(annot_candidates)]
annot_file <- if (length(annot_file_existing) > 0) annot_file_existing[1] else NA

if (!file.exists(ibd_file)) {
  stop("Cannot find IBD window file: ", ibd_file, call. = FALSE)
}
if (!file.exists(frac_file)) {
  stop("Cannot find fraction file: ", frac_file, call. = FALSE)
}
if (is.na(annot_file)) {
  stop(
    "Cannot find annotated IBD windows file.\n",
    "Tried:\n  ",
    paste(annot_candidates, collapse = "\n  "),
    "\nMake sure summarise_hmmibd_windows.R has been run and --quantile_cutoff matches.",
    call. = FALSE
  )
}

message("Using annotated file: ", annot_file)

combined_ibd <- readr::read_tsv(ibd_file, col_types = cols())
fraction_ibd <- readr::read_tsv(frac_file, col_types = cols())
annot_ibd    <- readr::read_tsv(annot_file, col_types = cols())

# Expect columns: chr, win_start, win_end, category, fraction
needed_ibd_cols  <- c("chr", "win_start", "win_end", "category", "fraction")
missing_ibd_cols <- setdiff(needed_ibd_cols, colnames(combined_ibd))
if (length(missing_ibd_cols) > 0) {
  stop("Missing columns in combined IBD file: ",
       paste(missing_ibd_cols, collapse = ", "), call. = FALSE)
}

# category order
category_order <- combined_ibd %>%
  filter(!is.na(category)) %>%
  pull(category) %>%
  unique() %>%
  sort()

# ─────────────────────────────────────────────────────────────
# Reference index → chromosome coordinates
# ─────────────────────────────────────────────────────────────

fai <- read.table(ref_index, stringsAsFactors = FALSE) %>%
  dplyr::rename(chr = V1, end_chr = V2) %>%
  dplyr::select(chr, end_chr)

if (!is.null(rm_chr_str) && nzchar(rm_chr_str)) {
  rm_chr <- strsplit(rm_chr_str, ",")[[1]] |> trimws()
  fai <- fai %>% dplyr::filter(!chr %in% rm_chr)
}

fai <- fai %>%
  mutate(chr_num = as.numeric(stringr::str_match(chr, pattern)[, groupid])) %>%
  arrange(chr_num) %>%
  mutate(tr_chr    = dplyr::lag(cumsum(end_chr), default = 0),
         fill_color = rep(c("grey95", "white"), length.out = n()))

transpose_chr <- fai %>% dplyr::select(chr = chr_num, tr_chr)
chr_labels    <- fai %>% mutate(chr = chr_num, mid = tr_chr + end_chr / 2)

# ─────────────────────────────────────────────────────────────
# 1) BOX PLOT OF PAIRWISE FRACTION IBD (with sensible fonts)
# ─────────────────────────────────────────────────────────────

fraction_ibd_plot <- fraction_ibd %>%
  filter(!is.na(fraction), !is.na(category)) %>%
  mutate(category = factor(category, levels = category_order))

g_box <- ggplot(fraction_ibd_plot,
                aes(x = category, y = fraction, fill = category)) +
  geom_boxplot(outlier.alpha = 0.2) +
  stat_summary(fun = mean, geom = "point", shape = 18,
               size = 1.8, color = "yellow") +
  theme_classic(base_size = 14) +  # moderate base font
  ylim(0, 0.25) +
  labs(x = NULL, y = "Pairwise fraction IBD") +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y  = element_text(size = 11),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )

ggsave(file.path(outdir, "ibd_fraction_boxplot.png"),
       g_box, width = 8, height = 6, dpi = 300)

# also save underlying table for convenience
readr::write_tsv(fraction_ibd_plot,
                 file.path(outdir, "fraction_ibd.tsv"))

# ─────────────────────────────────────────────────────────────
# 2) GENOME-WIDE FRACTION IBD PLOT
# ─────────────────────────────────────────────────────────────

combined_ibd_tr <- combined_ibd %>%
  filter(!is.na(fraction), !is.na(category),
         !is.na(chr), !is.na(win_start)) %>%
  mutate(chr = as.numeric(chr)) %>%
  left_join(transpose_chr, by = "chr") %>%
  mutate(pos_bp_ed = as.numeric(win_start) + tr_chr,
         category  = factor(category, levels = category_order))

signif_cutoff <- quantile(combined_ibd_tr$fraction,
                          probs = quantile_cutoff, na.rm = TRUE)

p_clean <- ggplot(combined_ibd_tr) +
  geom_rect(data = fai,
            aes(xmin = tr_chr, xmax = tr_chr + end_chr,
                ymin = -Inf, ymax = Inf, fill = fill_color),
            inherit.aes = FALSE, alpha = 0.5, color = NA) +
  scale_fill_identity() +
  geom_line(aes(x = pos_bp_ed, y = fraction, color = category),
            linewidth = 0.8) +
  facet_wrap(~ category, ncol = 1, scales = "free_x") +
  scale_x_continuous(breaks = chr_labels$mid, labels = chr_labels$chr) +
  scale_y_continuous(limits = c(0, 0.5)) +
  geom_hline(yintercept = signif_cutoff,
             linetype = "dashed", color = "red") +
  labs(x = "Chromosome", y = "IBD fraction", color = "Category") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_classic(base_size = 14) +
  theme(
    strip.text        = element_blank(),
    strip.background  = element_blank(),
    axis.text.x       = element_text(size = 11),
    axis.text.y       = element_text(size = 10),
    axis.title        = element_text(size = 16),
    legend.title      = element_text(size = 13),
    legend.text       = element_text(size = 11),
    legend.position   = "bottom"
  )

ggsave(file.path(outdir, "ibd_genomewide_fraction_cleaned.png"),
       p_clean, width = 12, height = 14, dpi = 300)

# ─────────────────────────────────────────────────────────────
# 3) CHROMOSOME PAINTING + DRUG RESISTANCE HIGHLIGHTING
# ─────────────────────────────────────────────────────────────

# (a) collapse high-IBD windows into contiguous segments
ibd_top <- combined_ibd_tr %>%
  filter(fraction > signif_cutoff) %>%
  transmute(category,
            chr   = as.integer(chr),
            start = as.numeric(win_start),
            end   = as.numeric(win_end))

collapse_intervals <- function(df) {
  if (nrow(df) == 0) return(df)
  dt <- as.data.table(df)
  setorder(dt, category, chr, start, end)
  dt[, lag_end := shift(end, type = "lag"), by = .(category, chr)]
  dt[, grp := cumsum(ifelse(is.na(lag_end) | start > lag_end, 1L, 0L)),
     by = .(category, chr)]
  out <- dt[, .(start = min(start), end = max(end)),
            by = .(category, chr, grp)]
  out[, grp := NULL]
  as_tibble(out)
}

ibd_segments <- collapse_intervals(ibd_top)

# (b) compute DR-gene overlapping segments from annotated file
#     (uses either PRODUCT text or GENE_ID)
res_gene_patterns <- c(
  # product-based keywords
  "Kelch13", "K13", "CRT", "chloroquine resistance transporter",
  "DHFR", "dihydrofolate reductase",
  "DHPS", "dihydropteroate synthetase", "PPPK-DHPS",
  "MDR1", "multidrug resistance protein 1",
  "UBP1", "ubiquitin carboxyl-terminal hydrolase 1",
  "coronin", "AP2-MU", "PX1",
  # gene IDs
  "PF3D7_0709000",  # CRT
  "PF3D7_1343700",  # K13
  "PF3D7_0417200",  # DHFR
  "PF3D7_0810800",  # DHPS
  "PF3D7_0104300",  # UBP1
  "PF3D7_0629500",  # AAT1
  "PF3D7_1224000",  # GCH1
  "PF3D7_0720700",  # PX1
  "PF3D7_1218300",  # AP2-MU
  "PF3D7_0523000"   # MDR1
)

dr_pattern_regex <- regex(paste(res_gene_patterns, collapse = "|"),
                          ignore_case = TRUE)

dr_segments <- annot_ibd %>%
  # make names lower case for easier handling
  rename_with(~tolower(.x)) %>%
  # standardise names after lowering
  rename(
    chr        = chr,
    win_start  = start,
    win_end    = end,
    gene_id    = gene_id,
    gene_name  = gene_name,
    gene_prod  = gene_product,
    category   = category
  ) %>%
  mutate(
    chr = as.integer(chr)
  ) %>%
  filter(!is.na(category), !is.na(chr), !is.na(win_start), !is.na(win_end)) %>%
  mutate(
    # detect DR genes via product OR ID/NAME
    text_for_match = paste(gene_id, gene_name, gene_prod),
    is_dr = str_detect(text_for_match, dr_pattern_regex)
  ) %>%
  filter(is_dr) %>%
  mutate(
    seg_start = win_start,
    seg_end   = win_end
  ) %>%
  transmute(category, chr, start = seg_start, end = seg_end) %>%
  distinct()

# (c) build chromosome “pill” scaffold
pill_x <- tibble(chr = sort(unique(fai$chr_num))) %>%
  mutate(x = chr)

chr_lens <- fai %>% transmute(chr = chr_num, len_bp = end_chr)
to_mb <- function(x) x / 1e6

pill_backbone <- tidyr::expand_grid(
  category = factor(category_order, levels = category_order),
  pill_x
) %>%
  left_join(chr_lens, by = "chr") %>%
  mutate(y0 = 0, y1 = to_mb(len_bp))

paint_segs <- ibd_segments %>%
  left_join(pill_x, by = "chr") %>%
  mutate(y0 = to_mb(start),
         y1 = to_mb(end))

dr_paint_segs <- dr_segments %>%
  left_join(pill_x, by = "chr") %>%
  mutate(y0 = to_mb(start),
         y1 = to_mb(end))

region_cols <- hue_pal()(length(category_order))
names(region_cols) <- category_order

paint_plot <- ggplot() +
  # backbone outline
  geom_segment(data = pill_backbone,
               aes(x = x, xend = x, y = y0, yend = y1),
               linewidth = 8.5, lineend = "round", color = "grey30") +
  # inner body
  geom_segment(data = pill_backbone,
               aes(x = x, xend = x, y = y0, yend = y1),
               linewidth = 6.0, lineend = "round", color = "white") +
  # high-IBD segments (per-category colour)
  geom_segment(data = paint_segs,
               aes(x = x, xend = x, y = y0, yend = y1, color = category),
               linewidth = 5, lineend = "butt", alpha = 0.9) +
  # DR-overlapping segments (bright red on top)
  geom_segment(data = dr_paint_segs,
               aes(x = x, xend = x, y = y0, yend = y1),
               linewidth = 5.4, lineend = "butt", color = "red", alpha = 0.95) +
  scale_color_manual(values = region_cols, guide = "none") +
  facet_wrap(~ factor(category, levels = category_order),
             ncol = 1, strip.position = "left") +
  # chromosome numbers above each pill
  geom_text(
    data = pill_backbone %>% group_by(category, chr, x) %>% slice_tail(),
    aes(x = x, y = y1 + 0.4, label = chr),
    size = 5
  ) +
  scale_y_continuous(
    limits = c(0, max(pill_backbone$y1, na.rm = TRUE) + 0.6),
    expand = expansion(mult = c(0, 0.005))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 14) +
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

message("Finished plotting IBD boxplot, genome-wide fraction, and chromosome painting (with DR highlights).")

