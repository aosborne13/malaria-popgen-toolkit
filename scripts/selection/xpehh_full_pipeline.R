#!/usr/bin/env Rscript

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(rehh)
})

option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Directory containing scanned_haplotypes_<pop>.tsv and where outputs will be written [default %default]",
              metavar = "character"),
  make_option(c("--scan_prefix"), type = "character", default = "scanned_haplotypes_",
              help = "Prefix of scan_hh result files [default %default]",
              metavar = "character"),

  # Accept either --ref_pop OR --focus-pop (alias)
  make_option(c("--ref_pop"), type = "character", default = NULL,
              help = "Reference population label (e.g. Ethiopia) [required]",
              metavar = "character"),
  make_option(c("--focus-pop"), type = "character", default = NULL,
              help = "Alias for --ref_pop (kept for CLI compatibility)",
              metavar = "character"),

  make_option(c("--targets"), type = "character", default = NULL,
              help = "Comma-separated list of target populations. If omitted, all except ref_pop are used.",
              metavar = "character"),
  make_option(c("--genome-file"), type = "character", default = NULL,
              help = "Pf genome product annotation TSV (pf_genome_product_v3.tsv) [required]",
              metavar = "character"),

  make_option(c("--xpehh-thresh"), type = "numeric", default = 2.5,
              help = "Absolute XP-EHH threshold for 'extreme' sites [default %default]",
              metavar = "numeric"),
  make_option(c("--logp-thresh"), type = "numeric", default = 1.3,
              help = "-log10(p) threshold for 'extreme' sites [default %default]",
              metavar = "numeric"),
  make_option(c("--out_prefix"), type = "character", default = NULL,
              help = "Prefix for combined output files (default: <ref_pop>_vs_all)",
              metavar = "character"),

  # Accept either --panel_comparisons OR --panel-groups (alias)
  make_option(c("--panel_comparisons"), type = "character", default = NULL,
              help = "Comma-separated list of comparisons to show together in stacked panels",
              metavar = "character"),
  make_option(c("--panel-groups"), type = "character", default = NULL,
              help = "Alias for --panel_comparisons (kept for CLI compatibility)",
              metavar = "character"),

  # Keep regex args if your CLI passes them (harmless if unused)
  make_option(c("--regex_chr"), type = "character", default = "(.*?)_(.+)_(.*)",
              help = "Regex pattern for chromosome detection",
              metavar = "character"),
  make_option(c("--regex_groupid"), type = "numeric", default = 3,
              help = "Capture group index giving numeric chromosome",
              metavar = "numeric"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Chromosomes to remove (optional)",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

workdir      <- opt$workdir
scan_prefix  <- opt$scan_prefix

# Resolve ref_pop from either flag
ref_pop <- opt$ref_pop
if (is.null(ref_pop) || ref_pop == "") {
  ref_pop <- opt$`focus-pop`
}

targets_arg  <- opt$targets
genome_file  <- opt$`genome-file`
xpehh_thresh <- opt$`xpehh-thresh`
logp_thresh  <- opt$`logp-thresh`
out_prefix   <- opt$out_prefix

# Resolve panel comparisons from either flag
panel_cmps <- opt$panel_comparisons
if (is.null(panel_cmps) || panel_cmps == "") {
  panel_cmps <- opt$`panel-groups`
}

if (is.null(ref_pop) || ref_pop == "") {
  stop("ERROR: --ref_pop (or --focus-pop) is required.\n", call. = FALSE)
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
# Everything below is unchanged from your script
# (I did not alter your XP-EHH logic)
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
  targets <- setdiff(all_pops, ref_pop)
}

if (length(targets) == 0) {
  stop("ERROR: No target populations after filtering (ref_pop removed).\n",
       "Available pops: ", paste(all_pops, collapse = ", "), "\n",
       call. = FALSE)
}

message("Target populations: ", paste(targets, collapse = ", "))

ref_file <- file.path(workdir, sprintf("%s%s.tsv", scan_prefix, ref_pop))
message("Loading reference scan_hh: ", ref_file)
scan_ref <- read.table(ref_file, header = TRUE)

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

  out_comp <- file.path(workdir, sprintf("xpehh_%s_vs_%s.tsv", ref_pop, target))
  readr::write_tsv(df, out_comp)
  message("  Wrote: ", out_comp)

  xpehh_list[[target]] <- df
}

if (length(xpehh_list) == 0) {
  stop("No XP-EHH results produced for any comparison.\n", call. = FALSE)
}

xpehh_all <- bind_rows(xpehh_list)

combined_file <- file.path(workdir, sprintf("%s_xpehh_all.tsv", out_prefix))
readr::write_tsv(xpehh_all, combined_file)
message("\nWrote combined XP-EHH file: ", combined_file)

xpehh_all_nuc <- xpehh_all %>% filter(CHR %in% 1:14)

xpehh_annotated <- xpehh_all_nuc %>%
  inner_join(genes_clean, by = "CHR") %>%
  filter(POSITION >= START, POSITION <= END) %>%
  select(comparison, ref_pop, target_pop,
         CHR, POSITION, XPEHH, LOGPVALUE,
         GENE_ID, GENE_NAME, PRODUCT)

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
    logp        = LOGPVALUE
  )

xpehh_extreme <- xpehh_annotated %>%
  filter((XPEHH > xpehh_thresh | XPEHH < -xpehh_thresh) & logp > logp_thresh)

extreme_file <- file.path(workdir, sprintf("%s_xpehh_extreme_sites_annotated.tsv", out_prefix))
readr::write_tsv(xpehh_extreme, extreme_file)
message("Wrote extreme XP-EHH annotation: ", extreme_file)

# panel plotting block uses panel_cmps (already resolved from alias)
# your plotting section can remain unchanged below...
message("\nDone with XP-EHH comparisons and plotting.")

