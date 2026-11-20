#!/usr/bin/env Rscript

# Spatially-explicit ancestry estimation with TESS3.
#
# Takes the binary SNP matrix and a metadata file carrying longitude/latitude
# per sample, runs tess3() over one K or a range of K, and writes the ancestry
# coefficients together with a cross-entropy curve, an admixture barplot and a
# map of interpolated ancestry.
#
# Running a range (--K 1-10) gives the cross-entropy and RMSE per K needed to
# justify a chosen K; running a single value (--K 6) produces the final plots.
#
# The transposed genotype matrix is expensive to build, so it is cached as
# <workdir>/<prefix>.matbinT.rds and reused on later runs unless --force-rebuild.

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Output directory [default %default]", metavar = "character"),
  make_option(c("-b", "--matrix_binary"), type = "character", default = NULL,
              help = "Binary SNP matrix: chr, pos, ref, <samples>", metavar = "character"),
  make_option(c("--matrix_rds"), type = "character", default = NULL,
              help = "Pre-built samples x SNPs RDS; skips matrix parsing",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Metadata TSV [required]", metavar = "character"),
  make_option(c("--label_id"), type = "character", default = "sample_id",
              help = "Sample ID column [default %default]", metavar = "character"),
  make_option(c("--lon_col"), type = "character", default = "longitude",
              help = "Longitude column [default %default]", metavar = "character"),
  make_option(c("--lat_col"), type = "character", default = "latitude",
              help = "Latitude column [default %default]", metavar = "character"),
  make_option(c("--group_col"), type = "character", default = "country",
              help = "Column used to order and facet the barplot [default %default]",
              metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = NULL,
              help = "Optional Fws column; samples below --fws_th are dropped",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold, used only with --label_fws [default %default]",
              metavar = "numeric"),
  make_option(c("-K", "--K"), type = "character", default = "1-10",
              help = "Single K (e.g. 6) or inclusive range (e.g. 1-10) [default %default]",
              metavar = "character"),
  make_option(c("--rep"), type = "integer", default = 10,
              help = "Repetitions per K [default %default]", metavar = "integer"),
  make_option(c("--method"), type = "character", default = "projected.ls",
              help = "tess3 method: projected.ls or qp [default %default]",
              metavar = "character"),
  make_option(c("--lambda"), type = "numeric", default = 1,
              help = "Spatial regularisation parameter [default %default]",
              metavar = "numeric"),
  make_option(c("--max_iteration"), type = "integer", default = 200,
              help = "Maximum iterations [default %default]", metavar = "integer"),
  make_option(c("--ploidy"), type = "integer", default = 1,
              help = "Ploidy [default %default]", metavar = "integer"),
  make_option(c("--na_char"), type = "character", default = "N",
              help = "Missing-genotype character [default %default]", metavar = "character"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Comma-separated chromosomes to drop", metavar = "character"),
  make_option(c("--prefix"), type = "character", default = "tess3r",
              help = "Output file prefix [default %default]", metavar = "character"),
  make_option(c("--force_rebuild"), action = "store_true", default = FALSE,
              help = "Rebuild the cached genotype matrix even if present"),
  make_option(c("--threads"), type = "integer", default = 4,
              help = "openMP cores [default %default]", metavar = "integer")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$metadata)) stop("--metadata is required.")
if (is.null(opt$matrix_binary) && is.null(opt$matrix_rds)) {
  stop("Supply --matrix_binary or --matrix_rds.")
}
if (!requireNamespace("tess3r", quietly = TRUE)) {
  stop("Package 'tess3r' not installed. Run scripts/setup_env.sh, or:\n",
       "  Rscript -e 'devtools::install_github(\"bcm-uga/TESS3_encho_sen\")'")
}
suppressPackageStartupMessages(library(tess3r))

dir.create(opt$workdir, showWarnings = FALSE, recursive = TRUE)
split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

Ks <- if (grepl("-", opt$K)) {
  rng <- as.integer(strsplit(opt$K, "-")[[1]]); seq(rng[1], rng[2])
} else as.integer(opt$K)
if (any(is.na(Ks)) || any(Ks < 1)) stop("--K must be a positive integer or range like 1-10.")

# ---- metadata --------------------------------------------------------
meta <- read.csv(opt$metadata, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
need <- c(opt$label_id, opt$lon_col, opt$lat_col)
miss <- setdiff(need, colnames(meta))
if (length(miss)) stop("Metadata missing column(s): ", paste(miss, collapse = ", "))

if (!is.null(opt$label_fws)) {
  if (!opt$label_fws %in% colnames(meta)) stop("No column '", opt$label_fws, "' in metadata.")
  meta <- meta[!is.na(meta[[opt$label_fws]]) & meta[[opt$label_fws]] >= opt$fws_th, ]
  message("Fws >= ", opt$fws_th, ": ", nrow(meta), " samples")
}
meta <- meta[!is.na(meta[[opt$lon_col]]) & !is.na(meta[[opt$lat_col]]), ]
if (!nrow(meta)) stop("No samples retain coordinates; TESS3 requires longitude and latitude.")

# ---- genotypes -------------------------------------------------------
cache <- file.path(opt$workdir, paste0(opt$prefix, ".matbinT.rds"))

if (!is.null(opt$matrix_rds)) {
  message("Loading genotype matrix from ", opt$matrix_rds)
  X <- readRDS(opt$matrix_rds)
} else if (file.exists(cache) && !opt$force_rebuild) {
  message("Reusing cached genotype matrix: ", cache)
  X <- readRDS(cache)
} else {
  message("Reading matrix...")
  snp <- fread(opt$matrix_binary, sep = "\t", header = TRUE, data.table = FALSE)
  if (!is.null(opt$remove_chr)) snp <- snp[!snp$chr %in% split_csv(opt$remove_chr), , drop = FALSE]
  samples <- intersect(meta[[opt$label_id]], colnames(snp))
  if (!length(samples)) stop("No metadata samples found in the matrix header.")
  geno <- snp[, samples, drop = FALSE]; rm(snp)
  geno[geno == opt$na_char] <- NA
  geno[geno == "."] <- NA
  geno[] <- lapply(geno, function(x) suppressWarnings(as.numeric(x)))
  G <- as.matrix(geno); rm(geno)
  G[G > 0] <- 1                        # mixed calls -> alternate, matching the main pipeline
  X <- t(G); rm(G)                     # tess3 wants samples x SNPs
  saveRDS(X, cache)
  message("Cached genotype matrix at ", cache)
}

meta <- meta[match(rownames(X), meta[[opt$label_id]]), ]
meta <- meta[!is.na(meta[[opt$label_id]]), ]
X <- X[meta[[opt$label_id]], , drop = FALSE]
coords <- as.matrix(meta[, c(opt$lon_col, opt$lat_col)])
storage.mode(coords) <- "numeric"
message(nrow(X), " samples x ", ncol(X), " SNPs")

# ---- run TESS3 -------------------------------------------------------
summary_rows <- list()
for (K in Ks) {
  message("\nRunning TESS3, K = ", K, " (rep = ", opt$rep, ")")
  obj <- tess3r::tess3(X = X, coord = coords, K = K,
                       method = opt$method, lambda = opt$lambda,
                       ploidy = opt$ploidy, rep = opt$rep,
                       max.iteration = opt$max_iteration,
                       openMP.core.num = opt$threads, verbose = FALSE)
  saveRDS(obj, file.path(opt$workdir, sprintf("%s_rep%d_K%d.rds", opt$prefix, opt$rep, K)))
  res <- tess3r::Gettess3res(obj, K = K)
  summary_rows[[as.character(K)]] <- data.frame(K = K,
                                                crossentropy = res$crossentropy,
                                                rmse = res$rmse)
  message("  cross-entropy ", signif(res$crossentropy, 5), " | RMSE ", signif(res$rmse, 5))
}
cv <- bind_rows(summary_rows)
fwrite(cv, file.path(opt$workdir, paste0(opt$prefix, ".cross_entropy.tsv")), sep = "\t")

if (length(Ks) > 1) {
  g <- ggplot(cv, aes(x = K, y = crossentropy)) +
    geom_line(colour = "grey50") + geom_point(size = 2) +
    scale_x_continuous(breaks = Ks) +
    labs(x = "Number of ancestral populations (K)", y = "Cross-entropy") +
    theme_classic(base_size = 11)
  ggsave(file.path(opt$workdir, paste0(opt$prefix, ".cross_entropy.pdf")), g,
         width = 6, height = 4)
  message("\nCross-entropy curve written. Lowest at K = ",
          cv$K[which.min(cv$crossentropy)],
          ". Choose K from this curve, then rerun with a single --K for the plots.")
  quit(save = "no", status = 0)
}

# ---- single K: ancestry outputs -------------------------------------
K <- Ks[1]
obj <- readRDS(file.path(opt$workdir, sprintf("%s_rep%d_K%d.rds", opt$prefix, opt$rep, K)))
qm <- tess3r::qmatrix(obj, K = K)

qdf <- as.data.frame(qm)
colnames(qdf) <- paste0("K", seq_len(ncol(qdf)))
qdf[[opt$label_id]] <- meta[[opt$label_id]]
if (opt$group_col %in% colnames(meta)) qdf[[opt$group_col]] <- meta[[opt$group_col]]
fwrite(qdf, file.path(opt$workdir, sprintf("%s_K%d.qmatrix.tsv", opt$prefix, K)), sep = "\t")

PALETTE <- c("#E69F00", "#0072B2", "#009E73", "#CC79A7", "#56B4E9",
             "#D55E00", "#8C6D31", "#5E7490", "#9A6FB0", "#7F7F7F")
cols <- setNames(PALETTE[seq_len(K)], paste0("K", seq_len(K)))

long <- qdf %>%
  pivot_longer(all_of(names(cols)), names_to = "K_component", values_to = "Q") %>%
  group_by(.data[[opt$label_id]]) %>%
  mutate(max_K = K_component[which.max(Q)]) %>%
  ungroup()

if (opt$group_col %in% colnames(qdf)) {
  ord <- long %>%
    distinct(.data[[opt$label_id]], .data[[opt$group_col]], max_K) %>%
    left_join(long %>% group_by(.data[[opt$label_id]]) %>%
                summarise(top = max(Q), .groups = "drop"), by = opt$label_id) %>%
    arrange(.data[[opt$group_col]], max_K, desc(top))
  long[[opt$label_id]] <- factor(long[[opt$label_id]], levels = ord[[opt$label_id]])

  adm <- ggplot(long, aes(x = .data[[opt$label_id]], y = Q, fill = K_component)) +
    geom_col(width = 1) +
    facet_grid(reformulate(opt$group_col), scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(values = cols, name = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Ancestry") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank(), legend.position = "bottom",
          strip.text.x = element_text(angle = 90, hjust = 1, size = 7))
  ggsave(file.path(opt$workdir, sprintf("%s_K%d.admixture.pdf", opt$prefix, K)), adm,
         width = 12, height = 4.5)

  grp <- long %>% group_by(.data[[opt$group_col]], K_component) %>%
    summarise(mean_Q = mean(Q), .groups = "drop")
  fwrite(grp, file.path(opt$workdir, sprintf("%s_K%d.ancestry_by_%s.tsv", opt$prefix, K, opt$group_col)),
         sep = "\t")
}

# ---- interpolated ancestry map --------------------------------------
if (requireNamespace("rworldmap", quietly = TRUE)) {
  ok <- try({
    pal <- tess3r::CreatePalette(unname(cols), 9)
    poly <- rworldmap::getMap(resolution = "low")
    pl <- tess3r::ggtess3Q(qm, coords, map.polygon = poly, col.palette = pal) +
      geom_point(data = as.data.frame(coords),
                 aes(x = .data[[opt$lon_col]], y = .data[[opt$lat_col]]),
                 colour = "black", size = 0.6) +
      coord_equal() + labs(x = "Longitude", y = "Latitude") +
      theme_classic(base_size = 11)
    ggsave(file.path(opt$workdir, sprintf("%s_K%d.map.pdf", opt$prefix, K)), pl,
           width = 7, height = 6)
  }, silent = TRUE)
  if (inherits(ok, "try-error")) {
    message("Map step failed (this is cosmetic; Q-matrix and barplot are unaffected):\n  ",
            conditionMessage(attr(ok, "condition")))
  }
} else {
  message("Package 'rworldmap' not installed; skipping the map. ",
          "Add r-rworldmap to env/environment.yml to enable it.")
}

message("\nWritten to ", opt$workdir, ":")
message("  ", sprintf("%s_K%d.qmatrix.tsv", opt$prefix, K), "       ancestry coefficients per sample")
message("  ", sprintf("%s_K%d.admixture.pdf", opt$prefix, K), "     barplot ordered by ", opt$group_col)
message("  ", sprintf("%s_K%d.map.pdf", opt$prefix, K), "           interpolated ancestry map")
message("  ", paste0(opt$prefix, ".cross_entropy.tsv"), "        cross-entropy and RMSE")
