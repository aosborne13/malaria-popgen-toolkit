#!/usr/bin/env Rscript

# Neighbour-joining tree from a binary SNP matrix.
#
# Distances are Manhattan (amap::Dist), the same metric pca-plot uses, so a tree
# and a PCoA built from the same matrix are describing the same distances. The
# tree itself is ape::nj(), with optional bootstrap support from ape::boot.phylo().
#
# Genome-wide by default; --region builds a tree from one window, which is what
# you want when asking whether isolates sharing a locus also share its haplotype.
#
# Outputs Newick plus a tip table carrying the metadata, so the tree can be
# plotted and annotated downstream without rejoining anything by hand.

options(scipen = 999)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(amap)
  library(ape)
})

option_list <- list(
  make_option(c("-d", "--workdir"), type = "character", default = ".",
              help = "Output directory [default %default]", metavar = "character"),
  make_option(c("-b", "--matrix_binary"), type = "character", default = NULL,
              help = "Binary SNP matrix: chr, pos, ref, <samples> [required]",
              metavar = "character"),
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "Metadata TSV [required]", metavar = "character"),
  make_option(c("--label_id"), type = "character", default = "sample_id",
              help = "Sample ID column [default %default]", metavar = "character"),
  make_option(c("--label_fws"), type = "character", default = NULL,
              help = "Optional Fws column; samples below --fws_th are dropped",
              metavar = "character"),
  make_option(c("--fws_th"), type = "numeric", default = 0.95,
              help = "Fws threshold, used only with --label_fws [default %default]",
              metavar = "numeric"),
  make_option(c("--region"), type = "character", default = NULL,
              help = "Restrict to one window, as chr:start-end", metavar = "character"),
  make_option(c("--remove_chr"), type = "character", default = NULL,
              help = "Comma-separated chromosomes to drop", metavar = "character"),
  make_option(c("--mixed_call"), type = "character", default = "alt",
              help = "How to treat mixed (0.5) calls: alt, missing, or keep [default %default]",
              metavar = "character"),
  make_option(c("--max_sample_missing"), type = "numeric", default = 0.2,
              help = "Drop samples missing more than this fraction [default %default]",
              metavar = "numeric"),
  make_option(c("--max_locus_missing"), type = "numeric", default = 0.2,
              help = "Drop loci missing in more than this fraction [default %default]",
              metavar = "numeric"),
  make_option(c("--min_shared"), type = "integer", default = 50,
              help = "Warn if any pair shares fewer than this many called loci [default %default]",
              metavar = "integer"),
  make_option(c("--max_loci"), type = "integer", default = NA,
              help = "Randomly thin to at most this many loci (speeds up bootstrap)",
              metavar = "integer"),
  make_option(c("--n_boot"), type = "integer", default = 0,
              help = "Bootstrap replicates; 0 skips bootstrapping [default %default]",
              metavar = "integer"),
  make_option(c("--na_char"), type = "character", default = "N",
              help = "Missing-genotype character [default %default]", metavar = "character"),
  make_option(c("--seed"), type = "integer", default = 1,
              help = "Random seed [default %default]", metavar = "integer"),
  make_option(c("--prefix"), type = "character", default = "nj_tree",
              help = "Output file prefix [default %default]", metavar = "character"),
  make_option(c("--plot"), action = "store_true", default = FALSE,
              help = "Also write a quick unrooted PDF of the tree")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$matrix_binary)) stop("--matrix_binary is required.")
if (is.null(opt$metadata))      stop("--metadata is required.")
if (!opt$mixed_call %in% c("alt", "missing", "keep")) {
  stop("--mixed_call must be alt, missing or keep.")
}

set.seed(opt$seed)
dir.create(opt$workdir, showWarnings = FALSE, recursive = TRUE)
split_csv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])
out <- function(...) file.path(opt$workdir, paste0(opt$prefix, ...))

# ---- metadata --------------------------------------------------------
meta <- read.csv(opt$metadata, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (!opt$label_id %in% colnames(meta)) {
  stop("No column '", opt$label_id, "' in metadata.")
}
meta[[opt$label_id]] <- as.character(meta[[opt$label_id]])
meta <- meta[!duplicated(meta[[opt$label_id]]), ]

if (!is.null(opt$label_fws)) {
  if (!opt$label_fws %in% colnames(meta)) stop("No column '", opt$label_fws, "' in metadata.")
  n0 <- nrow(meta)
  meta <- meta[!is.na(meta[[opt$label_fws]]) & meta[[opt$label_fws]] >= opt$fws_th, ]
  message("Fws >= ", opt$fws_th, ": ", nrow(meta), " of ", n0, " samples retained")
}

# ---- matrix ----------------------------------------------------------
message("Reading matrix...")
snp <- fread(opt$matrix_binary, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("chr", "pos") %in% colnames(snp))) {
  stop("Matrix must carry 'chr' and 'pos' columns.")
}

if (!is.null(opt$remove_chr)) {
  snp <- snp[!snp$chr %in% split_csv(opt$remove_chr), , drop = FALSE]
}

if (!is.null(opt$region)) {
  rr <- regmatches(opt$region, regexec("^(.+):([0-9]+)-([0-9]+)$", opt$region))[[1]]
  if (length(rr) != 4) stop("--region must look like chr:start-end")
  rchr <- rr[2]; rs <- as.numeric(rr[3]); re <- as.numeric(rr[4])
  snp <- snp[snp$chr == rchr & snp$pos >= rs & snp$pos <= re, , drop = FALSE]
  message("Region ", opt$region, ": ", nrow(snp), " loci")
  if (!nrow(snp)) stop("No loci in that region.")
}

samples <- intersect(meta[[opt$label_id]], colnames(snp))
if (length(samples) < 3) {
  stop("Need at least 3 samples shared between matrix and metadata; found ", length(samples))
}
geno <- snp[, samples, drop = FALSE]
pos_chr <- snp$chr; pos_bp <- snp$pos
rm(snp)

geno[geno == opt$na_char] <- NA
geno[geno == "."] <- NA
geno[] <- lapply(geno, function(x) suppressWarnings(as.numeric(x)))
G <- as.matrix(geno); rm(geno)

if (opt$mixed_call == "alt") {
  G[G > 0 & G < 1] <- 1
} else if (opt$mixed_call == "missing") {
  G[G > 0 & G < 1] <- NA
}

# ---- filtering -------------------------------------------------------
locus_miss <- rowMeans(is.na(G))
keep_l <- locus_miss <= opt$max_locus_missing
message("Loci: ", sum(keep_l), " of ", length(keep_l),
        " at <= ", opt$max_locus_missing, " missing")
G <- G[keep_l, , drop = FALSE]
pos_chr <- pos_chr[keep_l]; pos_bp <- pos_bp[keep_l]

samp_miss <- colMeans(is.na(G))
keep_s <- samp_miss <= opt$max_sample_missing
if (sum(keep_s) < 3) stop("Fewer than 3 samples survive --max_sample_missing.")
message("Samples: ", sum(keep_s), " of ", length(keep_s),
        " at <= ", opt$max_sample_missing, " missing")
G <- G[, keep_s, drop = FALSE]

# drop invariant loci - they carry no distance information
nvar <- apply(G, 1, function(r) length(unique(r[!is.na(r)])))
G <- G[nvar > 1, , drop = FALSE]
pos_chr <- pos_chr[nvar > 1]; pos_bp <- pos_bp[nvar > 1]
message("Variant loci retained: ", nrow(G))
if (nrow(G) < 10) stop("Fewer than 10 variant loci remain; a tree would not be meaningful.")

if (!is.na(opt$max_loci) && nrow(G) > opt$max_loci) {
  idx <- sort(sample.int(nrow(G), opt$max_loci))
  G <- G[idx, , drop = FALSE]; pos_chr <- pos_chr[idx]; pos_bp <- pos_bp[idx]
  message("Thinned to ", nrow(G), " loci (--max_loci)")
}

X <- t(G)   # samples x loci, which is what Dist() and boot.phylo() both want
n_mrk <- ncol(X)
message(nrow(X), " samples x ", n_mrk, " markers")

# ---- pairwise sharing diagnostic ------------------------------------
# A pair sharing very few called loci gets a distance estimated from almost
# nothing, and amap's scaling will inflate it rather than flag it.
fin <- !is.na(X)
shared <- fin %*% t(fin)
diag(shared) <- NA
min_shared <- min(shared, na.rm = TRUE)
message("Fewest jointly called loci in any pair: ", min_shared)
if (min_shared < opt$min_shared) {
  n_bad <- sum(shared < opt$min_shared, na.rm = TRUE) / 2
  warning(n_bad, " pair(s) share fewer than ", opt$min_shared,
          " called loci; their distances rest on very little data.",
          call. = FALSE)
}

# ---- distances and tree ---------------------------------------------
message("Computing Manhattan distances...")
D <- amap::Dist(X, method = "manhattan")

dm <- as.matrix(D)
fwrite(data.table(sample_id = rownames(dm), as.data.frame(dm)),
       out(".distance.tsv"), sep = "\t")

message("Building neighbour-joining tree...")
tree <- ape::nj(D)

if (opt$n_boot > 0) {
  message("Bootstrapping (", opt$n_boot, " replicates)...")
  FUN <- function(xx) ape::nj(amap::Dist(xx, method = "manhattan"))
  bp <- ape::boot.phylo(tree, X, FUN, B = opt$n_boot, rooted = FALSE,
                        trees = FALSE, quiet = TRUE)
  tree$node.label <- round(100 * bp / opt$n_boot, 1)
  message("Median bootstrap support: ",
          round(median(as.numeric(tree$node.label), na.rm = TRUE), 1), "%")
}

ape::write.tree(tree, file = out(".nwk"))

# ---- tip table -------------------------------------------------------
tips <- data.frame(sample_id = tree$tip.label, stringsAsFactors = FALSE)
names(tips)[1] <- opt$label_id
tips <- merge(tips, meta, by = opt$label_id, all.x = TRUE, sort = FALSE)
tips <- tips[match(tree$tip.label, tips[[opt$label_id]]), ]
fwrite(tips, out("_tips.tsv"), sep = "\t")

if (opt$plot) {
  pdf(out(".pdf"), width = 10, height = 10)
  plot(tree, type = "unrooted", show.tip.label = (ape::Ntip(tree) <= 120),
       cex = 0.5, no.margin = TRUE, lab4ut = "axial")
  if (opt$n_boot > 0) ape::nodelabels(tree$node.label, frame = "none", cex = 0.4, adj = c(1.1, -0.3))
  invisible(dev.off())
}

message("\nWritten to ", opt$workdir, ":")
message("  ", basename(out(".nwk")), "            Newick tree",
        if (opt$n_boot > 0) " with bootstrap node labels" else "")
message("  ", basename(out("_tips.tsv")), "       tip order with metadata joined")
message("  ", basename(out(".distance.tsv")), "   pairwise Manhattan distance matrix")
if (opt$plot) message("  ", basename(out(".pdf")), "            quick unrooted plot")
message("\nMarkers used: ", n_mrk,
        if (!is.null(opt$region)) paste0(" (region ", opt$region, ")") else " (genome-wide)")
