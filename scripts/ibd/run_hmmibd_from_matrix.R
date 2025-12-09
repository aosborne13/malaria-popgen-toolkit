#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

# ---- CLI ----
option_list <- list(
  make_option("--outdir", type="character"),
  make_option("--binary_matrix", type="character"),
  make_option("--metadata", type="character"),
  make_option("--category", type="character", default=NULL),
  make_option("--subgroup_col", type="character", default=NULL),
  make_option("--label_category", type="character", default="country"),
  make_option("--label_id", type="character", default="sample_id"),
  make_option("--label_fws", type="character", default="fws"),
  make_option("--fws_th", type="numeric", default=0.95),
  make_option("--maf", type="numeric", default=0.01),
  make_option("--na_char", type="character", default="N"),
  make_option("--threads", type="integer", default=4),
  make_option("--remove_chr", type="character", default=NULL),
  make_option("--hmmibd_bin", type="character", default="hmmIBD"),
  make_option("--skip_hmmibd", action="store_true", default=FALSE)
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)
setDTthreads(opt$threads)

# ---- Metadata ----
meta <- fread(opt$metadata)
stopifnot(opt$label_id %in% names(meta))
stopifnot(opt$label_category %in% names(meta))
stopifnot(opt$label_fws %in% names(meta))

if (!is.null(opt$subgroup_col)) {
  stopifnot(opt$subgroup_col %in% names(meta))
}

# ---- Matrix ----
message("Loading matrix...")
snp <- fread(opt$binary_matrix)
stopifnot(all(c("chr","pos","ref") %in% names(snp)[1:3]))

if (!is.null(opt$remove_chr)) {
  drop <- strsplit(opt$remove_chr, ",")[[1]]
  snp <- snp[!chr %in% drop]
}

geno_cols <- setdiff(names(snp), c("chr","pos","ref"))

# ---- Categories ----
cats <- if (is.null(opt$category)) unique(meta[[opt$label_category]]) else opt$category

for (cat in cats) {

  meta_cat <- meta |> filter(.data[[opt$label_category]] == cat,
                             .data[[opt$label_fws]] >= opt$fws_th)

  if (nrow(meta_cat) < 2) next

  subgroups <- if (is.null(opt$subgroup_col)) NA else unique(meta_cat[[opt$subgroup_col]])

  for (subg in subgroups) {

    tag <- if (is.na(subg)) cat else paste(cat, subg, sep="_")
    out_sub <- file.path(opt$outdir, tag)
    dir.create(out_sub, showWarnings=FALSE, recursive=TRUE)

    meta_use <- if (is.na(subg)) meta_cat else meta_cat |> filter(.data[[opt$subgroup_col]] == subg)

    samples <- intersect(meta_use[[opt$label_id]], geno_cols)
    if (length(samples) < 2) next

    g <- snp[, c("chr","pos",samples), drop=FALSE]
    g[g == opt$na_char] <- NA
    g[g == "N"] <- NA
    g[g == "."] <- NA
    g <- apply(g[,-(1:2)], 2, as.numeric)

    g[g == 0.5] <- 1
    g[is.na(g)] <- -1

    # MAF filter
    p <- rowMeans(g == 1, na.rm=TRUE)
    maf <- pmin(p, 1-p)
    keep <- which(maf >= opt$maf)

    if (length(keep) == 0) next
    g2 <- g[keep,,drop=FALSE]
    outmat <- cbind(snp$chr[keep], snp$pos[keep], g2)
    fout <- file.path(out_sub, "hmmIBD_input.txt")
    fwrite(outmat, fout, sep="\t")

    if (!opt$skip_hmmibd) {
      cmd <- sprintf("cd %s && %s -i hmmIBD_input.txt -o hmmIBD_out",
                     shQuote(out_sub), shQuote(opt$hmmibd_bin))
      system(cmd)
    }
  }
}

message("Done.")


