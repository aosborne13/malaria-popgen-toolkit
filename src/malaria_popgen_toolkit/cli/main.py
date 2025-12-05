# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.
Commands:
  - missense-drugres-af   (AF of missense drug-res variants with DP filter; flexible grouping)
  - hapmap-africa         (haplotype map for Africa)
  - hapmap-samerica       (haplotype map for South America)
  - hapmap-seasia         (haplotype map for Southeast Asia)
  - fws-dotplot           (Fws jittered dot plots by metadata groups)
  - pca                   (Distance-based PCA/PCoA from binary matrix or VCF)
  - hmmibd-matrix         (Prepare hmmIBD input from binary matrix and run hmmIBD via R)
"""

import argparse
import sys
import shutil

from malaria_popgen_toolkit.commands import (
    missense_drugres_af,
    haplotype_map_region,
    fws_dotplot,
    pca_plot,
    hmmibg_matrix,
)


def require_tool(name: str):
    """Exit with an error message if a required tool is missing."""
    if shutil.which(name) is None:
        sys.exit(f"ERROR: '{name}' not found in PATH. Please install it and try again.")


def main():
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # -------------------------------------------------------------------------
    # missense-drugres-af
    # -------------------------------------------------------------------------
    p1 = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes (DP filter)"
    )
    p1.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p1.add_argument("--ref", required=True, help="Reference FASTA file")
    p1.add_argument("--gff3", required=True, help="GFF3 annotation file")
    p1.add_argument("--metadata", required=True, help="TSV with 'sample_id' and grouping column(s)")
    p1.add_argument("--outdir", default="output_missense", help="Output directory")
    p1.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
    p1.add_argument(
        "--group-by",
        default="country",
        help="Metadata column to group by (e.g., country, region, site, year)"
    )

    # -------------------------------------------------------------------------
    # hapmap-africa / samerica / seasia
    # -------------------------------------------------------------------------
    for cmd, region in [
        ("hapmap-africa", "africa"),
        ("hapmap-samerica", "south_america"),
        ("hapmap-seasia", "southeast_asia"),
    ]:
        px = sub.add_parser(
            cmd,
            help=f"Plot haplotype pies for {region.replace('_', ' ').title()} from VCF + metadata"
        )
        px.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
        px.add_argument("--metadata", required=True, help="TSV with samples and country column")
        px.add_argument("--outdir", required=True, help="Output directory for map + CSVs")
        px.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP for hap calls")
        px.add_argument("--sample-col", default="sample_id", help="Sample column name in metadata")
        px.add_argument("--country-col", default="country", help="Country column in metadata")

    # -------------------------------------------------------------------------
    # fws-dotplot
    # -------------------------------------------------------------------------
    p5 = sub.add_parser(
        "fws-dotplot",
        help="Create jittered dot plots of Fws by metadata grouping columns"
    )
    p5.add_argument("--metadata", required=True, help="TSV with at least 'fws' and grouping columns")
    p5.add_argument("--outdir", default="fws_plots", help="Output directory for PDFs")
    p5.add_argument(
        "--group-by",
        action="append",
        help="Metadata column to group by (repeat for multiple). "
             "Default tries: region, country, year."
    )
    p5.add_argument("--width", type=float, default=10.0, help="Plot width (inches)")
    p5.add_argument("--height", type=float, default=6.0, help="Plot height (inches)")

    # -------------------------------------------------------------------------
    # PCA (distance-based PCoA)
    # -------------------------------------------------------------------------
    p6 = sub.add_parser(
        "pca",
        help="Distance-based PCA/PCoA from genotype matrix or VCF (Manhattan SNP differences)"
    )

    p6.add_argument("--matrix", help="Binary genotype matrix (.tsv) with samples as columns")
    p6.add_argument("--vcf", help="VCF file to convert into genotype matrix")
    p6.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p6.add_argument("--outdir", default="pca_plots", help="Output directory for PCA plots")
    p6.add_argument(
        "--sample-col",
        default="sample_id",
        help="Metadata sample column (default: 'sample_id')"
    )
    p6.add_argument(
        "--group-by",
        nargs="+",
        help="Metadata columns to color by. Default tries region, country, year."
    )
    p6.add_argument(
        "--pcs",
        nargs="+",
        help="PC pairs to plot, e.g. --pcs 1,2 1,3 (default: 1,2 and 1,3)"
    )
    p6.add_argument(
        "--max-sample-missing",
        type=float,
        help="Drop samples with >X fraction missing (e.g., 0.3)"
    )

    # -------------------------------------------------------------------------
    # hmmibd-matrix (IBD / hmmIBD from binary matrix)
    # -------------------------------------------------------------------------
    p7 = sub.add_parser(
        "hmmibd-matrix",
        help="Prepare hmmIBD input from a binary SNP matrix and run hmmIBD via R"
    )

    p7.add_argument(
        "--matrix",
        required=True,
        help="Binary SNP genotype matrix (0/0.5/1/N) with columns: chr,pos,ref,<samples>",
    )
    p7.add_argument(
        "--metadata",
        required=True,
        help="Metadata TSV with sample IDs, category column, and Fws column",
    )
    p7.add_argument(
        "--workdir",
        required=True,
        help="Working/output directory for hmmIBD inputs/outputs",
    )

    p7.add_argument(
        "--category-col",
        default="country",
        help="Metadata column defining category/group (e.g. country, region, year)",
    )
    p7.add_argument(
        "--category",
        default=None,
        help="Category name (e.g. 'Ethiopia'). If omitted, the R script can loop over all categories.",
    )
    p7.add_argument(
        "--sample-col",
        default="sample_id",
        help="Sample ID column in metadata (must match matrix column names)",
    )
    p7.add_argument(
        "--fws-col",
        default="fws",
        help="Metadata column containing Fws values",
    )
    p7.add_argument(
        "--fws-th",
        type=float,
        default=0.95,
        help="Fws threshold (keep samples with Fws ≥ this value)",
    )
    p7.add_argument(
        "--maf",
        type=float,
        default=0.01,
        help="Minor allele frequency threshold for SNPs before hmmIBD",
    )
    p7.add_argument(
        "--na-char",
        default="NA",
        help="Character used for missing genotypes in matrix (also '.' and 'N' treated as missing)",
    )
    p7.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for data.table inside the R script",
    )
    p7.add_argument(
        "--exclude-chr",
        default="Pf3D7_API_v3,Pf3D7_MIT_v3",
        help="Comma-separated list of chromosomes to drop (e.g. Pf3D7_API_v3,Pf3D7_MIT_v3)",
    )
    p7.add_argument(
        "--regex-chr",
        default="(.*?)_(.+)_(.*)",
        help="Regex to parse chromosome names (default matches Pf3D7_01_v3 style)",
    )
    p7.add_argument(
        "--regex-group",
        type=int,
        default=3,
        help="Capture group index in regex_chr that contains the numeric chromosome",
    )
    p7.add_argument(
        "--hmmibd-bin",
        default="hmmIBD",
        help="Path/name of hmmIBD binary (e.g. 'hmmIBD' if it’s on PATH)",
    )
    p7.add_argument(
        "--skip-hmmibd",
        action="store_true",
        help="Only prepare hmmIBD input files; do NOT actually run hmmIBD",
    )

    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # HANDLERS
    # -------------------------------------------------------------------------
    if args.command == "missense-drugres-af":
        require_tool("bcftools")
        missense_drugres_af.run(
            vcf=args.vcf,
            ref_fasta=args.ref,
            gff3=args.gff3,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            group_by=args.group_by,
        )

    elif args.command in ("hapmap-africa", "hapmap-samerica", "hapmap-seasia"):
        require_tool("bcftools")
        region = {
            "hapmap-africa": "africa",
            "hapmap-samerica": "south_america",
            "hapmap-seasia": "southeast_asia",
        }[args.command]
        haplotype_map_region.run(
            region=region,
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )

    elif args.command == "fws-dotplot":
        fws_dotplot.run(
            metadata_path=args.metadata,
            outdir=args.outdir,
            group_by=args.group_by,
            width=args.width,
            height=args.height,
        )

    elif args.command == "pca":
        pca_plot.run(
            matrix=args.matrix,
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            group_by=args.group_by,
            pcs=args.pcs,
            max_sample_missing=args.max_sample_missing,
        )

    elif args.command == "hmmibd-matrix":
        # Ensure Rscript and hmmIBD are available
        require_tool("Rscript")
        require_tool(args.hmmibd_bin)
        hmmibg_matrix.run(
            matrix=args.matrix,
            metadata_path=args.metadata,
            workdir=args.workdir,
            category_col=args.category_col,
            category=args.category,
            sample_col=args.sample_col,
            fws_col=args.fws_col,
            fws_th=args.fws_th,
            maf=args.maf,
            na_char=args.na_char,
            threads=args.threads,
            exclude_chr=args.exclude_chr,
            regex_chr=args.regex_chr,
            regex_group=args.regex_group,
            hmmibd_bin=args.hmmibd_bin,
            skip_hmmibd=args.skip_hmmibd,
        )


if __name__ == "__main__":
    main()





