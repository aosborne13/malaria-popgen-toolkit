# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.

Commands:
  - missense-drugres-af   Allele frequencies of missense drug-resistance variants
  - hapmap-africa         Haplotype map for Africa
  - hapmap-samerica       Haplotype map for South America
  - hapmap-seasia         Haplotype map for Southeast Asia
  - fws-dotplot           Fws jittered dot plots by metadata groups
  - pca                   Distance-based PCA / PCoA from binary matrix
  - hmmibd-matrix         Prepare and run hmmIBD from a binary SNP matrix
  - dataset-stats         Quick sample/variant counts for VCF or matrix
  - hmmibd-ibdplot        Plot IBD summaries (boxplot, genome-wide, chromosome painting)
"""

import argparse
import sys
import shutil

from malaria_popgen_toolkit.commands import (
    missense_drugres_af,
    haplotype_map_region,
    fws_dotplot,
    pca_plot,
    hmmibd_matrix,
    dataset_stats,
    ibd_plot,
)


def require_tool(name: str):
    if shutil.which(name) is None:
        sys.exit(f"ERROR: Required tool '{name}' not found in PATH.")


def main():
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ------------------------------------------------------------------
    # missense-drugres-af
    # ------------------------------------------------------------------
    p1 = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes"
    )
    p1.add_argument("--vcf", required=True)
    p1.add_argument("--ref", required=True)
    p1.add_argument("--gff3", required=True)
    p1.add_argument("--metadata", required=True)
    p1.add_argument("--outdir", default="missense_af")
    p1.add_argument("--min-dp", type=int, default=5)
    p1.add_argument("--group-by", default="country")

    # ------------------------------------------------------------------
    # Haplotype maps
    # ------------------------------------------------------------------
    for cmd, region in [
        ("hapmap-africa", "africa"),
        ("hapmap-samerica", "south_america"),
        ("hapmap-seasia", "southeast_asia"),
    ]:
        p = sub.add_parser(cmd, help=f"Haplotype map for {region}")
        p.add_argument("--vcf", required=True)
        p.add_argument("--metadata", required=True)
        p.add_argument("--outdir", required=True)
        p.add_argument("--min-dp", type=int, default=5)
        p.add_argument("--sample-col", default="sample_id")
        p.add_argument("--country-col", default="country")

    # ------------------------------------------------------------------
    # Fws dotplot
    # ------------------------------------------------------------------
    p5 = sub.add_parser("fws-dotplot", help="Fws jittered dot plots")
    p5.add_argument("--metadata", required=True)
    p5.add_argument("--outdir", default="fws_plots")
    p5.add_argument("--group-by", action="append",
                    help="Metadata column(s) to group by; repeat for multiple. "
                         "If omitted, tries region, country, year.")
    p5.add_argument("--width", type=float, default=10)
    p5.add_argument("--height", type=float, default=6)

    # ------------------------------------------------------------------
    # PCA / PCoA
    # ------------------------------------------------------------------
    p6 = sub.add_parser("pca", help="Distance-based PCA / PCoA")
    p6.add_argument("--matrix", required=True,
                    help="Binary genotype matrix (.tsv) with samples as columns")
    p6.add_argument("--metadata", required=True)
    p6.add_argument("--outdir", default="pca_plots")
    p6.add_argument("--sample-col", default="sample_id")
    p6.add_argument("--group-by", nargs="+",
                    help="Metadata columns to color by (e.g. country region year)")
    p6.add_argument("--pcs", nargs="+",
                    help="PC pairs to plot, e.g. --pcs 1,2 1,3 (default: 1,2 and 1,3)")
    p6.add_argument("--max-sample-missing", type=float,
                    help="Drop samples with >X fraction missing (e.g. 0.3)")

    # ------------------------------------------------------------------
    # hmmIBD (matrix -> input + run)
    # ------------------------------------------------------------------
    p7 = sub.add_parser(
        "hmmibd-matrix",
        help="Prepare input and run hmmIBD from a binary SNP matrix"
    )
    p7.add_argument("--matrix", required=True,
                    help="Binary SNP matrix with columns: chr,pos,ref,<samples>")
    p7.add_argument("--metadata", required=True,
                    help="Metadata TSV with sample IDs, category, and Fws column")
    p7.add_argument("--outdir", required=True,
                    help="Working/output directory for hmmIBD inputs/outputs")
    p7.add_argument("--category-col", default="country")
    p7.add_argument("--category", default=None,
                    help="Category name or comma-separated list (e.g. Ethiopia,Kenya). "
                         "If omitted, run all categories.")
    p7.add_argument("--subgroup-col", default=None,
                    help="Optional second-level grouping (e.g. year) to split within category")
    p7.add_argument("--sample-col", default="sample_id")
    p7.add_argument("--fws-col", default="fws")
    p7.add_argument("--fws-th", type=float, default=0.95)
    p7.add_argument("--maf", type=float, default=0.01)
    p7.add_argument("--threads", type=int, default=4)
    p7.add_argument("--hmmibd-bin", default="hmmIBD")
    p7.add_argument("--skip-hmmibd", action="store_true")
    p7.add_argument("--na-char", default="N",
                    help="Character used for missing genotypes in matrix "
                         "(also '.' and 'N' treated as missing)")
    p7.add_argument("--exclude-chr", default="Pf3D7_API_v3,Pf3D7_MIT_v3",
                    help="Comma-separated list of chromosomes to drop "
                         "(e.g. Pf3D7_API_v3,Pf3D7_MIT_v3)")
    p7.add_argument("--regex-chr", default="(.*?)_(.+)_(.*)",
                    help="Regex to parse chromosome names, default matches Pf3D7_01_v3")
    p7.add_argument("--regex-group", type=int, default=3,
                    help="Capture group index in regex_chr that contains numeric chromosome")

    # ------------------------------------------------------------------
    # dataset-stats
    # ------------------------------------------------------------------
    p8 = sub.add_parser(
        "dataset-stats",
        help="Quickly report #samples and #variants from a VCF or binary matrix"
    )
    p8.add_argument("--vcf", help="VCF/BCF (bgzipped) file")
    p8.add_argument("--matrix", help="Binary matrix (.tsv) with samples as columns")

    # ------------------------------------------------------------------
    # hmmibd-ibdplot: plotting the IBD summaries
    # ------------------------------------------------------------------
    p9 = sub.add_parser(
        "hmmibd-ibdplot",
        help="Plot IBD summaries (boxplot, genome-wide, chromosome painting) from hmmIBD outputs"
    )
    p9.add_argument("--workdir", required=True,
                    help="Directory containing <suffix>_hmmIBD_ibd_win50kb.tsv "
                         "and <suffix>_hmmIBD_fraction.tsv")
    p9.add_argument("--metadata", required=True,
                    help="Metadata TSV used for grouping (e.g. Fws file)")
    p9.add_argument("--ref-fai", required=True,
                    help="Reference FASTA .fai index (Pf3D7)")
    p9.add_argument("--gene-file", required=True,
                    help="Pf3D7 genome product TSV annotation")
    p9.add_argument("--suffix", required=True,
                    help="Prefix used in hmmIBD summarisation, e.g. '13_08_2025'")
    p9.add_argument("--group-var", default="collection",
                    help="Metadata column used as category (e.g. collection, year) "
                         "[default: collection]")
    p9.add_argument("--remove-chr", default="Pf3D7_API_v3,Pf3D7_MIT_v3",
                    help="Comma-separated chromosomes to remove from FAI [default: Pf3D7_API_v3,Pf3D7_MIT_v3]")
    p9.add_argument("--outdir", default="win_50kb",
                    help="Subdirectory inside workdir for plots/tables [default: win_50kb]")

    args = parser.parse_args()

    # ==========================
    # Dispatch
    # ==========================

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
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            group_by=args.group_by,
            pcs=args.pcs,
            max_sample_missing=args.max_sample_missing,
        )

    elif args.command == "hmmibd-matrix":
        require_tool("Rscript")
        require_tool(args.hmmibd_bin)
        hmmibd_matrix.run(
            matrix=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            category_col=args.category_col,
            category=args.category,
            subgroup_col=args.subgroup_col,
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

    elif args.command == "dataset-stats":
        dataset_stats.run(
            vcf=args.vcf,
            matrix=args.matrix,
        )

    elif args.command == "hmmibd-ibdplot":
        require_tool("Rscript")
        ibd_plot.run(
            workdir=args.workdir,
            metadata_path=args.metadata,
            ref_fai=args.ref_fai,
            gene_file=args.gene_file,
            suffix=args.suffix,
            group_var=args.group_var,
            remove_chr=args.remove_chr,
            outdir=args.outdir,
        )


if __name__ == "__main__":
    main()


