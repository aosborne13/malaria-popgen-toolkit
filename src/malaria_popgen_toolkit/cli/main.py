# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.

Commands:
  - dataset-stats          Report number of samples and SNPs (VCF or matrix)
  - missense-drugres-af    Allele frequencies of missense drug-resistance variants
  - hapmap-africa          Haplotype map for Africa
  - hapmap-samerica        Haplotype map for South America
  - hapmap-seasia          Haplotype map for Southeast Asia
  - fws-dotplot            Fws jittered dot plots by metadata groups
  - pca                    Distance-based PCA / PCoA from binary matrix
  - hmmibd-matrix          Prepare and run hmmIBD from a binary SNP matrix
"""

import argparse
import sys
import shutil

from malaria_popgen_toolkit.commands import (
    dataset_stats,
    missense_drugres_af,
    haplotype_map_region,
    fws_dotplot,
    pca_plot,
    hmmibd_matrix,
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
    # dataset-stats
    # ------------------------------------------------------------------
    p0 = sub.add_parser(
        "dataset-stats",
        help="Report number of samples and SNPs from a VCF or matrix"
    )
    p0.add_argument("--vcf", help="Input VCF (bgzipped)")
    p0.add_argument("--matrix", help="Genotype matrix (.tsv / .bin)")

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
    p5.add_argument("--group-by", action="append")
    p5.add_argument("--width", type=float, default=10)
    p5.add_argument("--height", type=float, default=6)

    # ------------------------------------------------------------------
    # PCA / PCoA
    # ------------------------------------------------------------------
    p6 = sub.add_parser("pca", help="Distance-based PCA / PCoA")
    p6.add_argument("--matrix", required=True)
    p6.add_argument("--metadata", required=True)
    p6.add_argument("--outdir", default="pca_plots")
    p6.add_argument("--sample-col", default="sample_id")
    p6.add_argument("--group-by", nargs="+")
    p6.add_argument("--pcs", nargs="+")
    p6.add_argument("--max-sample-missing", type=float)

    # ------------------------------------------------------------------
    # hmmIBD
    # ------------------------------------------------------------------
    p7 = sub.add_parser(
        "hmmibd-matrix",
        help="Prepare input and run hmmIBD from a binary SNP matrix"
    )
    p7.add_argument("--matrix", required=True)
    p7.add_argument("--metadata", required=True)
    p7.add_argument("--workdir", required=True)
    p7.add_argument("--category-col", default="country")
    p7.add_argument("--category", default=None)
    p7.add_argument("--sample-col", default="sample_id")
    p7.add_argument("--fws-col", default="fws")
    p7.add_argument("--fws-th", type=float, default=0.95)
    p7.add_argument("--maf", type=float, default=0.01)
    p7.add_argument("--threads", type=int, default=4)
    p7.add_argument("--hmmibd-bin", default="hmmIBD")
    p7.add_argument("--skip-hmmibd", action="store_true")

    args = parser.parse_args()

    # ==========================
    # Dispatch
    # ==========================

    if args.command == "dataset-stats":
        if args.vcf:
            require_tool("bcftools")
        dataset_stats.run(vcf=args.vcf, matrix=args.matrix)

    elif args.command == "missense-drugres-af":
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
            workdir=args.workdir,
            category_col=args.category_col,
            category=args.category,
            sample_col=args.sample_col,
            fws_col=args.fws_col,
            fws_th=args.fws_th,
            maf=args.maf,
            threads=args.threads,
            hmmibd_bin=args.hmmibd_bin,
            skip_hmmibd=args.skip_hmmibd,
        )


if __name__ == "__main__":
    main()







