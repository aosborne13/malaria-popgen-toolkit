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
    p1 = sub.add_parser("missense-drugres-af")
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
        p = sub.add_parser(cmd)
        p.add_argument("--vcf", required=True)
        p.add_argument("--metadata", required=True)
        p.add_argument("--outdir", required=True)
        p.add_argument("--min-dp", type=int, default=5)
        p.add_argument("--sample-col", default="sample_id")
        p.add_argument("--country-col", default="country")

    # ------------------------------------------------------------------
    # Fws dotplot
    # ------------------------------------------------------------------
    p5 = sub.add_parser("fws-dotplot")
    p5.add_argument("--metadata", required=True)
    p5.add_argument("--outdir", default="fws_plots")
    p5.add_argument("--group-by", action="append")
    p5.add_argument("--width", type=float, default=10)
    p5.add_argument("--height", type=float, default=6)

    # ------------------------------------------------------------------
    # PCA / PCoA
    # ------------------------------------------------------------------
    p6 = sub.add_parser("pca")
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
    p7 = sub.add_parser("hmmibd-matrix")
    p7.add_argument("--matrix", required=True)
    p7.add_argument("--metadata", required=True)
    p7.add_argument("--outdir", required=True)
    p7.add_argument("--category-col", default="country")
    p7.add_argument("--category", default=None)
    p7.add_argument("--sample-col", default="sample_id")
    p7.add_argument("--fws-col", default="fws")
    p7.add_argument("--fws-th", type=float, default=0.95)
    p7.add_argument("--maf", type=float, default=0.01)
    p7.add_argument("--threads", type=int, default=4)
    p7.add_argument("--exclude-chr", default="Pf3D7_API_v3,Pf3D7_MIT_v3")
    p7.add_argument("--hmmibd-bin", default="hmmIBD")
    p7.add_argument("--skip-hmmibd", action="store_true")

    args = parser.parse_args()

    if args.command == "missense-drugres-af":
        require_tool("bcftools")
        missense_drugres_af.run(**vars(args))

    elif args.command in ("hapmap-africa", "hapmap-samerica", "hapmap-seasia"):
        require_tool("bcftools")
        region = {
            "hapmap-africa": "africa",
            "hapmap-samerica": "south_america",
            "hapmap-seasia": "southeast_asia",
        }[args.command]
        haplotype_map_region.run(region=region, **vars(args))

    elif args.command == "fws-dotplot":
        fws_dotplot.run(**vars(args))

    elif args.command == "pca":
        pca_plot.run(**vars(args))

    elif args.command == "hmmibd-matrix":
        require_tool("Rscript")
        require_tool(args.hmmibd_bin)
        hmmibd_matrix.run(
            matrix=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            category_col=args.category_col,
            category=args.category,
            sample_col=args.sample_col,
            fws_col=args.fws_col,
            fws_th=args.fws_th,
            maf=args.maf,
            threads=args.threads,
            exclude_chr=args.exclude_chr,
            hmmibd_bin=args.hmmibd_bin,
            skip_hmmibd=args.skip_hmmibd,
        )


if __name__ == "__main__":
    main()




