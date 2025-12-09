# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.
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
    # hmmIBD
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "hmmibd-matrix",
        help="Prepare input and run hmmIBD from a binary SNP matrix"
    )

    p.add_argument("--matrix", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", required=True)

    p.add_argument("--category-col", default="country")
    p.add_argument("--category", default=None)
    p.add_argument("--subgroup-col", default=None)

    p.add_argument("--sample-col", default="sample_id")
    p.add_argument("--fws-col", default="fws")
    p.add_argument("--fws-th", type=float, default=0.95)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--na-char", default="N")
    p.add_argument("--threads", type=int, default=4)

    p.add_argument(
        "--exclude-chr",
        default="Pf3D7_API_v3,Pf3D7_MIT_v3",
        help="Highly recommended to exclude organellar chromosomes"
    )

    p.add_argument("--hmmibd-bin", default="hmmIBD")
    p.add_argument("--skip-hmmibd", action="store_true")

    args = parser.parse_args()

    if args.command == "hmmibd-matrix":
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
            hmmibd_bin=args.hmmibd_bin,
            skip_hmmibd=args.skip_hmmibd,
        )


if __name__ == "__main__":
    main()


