# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.
Currently supports: missense-drugres-af
"""

import argparse
import sys
import shutil
from malaria_popgen_toolkit.commands import missense_drugres_af


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
    # Command: missense-drugres-af
    # -------------------------------------------------------------------------
    p1 = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes"
    )
    p1.add_argument("vcf", help="Input VCF (bgzipped)")
    p1.add_argument("--ref", required=True, help="Reference FASTA file")
    p1.add_argument("--gff3", required=True, help="GFF3 annotation file")
    p1.add_argument("--metadata", required=True, help="TSV with 'sample_id' and grouping column")
    p1.add_argument("--outdir", default="output_missense", help="Output directory")
    p1.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
    p1.add_argument(
        "--group-by",
        default="country",
        help="Metadata column to group by (e.g., country, region, site, year)"
    )

    # Aliases for backwards compatibility
    for alias in ("run", "missense-af"):
        pa = sub.add_parser(alias, help="Alias of 'missense-drugres-af'")
        pa.add_argument("vcf", help="Input VCF (bgzipped)")
        pa.add_argument("--ref", required=True, help="Reference FASTA file")
        pa.add_argument("--gff3", required=True, help="GFF3 annotation file")
        pa.add_argument("--metadata", required=True, help="TSV with 'sample_id' and grouping column")
        pa.add_argument("--outdir", default="output_missense", help="Output directory")
        pa.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
        pa.add_argument(
            "--group-by",
            default="country",
            help="Metadata column to group by (e.g., country, region, site, year)"
        )

    args = parser.parse_args()
    require_tool("bcftools")

    if args.command in ("missense-drugres-af", "run", "missense-af"):
        missense_drugres_af.run(
            vcf=args.vcf,
            ref_fasta=args.ref,
            gff3=args.gff3,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            group_by=args.group_by,
        )


if __name__ == "__main__":
    main()

