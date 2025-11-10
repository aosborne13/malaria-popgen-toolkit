# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.
Commands:
  - missense-drugres-af
  - hapmap-africa
  - hapmap-samerica
  - hapmap-seasia
"""

import argparse
import sys
import shutil
from malaria_popgen_toolkit.commands import missense_drugres_af, haplotype_map_region


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
    # Command: missense-drugres-af (VCF + DP filtering)
    # -------------------------------------------------------------------------
    p1 = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes"
    )
    p1.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
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

    # -------------------------------------------------------------------------
    # Command: hapmap-africa (VCF + DP filtering)
    # -------------------------------------------------------------------------
    p2 = sub.add_parser(
        "hapmap-africa",
        help="Plot haplotype pies for Africa (VCF + DP filter) and export per-gene haplotype frequencies"
    )
    p2.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p2.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p2.add_argument("--outdir", default="hapmap_africa_output", help="Output directory")
    p2.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
    p2.add_argument("--sample-col", default="sample_id", help="Metadata column with sample IDs")
    p2.add_argument("--country-col", default="country", help="Metadata column with country")

    # -------------------------------------------------------------------------
    # Command: hapmap-samerica (South America; VCF + DP filtering)
    # -------------------------------------------------------------------------
    p3 = sub.add_parser(
        "hapmap-samerica",
        help="Plot haplotype pies for South America (VCF + DP filter) and export per-gene haplotype frequencies"
    )
    p3.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p3.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p3.add_argument("--outdir", default="hapmap_samerica_output", help="Output directory")
    p3.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
    p3.add_argument("--sample-col", default="sample_id", help="Metadata column with sample IDs")
    p3.add_argument("--country-col", default="country", help="Metadata column with country")

    # -------------------------------------------------------------------------
    # Command: hapmap-seasia (Southeast Asia; VCF + DP filtering)
    # -------------------------------------------------------------------------
    p4 = sub.add_parser(
        "hapmap-seasia",
        help="Plot haplotype pies for Southeast Asia (VCF + DP filter) and export per-gene haplotype frequencies"
    )
    p4.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p4.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p4.add_argument("--outdir", default="hapmap_seasia_output", help="Output directory")
    p4.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
    p4.add_argument("--sample-col", default="sample_id", help="Metadata column with sample IDs")
    p4.add_argument("--country-col", default="country", help="Metadata column with country")

    # -------------------------------------------------------------------------
    # Parse and dispatch
    # -------------------------------------------------------------------------
    args = parser.parse_args()

    # All current commands query VCFs via bcftools
    if args.command in ("missense-drugres-af", "hapmap-africa", "hapmap-samerica", "hapmap-seasia"):
        require_tool("bcftools")

    if args.command == "missense-drugres-af":
        missense_drugres_af.run(
            vcf=args.vcf,
            ref_fasta=args.ref,
            gff3=args.gff3,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            group_by=args.group_by,
        )

    elif args.command == "hapmap-africa":
        haplotype_map_region.run(
            region="africa",
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )

    elif args.command == "hapmap-samerica":
        haplotype_map_region.run(
            region="south_america",
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )

    elif args.command == "hapmap-seasia":
        haplotype_map_region.run(
            region="southeast_asia",
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )


if __name__ == "__main__":
    main()


