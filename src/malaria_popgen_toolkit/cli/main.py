# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.
Commands:
  - missense-drugres-af   (AF of missense drug-res variants with DP filter; flexible grouping)
  - hapmap-africa         (haplotype map for Africa)
  - hapmap-samerica       (haplotype map for South America)
  - hapmap-seasia         (haplotype map for Southeast Asia)
  - fws-dotplot           (Fws jittered dot plots by metadata groups)
"""

import argparse
import sys
import shutil

from malaria_popgen_toolkit.commands import (
    missense_drugres_af,
    haplotype_map_region,
    fws_dotplot,
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
    # hapmap-africa / -samerica / -seasia
    # -------------------------------------------------------------------------
    p2 = sub.add_parser(
        "hapmap-africa",
        help="Plot haplotype pies for Africa from VCF + metadata"
    )
    p2.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p2.add_argument("--metadata", required=True, help="TSV with samples and country column")
    p2.add_argument("--outdir", required=True, help="Output directory for map + CSVs")
    p2.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP for hap calls")
    p2.add_argument("--sample-col", default="sample_id", help="Sample column name in metadata")
    p2.add_argument("--country-col", default="country", help="Country column name in metadata")

    p3 = sub.add_parser(
        "hapmap-samerica",
        help="Plot haplotype pies for South America from VCF + metadata"
    )
    p3.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p3.add_argument("--metadata", required=True, help="TSV with samples and country column")
    p3.add_argument("--outdir", required=True, help="Output directory for map + CSVs")
    p3.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP for hap calls")
    p3.add_argument("--sample-col", default="sample_id", help="Sample column name in metadata")
    p3.add_argument("--country-col", default="country", help="Country column name in metadata")

    p4 = sub.add_parser(
        "hapmap-seasia",
        help="Plot haplotype pies for Southeast Asia from VCF + metadata"
    )
    p4.add_argument("--vcf", required=True, help="Input VCF (bgzipped)")
    p4.add_argument("--metadata", required=True, help="TSV with samples and country column")
    p4.add_argument("--outdir", required=True, help="Output directory for map + CSVs")
    p4.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP for hap calls")
    p4.add_argument("--sample-col", default="sample_id", help="Sample column name in metadata")
    p4.add_argument("--country-col", default="country", help="Country column name in metadata")

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
        help="Metadata column to group by (repeat for multiple). Default tries: region, country, year."
    )
    p5.add_argument("--width", type=float, default=10.0, help="Plot width (inches)")
    p5.add_argument("--height", type=float, default=6.0, help="Plot height (inches)")

    args = parser.parse_args()

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


if __name__ == "__main__":
    main()


