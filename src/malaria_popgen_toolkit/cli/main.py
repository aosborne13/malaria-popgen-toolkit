# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.
Currently supports:
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
    if shutil.which(name) is None:
        sys.exit(f"ERROR: '{name}' not found in PATH. Please install it and try again.")


def main():
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # missense-drugres-af
    p1 = sub.add_parser("missense-drugres-af", help="Compute missense allele frequencies in drug-resistance genes")
    p1.add_argument("vcf", help="Input VCF (bgzipped)")
    p1.add_argument("--ref", required=True, help="Reference FASTA file")
    p1.add_argument("--gff3", required=True, help="GFF3 annotation file")
    p1.add_argument("--metadata", required=True, help="TSV with 'sample_id' and grouping column")
    p1.add_argument("--outdir", default="output_missense", help="Output directory")
    p1.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
    p1.add_argument("--group-by", default="country", help="Metadata column to group by (e.g., country, region, site, year)")
    for alias in ("run", "missense-af"):
        pa = sub.add_parser(alias, help="Alias of 'missense-drugres-af'")
        pa.add_argument("vcf", help="Input VCF (bgzipped)")
        pa.add_argument("--ref", required=True, help="Reference FASTA file")
        pa.add_argument("--gff3", required=True, help="GFF3 annotation file")
        pa.add_argument("--metadata", required=True, help="TSV with 'sample_id' and grouping column")
        pa.add_argument("--outdir", default="output_missense", help="Output directory")
        pa.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample read depth (DP)")
        pa.add_argument("--group-by", default="country", help="Metadata column to group by (e.g., country, region, site, year)")

    # hapmap-africa
    p2 = sub.add_parser("hapmap-africa", help="Plot haplotype pies for Africa and export per-gene haplotype frequencies")
    p2.add_argument("--matrix", required=True, help="Binary matrix TSV (chr, pos, ref, then sample columns)")
    p2.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p2.add_argument("--outdir", default="hapmap_africa_output", help="Output directory")
    p2.add_argument("--sample-col", default="sample_id", help="Metadata column with sample IDs")
    p2.add_argument("--country-col", default="country", help="Metadata column with country")

    # hapmap-samerica (South America)
    p3 = sub.add_parser("hapmap-samerica", help="Plot haplotype pies for South America and export per-gene frequencies")
    p3.add_argument("--matrix", required=True, help="Binary matrix TSV (chr, pos, ref, then sample columns)")
    p3.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p3.add_argument("--outdir", default="hapmap_samerica_output", help="Output directory")
    p3.add_argument("--sample-col", default="sample_id", help="Metadata column with sample IDs")
    p3.add_argument("--country-col", default="country", help="Metadata column with country")

    # hapmap-seasia (Southeast Asia)
    p4 = sub.add_parser("hapmap-seasia", help="Plot haplotype pies for Southeast Asia and export per-gene frequencies")
    p4.add_argument("--matrix", required=True, help="Binary matrix TSV (chr, pos, ref, then sample columns)")
    p4.add_argument("--metadata", required=True, help="TSV with sample metadata")
    p4.add_argument("--outdir", default="hapmap_seasia_output", help="Output directory")
    p4.add_argument("--sample-col", default="sample_id", help="Metadata column with sample IDs")
    p4.add_argument("--country-col", default="country", help="Metadata column with country")

    args = parser.parse_args()

    if args.command in ("missense-drugres-af", "run", "missense-af"):
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

    elif args.command == "hapmap-africa":
        haplotype_map_region.run(
            region="africa",
            matrix_path=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )

    elif args.command == "hapmap-samerica":
        haplotype_map_region.run(
            region="south_america",
            matrix_path=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )

    elif args.command == "hapmap-seasia":
        haplotype_map_region.run(
            region="southeast_asia",
            matrix_path=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )


if __name__ == "__main__":
    main()

