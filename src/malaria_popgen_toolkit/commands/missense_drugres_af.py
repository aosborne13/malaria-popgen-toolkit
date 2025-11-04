# SPDX-License-Identifier: Apache-2.0
import argparse
import shutil
import sys

# import the renamed command
from malaria_popgen_toolkit.commands import missense_drugres_af

def require_tool(name: str):
    if shutil.which(name) is None:
        sys.exit(f"ERROR: '{name}' not found in PATH. Please install it and try again.")

def main():
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # Preferred, descriptive name
    p1 = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes by country"
    )
    p1.add_argument("vcf", help="Input VCF (bgzipped)")
    p1.add_argument("--ref", required=True, help="Reference FASTA")
    p1.add_argument("--gff3", required=True, help="GFF3 annotation")
    p1.add_argument("--metadata", required=True, help="TSV with 'sample_id' and 'country'")
    p1.add_argument("--outdir", default="output_missense", help="Output directory")
    p1.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP")

    # Back-compat aliases (all route to the same function)
    for alias in ("run", "missense-af"):
        pa = sub.add_parser(alias, help=f"Alias of 'missense-drugres-af'")
        pa.add_argument("vcf", help="Input VCF (bgzipped)")
        pa.add_argument("--ref", required=True, help="Reference FASTA")
        pa.add_argument("--gff3", required=True, help="GFF3 annotation")
        pa.add_argument("--metadata", required=True, help="TSV with 'sample_id' and 'country'")
        pa.add_argument("--outdir", default="output_missense", help="Output directory")
        pa.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP")

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
        )
