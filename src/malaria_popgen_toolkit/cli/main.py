# SPDX-License-Identifier: Apache-2.0
import argparse
import shutil
import sys
from malaria_popgen_toolkit.commands import missense_af

def require_tool(name: str):
    if shutil.which(name) is None:
        sys.exit(f"ERROR: '{name}' not found in PATH. Please install it and try again.")

def main():
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # Primary command name
    p1 = sub.add_parser("missense-af", help="Compute missense allele frequencies by country")
    p1.add_argument("vcf", help="Input VCF (bgzipped)")
    p1.add_argument("--ref", required=True, help="Reference FASTA")
    p1.add_argument("--gff3", required=True, help="GFF3 annotation")
    p1.add_argument("--metadata", required=True, help="TSV with 'sample_id' and 'country'")
    p1.add_argument("--outdir", default="output_missense", help="Output directory")
    p1.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP")

    # Back-compat alias so README 'run' keeps working (routes to the same function)
    p2 = sub.add_parser("run", help="Alias of 'missense-af' for early README examples")
    p2.add_argument("vcf", help="Input VCF (bgzipped)")
    p2.add_argument("--ref", required=True, help="Reference FASTA")
    p2.add_argument("--gff3", required=True, help="GFF3 annotation")
    p2.add_argument("--metadata", required=True, help="TSV with 'sample_id' and 'country'")
    p2.add_argument("--outdir", default="output_missense", help="Output directory")
    p2.add_argument("--min-dp", type=int, default=5, help="Minimum per-sample DP")

    args = parser.parse_args()

    # external tool requirement
    require_tool("bcftools")

    if args.command in ("missense-af", "run"):
        missense_af.run(
            vcf=args.vcf,
            ref_fasta=args.ref,
            gff3=args.gff3,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
        )
