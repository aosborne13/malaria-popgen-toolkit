# SPDX-License-Identifier: Apache-2.0
"""
Report basic dataset statistics:
- Number of samples
- Number of variants (SNPs)

Supports:
- bgzipped VCF via bcftools
- binary / TSV genotype matrix
"""

from __future__ import annotations
import subprocess
import sys
import pandas as pd


def _run(cmd: list[str]) -> int:
    p = subprocess.run(cmd, shell=False, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"[ERROR] {' '.join(cmd)}\n{p.stderr}")
    try:
        return int(p.stdout.strip())
    except ValueError:
        sys.exit(f"[ERROR] Unexpected output from command: {' '.join(cmd)}")


def stats_vcf(vcf: str):
    # number of samples
    nsamples = _run(["bash", "-c", f"bcftools query -l {vcf} | wc -l"])

    # number of variants (non-header lines)
    nvars = _run(["bash", "-c", f"bcftools view -H {vcf} | wc -l"])

    print("Input type: VCF")
    print(f"Samples:  {nsamples}")
    print(f"Variants: {nvars}")


def stats_matrix(matrix: str):
    # read header only
    with open(matrix) as fh:
        header = fh.readline().rstrip("\n").split("\t")

    # first three columns assumed: chrom, pos, ref
    nsamples = len(header) - 3

    # line count minus header
    nvars = sum(1 for _ in open(matrix)) - 1

    print("Input type: matrix")
    print(f"Samples:  {nsamples}")
    print(f"Variants: {nvars}")


def run(vcf: str | None = None, matrix: str | None = None):
    if not vcf and not matrix:
        sys.exit("ERROR: Must provide --vcf or --matrix")

    if vcf and matrix:
        sys.exit("ERROR: Provide only one of --vcf or --matrix")

    if vcf:
        stats_vcf(vcf)
    else:
        stats_matrix(matrix)
