# SPDX-License-Identifier: Apache-2.0
"""
Quick dataset sanity checks.
Reports number of samples and SNPs/positions
from either a VCF or a binary genotype matrix.
"""

import subprocess
import pandas as pd
import sys


def _run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"[ERROR] {' '.join(cmd)}\n{p.stderr}")
    return p.stdout.strip()


def _vcf_stats(vcf):
    samples = _run(["bcftools", "query", "-l", vcf]).splitlines()
    n_samples = len(samples)

    n_snps = int(_run(["bcftools", "view", "-H", vcf, "|", "wc", "-l"]))
    return n_samples, n_snps


def _matrix_stats(matrix):
    df = pd.read_csv(matrix, sep="\t", nrows=5)
    n_samples = df.shape[1] - 3  # chr, pos, ref
    n_snps = sum(1 for _ in open(matrix)) - 1  # header
    return n_samples, n_snps


def run(vcf=None, matrix=None):
    if bool(vcf) == bool(matrix):
        sys.exit("Provide exactly ONE of --vcf or --matrix")

    if vcf:
        n_samples, n_snps = _vcf_stats(vcf)
        print("Input type: VCF")

    else:
        n_samples, n_snps = _matrix_stats(matrix)
        print("Input type: matrix")

    print(f"Samples:  {n_samples}")
    print(f"Variants: {n_snps}")
