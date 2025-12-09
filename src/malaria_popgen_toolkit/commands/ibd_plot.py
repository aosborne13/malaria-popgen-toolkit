# SPDX-License-Identifier: Apache-2.0
"""
ibd_plot: Python wrapper for the hmmIBD IBD plotting R script.

This calls:
  scripts/ibd/plot_hmmibd_ibd.R

which expects:
  --workdir       directory with *_hmmIBD_ibd_win50kb.tsv and *_hmmIBD_fraction.tsv
  --metadata      metadata TSV with grouping column (e.g. collection/year)
  --ref_fai       Pf3D7 FASTA .fai index
  --gene_file     Pf3D7 genome product TSV
  --suffix        prefix used for the hmmIBD summary files
  --group_var     metadata column used as 'category'
  --remove_chr    chromosomes to drop from FAI
  --output_dir    subdirectory under workdir for plots/tables
"""

from __future__ import annotations

import os
import sys
import subprocess
from pathlib import Path
from shutil import which


def _require_rscript() -> None:
    if which("Rscript") is None:
        sys.exit("ERROR: 'Rscript' not found in PATH. Install R and ensure Rscript is available.")


def run(
    workdir: str,
    metadata_path: str,
    ref_fai: str,
    gene_file: str,
    suffix: str,
    group_var: str = "collection",
    remove_chr: str = "Pf3D7_API_v3,Pf3D7_MIT_v3",
    outdir: str = "win_50kb",
) -> None:
    """
    Thin wrapper that forwards arguments to plot_hmmibd_ibd.R.
    """

    _require_rscript()

    # Locate the R script relative to this file:
    #   src/malaria_popgen_toolkit/commands/ibd_plot.py
    #   -> repo_root = parents[3]
    repo_root = Path(__file__).resolve().parents[3]
    r_script = repo_root / "scripts" / "ibd" / "plot_hmmibd_ibd.R"

    if not r_script.exists():
        sys.exit(f"ERROR: R script not found at: {r_script}")

    # Ensure workdir exists (R will create outdir inside it)
    os.makedirs(workdir, exist_ok=True)

    cmd = [
        "Rscript",
        str(r_script),
        "--workdir",
        workdir,
        "--metadata",
        metadata_path,
        "--ref_fai",
        ref_fai,
        "--gene_file",
        gene_file,
        "--suffix",
        suffix,
        "--group_var",
        group_var,
        "--remove_chr",
        remove_chr,
        "--output_dir",  # R script uses --output_dir internally
        outdir,
    ]

    print("[hmmibd-ibdplot] Running:")
    print("  " + " ".join(cmd) + "\n")

    subprocess.check_call(cmd)
