# SPDX-License-Identifier: Apache-2.0
"""
hmmibd-matrix: prepare and run hmmIBD from a binary SNP matrix.

This wraps scripts/ibd/run_hmmibd_from_matrix.R

It:
- subsets the matrix by category / subgroup (done inside the R script),
- filters samples by Fws >= threshold,
- writes per-category (or category_subgroup) hmmIBD input files,
- runs hmmIBD for each group (unless --skip-hmmibd is used).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from shutil import which


def _require_rscript() -> None:
    if which("Rscript") is None:
        sys.exit("ERROR: 'Rscript' not found in PATH. Install R and ensure Rscript is available.")


def _require_hmmibd(bin_name: str) -> None:
    if which(bin_name) is None:
        sys.exit(f"ERROR: hmmIBD binary '{bin_name}' not found in PATH.")


def run(
    matrix: str,
    metadata_path: str,
    outdir: str,
    category_col: str = "country",
    category: str | None = None,
    subgroup_col: str | None = None,
    sample_col: str = "sample_id",
    fws_col: str = "fws",
    fws_th: float = 0.95,
    maf: float = 0.01,
    na_char: str = "N",
    threads: int = 4,
    exclude_chr: str | None = "Pf3D7_API_v3,Pf3D7_MIT_v3",
    # These are accepted to match the CLI, but currently NOT passed to the R script.
    regex_chr: str = "(.*?)_(.+)_(.*)",
    regex_group: int = 3,
    hmmibd_bin: str = "hmmIBD",
    skip_hmmibd: bool = False,
) -> None:
    """
    Parameters
    ----------
    matrix : str
        Path to binary SNP matrix (chr, pos, ref, <samples...>).
    metadata_path : str
        TSV with at least sample_col, category_col, fws_col.
    outdir : str
        Output directory for hmmIBD inputs and outputs.
    category_col : str
        Metadata column used to define categories (e.g. country).
    category : str or None
        Single category name or comma-separated list. If None or 'all',
        all categories present in metadata are used.
    subgroup_col : str or None
        Optional metadata column for subgrouping within category (e.g. year).
        If provided, groups are "<category>_<subgroup>".
    sample_col : str
        Metadata column containing sample IDs matching matrix columns.
    fws_col : str
        Metadata column containing Fws values.
    fws_th : float
        Minimum Fws threshold; only samples with Fws >= this value are kept.
    maf : float
        Minor allele frequency threshold for SNPs.
    na_char : str
        Character representing missing genotypes in the matrix (also '.', 'NA' treated as missing).
    threads : int
        Number of threads for the R summarisation (not hmmIBD itself).
    exclude_chr : str or None
        Comma-separated chromosome names to exclude (e.g. 'Pf3D7_API_v3,Pf3D7_MIT_v3').
    regex_chr : str
        Reserved for future use (accepts CLI arg but not used here).
    regex_group : int
        Reserved for future use (accepts CLI arg but not used here).
    hmmibd_bin : str
        Name/path of hmmIBD binary.
    skip_hmmibd : bool
        If True, only prepare input files; do not actually run hmmIBD.
    """
    _require_rscript()
    _require_hmmibd(hmmibd_bin)

    # Locate the R script relative to this file:
    # .../malaria-popgen-toolkit/src/malaria_popgen_toolkit/commands/hmmibd_matrix.py
    # repo root = parents[3]
    repo_root = Path(__file__).resolve().parents[3]
    r_script = repo_root / "scripts" / "ibd" / "run_hmmibd_from_matrix.R"

    if not r_script.exists():
        sys.exit(f"ERROR: run_hmmibd_from_matrix.R not found at: {r_script}")

    cmd = [
        "Rscript",
        str(r_script),
        "--outdir", outdir,
        "--binary_matrix", matrix,
        "--metadata", metadata_path,
        "--label_category", category_col,
        "--label_id", sample_col,
        "--label_fws", fws_col,
        "--fws_th", str(fws_th),
        "--maf", str(maf),
        "--na_char", na_char,
        "--threads", str(threads),
        "--hmmibd_bin", hmmibd_bin,
    ]

    # These are passed as --remove_chr to the R script
    if exclude_chr:
        cmd.extend(["--remove_chr", exclude_chr])

    # Category handling: None or 'all' => let the R script auto-detect all categories
    if category and category.lower() != "all":
        cmd.extend(["--category", category])

    # Optional subgrouping (e.g. year)
    if subgroup_col:
        cmd.extend(["--subgroup_col", subgroup_col])

    # Optional flag to skip actually running hmmIBD (R script handles it)
    if skip_hmmibd:
        cmd.append("--skip_hmmibd")

    print("[hmmibd-matrix] Running:")
    print("  " + " ".join(cmd) + "\n")

    subprocess.check_call(cmd)


