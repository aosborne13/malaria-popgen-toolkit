# SPDX-License-Identifier: Apache-2.0
"""
hmmibd-matrix: prepare hmmIBD input from a binary matrix (+ run hmmIBD via R).

This is a Python wrapper around scripts/ibd/run_hmmibd_from_matrix.R
so it can be called from the malaria-pipeline CLI.
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
    matrix: str,
    metadata_path: str,
    workdir: str,
    category_col: str = "country",
    category: str | None = None,
    sample_col: str = "sample_id",
    fws_col: str = "fws",
    fws_th: float = 0.95,
    maf: float = 0.01,
    na_char: str = "NA",
    threads: int = 4,
    exclude_chr: str | None = "Pf3D7_API_v3,Pf3D7_MIT_v3",
    regex_chr: str = "(.*?)_(.+)_(.*)",
    regex_group: int = 3,
    hmmibd_bin: str = "hmmIBD",
    skip_hmmibd: bool = False,
) -> None:
    """
    Orchestrate the Rscript call. All heavy lifting happens in R.

    Parameters mirror the CLI args defined in cli/main.py.
    """

    _require_rscript()

    # Resolve repo root (…/malaria-popgen-toolkit/)
    # __file__ = …/src/malaria_popgen_toolkit/commands/hmmibd_matrix.py
    repo_root = Path(__file__).resolve().parents[3]
    r_script = repo_root / "scripts" / "ibd" / "run_hmmibd_from_matrix.R"

    if not r_script.exists():
        sys.exit(f"ERROR: R script not found at: {r_script}")

    os.makedirs(workdir, exist_ok=True)

    cmd = [
        "Rscript",
        str(r_script),
        "--workdir",
        workdir,
        "--binary_matrix",
        matrix,
        "--metadata",
        metadata_path,
        "--label_category",
        category_col,
        "--label_id",
        sample_col,
        "--label_fws",
        fws_col,
        "--fws_th",
        str(fws_th),
        "--maf",
        str(maf),
        "--na_char",
        na_char,
        "--threads",
        str(threads),
        "--regex_chr",
        regex_chr,
        "--regex_groupid",
        str(regex_group),
        "--hmmibd_bin",
        hmmibd_bin,
    ]

    if category:
        cmd.extend(["--category", category])

    if exclude_chr:
        cmd.extend(["--remove_chr", exclude_chr])

    if skip_hmmibd:
        cmd.append("--skip_hmmibd")

    print("[hmmibd-matrix] Running R pipeline:")
    print("  " + " ".join(cmd) + "\n")

    subprocess.check_call(cmd)

