# SPDX-License-Identifier: Apache-2.0
"""
hmmibd-matrix: prepare hmmIBD input from a binary matrix and run hmmIBD.
"""

from __future__ import annotations
import os
import sys
import subprocess
from pathlib import Path
from shutil import which


def _require_rscript():
    if which("Rscript") is None:
        sys.exit("ERROR: Rscript not found in PATH.")


def run(
    matrix: str,
    metadata_path: str,
    outdir: str,
    category_col: str = "country",
    category: str | None = None,
    sample_col: str = "sample_id",
    fws_col: str = "fws",
    fws_th: float = 0.95,
    maf: float = 0.01,
    threads: int = 4,
    exclude_chr: str = "Pf3D7_API_v3,Pf3D7_MIT_v3",
    hmmibd_bin: str = "hmmIBD",
    skip_hmmibd: bool = False,
):
    _require_rscript()

    repo_root = Path(__file__).resolve().parents[3]
    r_script = repo_root / "scripts" / "ibd" / "run_hmmibd_from_matrix.R"

    if not r_script.exists():
        sys.exit(f"ERROR: R script not found at {r_script}")

    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "Rscript", str(r_script),
        "--outdir", outdir,
        "--binary_matrix", matrix,
        "--metadata", metadata_path,
        "--label_category", category_col,
        "--label_id", sample_col,
        "--label_fws", fws_col,
        "--fws_th", str(fws_th),
        "--maf", str(maf),
        "--threads", str(threads),
        "--hmmibd_bin", hmmibd_bin,
    ]

    if category:
        cmd += ["--category", category]

    if exclude_chr:
        cmd += ["--remove_chr", exclude_chr]

    if skip_hmmibd:
        cmd.append("--skip_hmmibd")

    print("[hmmibd-matrix] Running:\n  " + " ".join(cmd))
    subprocess.check_call(cmd)


