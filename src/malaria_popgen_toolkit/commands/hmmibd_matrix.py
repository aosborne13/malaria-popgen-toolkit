# SPDX-License-Identifier: Apache-2.0
"""
Run hmmIBD from a genotype matrix using the R script
'run_hmmibd_from_matrix.R'. Produces IBD segments, windowed fractions,
and genomewide summary files compatible with the toolkit's plotting scripts.
"""

import subprocess
import os
import sys
from shutil import which


def require_rscript():
    if which("Rscript") is None:
        sys.exit("ERROR: Rscript not found in PATH")


def run(
    matrix,
    metadata,
    outdir,
    category=None,
    label_category="country",
    label_id="sample",
    label_fws="fws",
    fws_th=0.95,
    maf=0.01,
    rm_chr=None,
    threads=4,
    na_char="NA",
    regex_chr="(.*?)_(.+)_(.*)",
    regex_groupid=3,
    script_path=None,
):
    require_rscript()

    if script_path is None:
        script_path = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "scripts",
            "ibd",
            "run_hmmibd_from_matrix.R"
        )

    if not os.path.exists(script_path):
        sys.exit(f"ERROR: Cannot find R script at: {script_path}")

    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "Rscript",
        script_path,
        "--workdir", outdir,
        "--binary_matrix", matrix,
        "--metadata", metadata,
        "--label_category", label_category,
        "--label_id", label_id,
        "--label_fws", label_fws,
        "--fws_th", str(fws_th),
        "--maf", str(maf),
        "--threads", str(threads),
        "--na_char", na_char,
        "--regex_chr", regex_chr,
        "--regex_groupid", str(regex_groupid),
    ]

    if category:
        cmd += ["--category", category]

    if rm_chr:
        cmd += ["--remove_chr", rm_chr]

    print("Running hmmIBD pipeline…")
    print(" ".join(cmd), "\n")

    subprocess.run(cmd, check=True)
