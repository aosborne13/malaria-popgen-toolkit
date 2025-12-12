# SPDX-License-Identifier: Apache-2.0
"""
Wrapper around scripts/ibd/plot_hmmibd_ibd.R

Exposes IBD plots via:
    malaria-pipeline hmmibd-ibdplots ...
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def run(
    workdir: str,
    ref_index: str,
    gene_product: str,
    suffix: str,
    window_size: int = 50000,
    quantile_cutoff: float = 0.95,
    remove_chr: Optional[str] = "Pf3D7_API_v3,Pf3D7_MIT_v3",
    regex_chr: str = "(.*?)_(.+)_(.*)",
    regex_groupid: int = 3,
    outdir: Optional[str] = None,
    rscript_bin: str = "Rscript",
) -> None:
    """
    Call Rscript plot_hmmibd_ibd.R with the provided options.
    """

    repo = _repo_root()
    script = repo / "scripts" / "ibd" / "plot_hmmibd_ibd.R"

    if not script.exists():
        raise FileNotFoundError(
            f"Could not find plot_hmmibd_ibd.R at {script}. "
            "Adjust path in hmmibd_ibdplots.py if your layout is different."
        )

    cmd = [
        rscript_bin,
        str(script),
        "--workdir",
        workdir,
        "--ref_index",
        ref_index,
        "--gene_product",
        gene_product,
        "--suffix",
        suffix,
        "--window_size",
        str(window_size),
        "--quantile_cutoff",
        str(quantile_cutoff),
        "--regex_chr",
        regex_chr,
        "--regex_groupid",
        str(regex_groupid),
    ]

    if remove_chr:
        cmd.extend(["--remove_chr", remove_chr])
    if outdir:
        cmd.extend(["--outdir", outdir])

    print("Running hmmIBD plotting with command:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True)
