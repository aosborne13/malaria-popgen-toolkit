# SPDX-License-Identifier: Apache-2.0
"""
Wrapper around scripts/ibd/summarise_hmmibd_windows.R

Exposes the window summarisation & annotation via:
    malaria-pipeline hmmibd-summary ...
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional


def _repo_root() -> Path:
    """
    Best-effort guess of the repo root assuming an editable install.
    This file is at:
        <repo>/src/malaria_popgen_toolkit/commands/hmmibd_summary.py
    So repo root is parents[3].
    """
    return Path(__file__).resolve().parents[3]


def run(
    workdir: str,
    ref_index: str,
    gene_product: str,
    suffix: str,
    window_size: int = 50000,
    maf: float = 0.01,
    quantile_cutoff: float = 0.95,
    remove_chr: Optional[str] = "Pf3D7_API_v3,Pf3D7_MIT_v3",
    regex_chr: str = "(.*?)_(.+)_(.*)",
    regex_groupid: int = 3,
    rscript_bin: str = "Rscript",
) -> None:
    """
    Call Rscript summarise_hmmibd_windows.R with the provided options.
    """

    repo = _repo_root()
    script = repo / "scripts" / "ibd" / "summarise_hmmibd_windows.R"

    if not script.exists():
        raise FileNotFoundError(
            f"Could not find summarise_hmmibd_windows.R at {script}. "
            "Adjust path in hmmibd_summary.py if your layout is different."
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
        "--maf",
        str(maf),
        "--quantile_cutoff",
        str(quantile_cutoff),
        "--regex_chr",
        regex_chr,
        "--regex_groupid",
        str(regex_groupid),
    ]

    if remove_chr:
        cmd.extend(["--remove_chr", remove_chr])

    print("Running hmmibd window summary with command:\n  " + " ".join(cmd))
    subprocess.run(cmd, check=True)
