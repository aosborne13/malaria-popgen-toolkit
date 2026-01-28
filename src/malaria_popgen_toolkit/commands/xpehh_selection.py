# malaria_popgen_toolkit/commands/xpehh_selection.py

from __future__ import annotations

import subprocess
from pathlib import Path


def run(
    workdir,
    genome_file,
    focus_pop,
    min_abs_xpehh=2.0,
    min_logp=1.3,
    remove_chr=None,
    regex_chr="(.*?)_(.+)_(.*)",
    regex_groupid=3,
    panel_groups=None,
):
    """
    Run XP-EHH comparisons using rehh from scanned_haplotypes_<group>.tsv files.

    - Automatically compares focus_pop vs all other groups found in workdir.
    - Inputs remain file paths (workdir + genome_file).
    """

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    genome_file = Path(genome_file)

    # Fail fast with clear errors
    if not genome_file.exists():
        raise FileNotFoundError(f"genome_file not found: {genome_file}")

    r_script = Path(__file__).parents[2] / "scripts/selection/xpehh_full_pipeline.R"
    if not r_script.exists():
        raise FileNotFoundError(f"R script not found: {r_script}")

    cmd = [
        "Rscript",
        str(r_script),
        "--workdir",
        str(workdir),
        "--genome-file",
        str(genome_file),
        "--focus-pop",
        str(focus_pop),
        "--min-abs-xpehh",
        str(min_abs_xpehh),
        "--min-logp",
        str(min_logp),
        "--regex_chr",
        str(regex_chr),
        "--regex_groupid",
        str(regex_groupid),
    ]

    if remove_chr:
        cmd += ["--remove_chr", str(remove_chr)]

    if panel_groups:
        cmd += ["--panel-groups", str(panel_groups)]

    subprocess.run(cmd, check=True)

