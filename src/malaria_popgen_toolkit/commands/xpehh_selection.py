# malaria_popgen_toolkit/commands/xpehh_selection.py

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

    Automatically compares focus_pop vs all other groups found in workdir.
    """

    workdir = Path(workdir)

    cmd = [
        "Rscript",
        str(Path(__file__).parents[2] / "scripts/selection/xpehh_full_pipeline.R"),
        "--workdir", str(workdir),
        "--genome-file", genome_file,
        "--focus-pop", focus_pop,
        "--min-abs-xpehh", str(min_abs_xpehh),
        "--min-logp", str(min_logp),
        "--regex_chr", regex_chr,
        "--regex_groupid", str(regex_groupid),
    ]

    if remove_chr:
        cmd += ["--remove_chr", remove_chr]

    if panel_groups:
        cmd += ["--panel-groups", panel_groups]

    subprocess.run(cmd, check=True)
