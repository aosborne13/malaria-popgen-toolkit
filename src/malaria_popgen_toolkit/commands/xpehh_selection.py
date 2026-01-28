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

    NOTE:
      - CLI uses --focus-pop; the R script accepts --ref_pop (and we support both).
      - This wrapper passes --ref_pop for compatibility with the R script you pasted.
    """
    workdir = Path(workdir)

    cmd = [
        "Rscript",
        str(Path(__file__).parents[2] / "scripts/selection/xpehh_full_pipeline.R"),
        "--workdir", str(workdir),
        "--genome-file", str(genome_file),
        "--ref_pop", str(focus_pop),
        "--xpehh-thresh", str(min_abs_xpehh),
        "--logp-thresh", str(min_logp),
        "--regex_chr", str(regex_chr),
        "--regex_groupid", str(regex_groupid),
    ]

    if remove_chr:
        cmd += ["--remove_chr", str(remove_chr)]

    if panel_groups:
        # Your R script uses --panel_comparisons; keep CLI arg name but map it here
        cmd += ["--panel_comparisons", str(panel_groups)]

    subprocess.run(cmd, check=True)


