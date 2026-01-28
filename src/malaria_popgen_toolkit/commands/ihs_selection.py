# malaria_popgen_toolkit/commands/ihs_selection.py

from __future__ import annotations

import subprocess
from pathlib import Path


def run(
    workdir,
    matrix_binary,
    metadata_path,
    genome_file,
    label_category="country",
    subgroup_col=None,
    label_id="sample_id",
    label_fws="fws",
    category=None,
    focus_pop=None,
    fws_th=0.95,
    maf=0.01,
    rehh_min_perc_hap=80,
    rehh_min_perc_mrk=70,
    na_char="NA",
    forced_recode=False,
    forced_mixed=False,
    remove_chr=None,
    regex_chr="(.*?)_(.+)_(.*)",
    regex_groupid=3,
    threads=4,
    min_maf_ihs=0.0,
    freqbin=0.05,
    ihs_thresh=2.0,
    logp_thresh=5.0,
):
    """
    Run iHS selection scan using rehh from a binary SNP matrix.

    NOTE:
      - The SNP-level annotation file (--annotation) has been removed.
      - The R script constructs the marker map directly from the matrix
        (chr/pos/ref) and uses the genome product TSV for gene interval mapping.

    Inputs are still file paths:
      - matrix_binary (binary SNP matrix TSV)
      - metadata_path (metadata TSV)
      - genome_file (genome product TSV)
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    matrix_binary = Path(matrix_binary)
    metadata_path = Path(metadata_path)
    genome_file = Path(genome_file)

    # Fail fast with clear errors (helps debugging cache/species resolution)
    for pth, label in (
        (matrix_binary, "matrix_binary"),
        (metadata_path, "metadata_path"),
        (genome_file, "genome_file"),
    ):
        if not pth.exists():
            raise FileNotFoundError(f"{label} not found: {pth}")

    r_script = Path(__file__).parents[2] / "scripts/selection/ihs_full_pipeline.R"
    if not r_script.exists():
        raise FileNotFoundError(f"R script not found: {r_script}")

    cmd = [
        "Rscript",
        str(r_script),
        "--workdir",
        str(workdir),
        "--matrix_binary",
        str(matrix_binary),
        "--metadata",
        str(metadata_path),
        "--genome-file",
        str(genome_file),
        "--label_category",
        str(label_category),
        "--label_id",
        str(label_id),
        "--label_fws",
        str(label_fws),
        "--fws_th",
        str(fws_th),
        "--maf",
        str(maf),
        "--rehh_min_perc_hap",
        str(rehh_min_perc_hap),
        "--rehh_min_perc_mrk",
        str(rehh_min_perc_mrk),
        "--na_char",
        str(na_char),
        "--regex_chr",
        str(regex_chr),
        "--regex_groupid",
        str(regex_groupid),
        "--threads",
        str(threads),
        "--min-maf-ihs",
        str(min_maf_ihs),
        "--freqbin",
        str(freqbin),
        "--ihs-thresh",
        str(ihs_thresh),
        "--logp-thresh",
        str(logp_thresh),
    ]

    if category:
        cmd += ["--category", str(category)]

    if subgroup_col:
        cmd += ["--subgroup_col", str(subgroup_col)]

    if focus_pop:
        cmd += ["--focus-pop", str(focus_pop)]

    if remove_chr:
        cmd += ["--remove_chr", str(remove_chr)]

    if forced_recode:
        cmd.append("--forced_recode")

    if forced_mixed:
        cmd.append("--forced_mixed")

    subprocess.run(cmd, check=True)
