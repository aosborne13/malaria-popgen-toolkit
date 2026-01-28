# malaria_popgen_toolkit/commands/ihs_selection.py

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
      - The R script now constructs the marker map directly from the matrix
        (chr/pos/ref) and uses the genome product TSV for gene interval mapping.
    """

    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "Rscript",
        str(Path(__file__).parents[2] / "scripts/selection/ihs_full_pipeline.R"),
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
