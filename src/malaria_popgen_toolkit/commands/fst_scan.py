# malaria_popgen_toolkit/commands/fst_scan.py

import subprocess
from pathlib import Path


def run(
    workdir,
    matrix_binary,
    metadata,
    pop1,
    pop2=None,
    label_id="sample_id",
    label_category="country",
    label_fws=None,
    fws_th=0.95,
    estimator="weir-cockerham",
    max_missing=0.2,
    maf=0.01,
    mixed_as="alt",
    na_char="N",
    remove_chr=None,
    gene_file=None,
    fst_thresh=0.95,
    out_prefix=None,
    threads=4,
):
    """
    Genome-wide pairwise Fst between two groups defined in the metadata.
    """

    workdir = Path(workdir)

    cmd = [
        "Rscript",
        str(Path(__file__).parents[2] / "scripts/popgen/fst_scan.R"),
        "--workdir", str(workdir),
        "--matrix_binary", str(matrix_binary),
        "--metadata", str(metadata),
        "--pop1", str(pop1),
        "--label_id", str(label_id),
        "--label_category", str(label_category),
        "--estimator", str(estimator),
        "--max_missing", str(max_missing),
        "--maf", str(maf),
        "--mixed_as", str(mixed_as),
        "--na_char", str(na_char),
        "--fst_thresh", str(fst_thresh),
        "--threads", str(threads),
    ]

    if pop2:
        cmd += ["--pop2", str(pop2)]
    if label_fws:
        cmd += ["--label_fws", str(label_fws), "--fws_th", str(fws_th)]
    if remove_chr:
        cmd += ["--remove_chr", str(remove_chr)]
    if gene_file:
        cmd += ["--gene-file", str(gene_file)]
    if out_prefix:
        cmd += ["--out_prefix", str(out_prefix)]

    subprocess.run(cmd, check=True)
