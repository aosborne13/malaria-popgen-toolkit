# malaria_popgen_toolkit/commands/nj_tree.py

import subprocess
from pathlib import Path


def run(
    workdir,
    matrix_binary,
    metadata,
    label_id="sample_id",
    label_fws=None,
    fws_th=0.95,
    region=None,
    remove_chr=None,
    mixed_call="alt",
    max_sample_missing=0.2,
    max_locus_missing=0.2,
    min_shared=50,
    max_loci=None,
    n_boot=0,
    na_char="N",
    seed=1,
    prefix="nj_tree",
    plot=False,
):
    """
    Neighbour-joining tree from a binary SNP matrix.

    Manhattan distances via amap::Dist (the metric pca-plot also uses), tree via
    ape::nj, optional bootstrap via ape::boot.phylo.
    """

    workdir = Path(workdir)

    cmd = [
        "Rscript",
        str(Path(__file__).resolve().parents[3] / "scripts/nj_tree.R"),
        "--workdir", str(workdir),
        "--matrix_binary", str(matrix_binary),
        "--metadata", str(metadata),
        "--label_id", str(label_id),
        "--mixed_call", str(mixed_call),
        "--max_sample_missing", str(max_sample_missing),
        "--max_locus_missing", str(max_locus_missing),
        "--min_shared", str(min_shared),
        "--n_boot", str(n_boot),
        "--na_char", str(na_char),
        "--seed", str(seed),
        "--prefix", str(prefix),
    ]

    if label_fws:
        cmd += ["--label_fws", str(label_fws), "--fws_th", str(fws_th)]
    if region:
        cmd += ["--region", str(region)]
    if remove_chr:
        cmd += ["--remove_chr", str(remove_chr)]
    if max_loci:
        cmd += ["--max_loci", str(max_loci)]
    if plot:
        cmd += ["--plot"]

    subprocess.run(cmd, check=True)
