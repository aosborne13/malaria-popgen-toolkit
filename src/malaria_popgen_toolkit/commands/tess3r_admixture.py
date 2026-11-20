# malaria_popgen_toolkit/commands/tess3r_admixture.py

import subprocess
from pathlib import Path


def run(
    workdir,
    metadata,
    matrix_binary=None,
    matrix_rds=None,
    label_id="sample_id",
    lon_col="longitude",
    lat_col="latitude",
    group_col="country",
    label_fws=None,
    fws_th=0.95,
    K="1-10",
    rep=10,
    method="projected.ls",
    lambda_=1.0,
    max_iteration=200,
    ploidy=1,
    na_char="N",
    remove_chr=None,
    prefix="tess3r",
    force_rebuild=False,
    threads=4,
):
    """
    Spatially-explicit ancestry estimation with TESS3.

    A range of K (e.g. "1-10") writes the cross-entropy curve used to choose K;
    a single K (e.g. "6") writes the ancestry coefficients, barplot and map.
    """

    workdir = Path(workdir)

    cmd = [
        "Rscript",
        str(Path(__file__).resolve().parents[3] / "scripts/tess3r_admixture.R"),
        "--workdir", str(workdir),
        "--metadata", str(metadata),
        "--label_id", str(label_id),
        "--lon_col", str(lon_col),
        "--lat_col", str(lat_col),
        "--group_col", str(group_col),
        "--K", str(K),
        "--rep", str(rep),
        "--method", str(method),
        "--lambda", str(lambda_),
        "--max_iteration", str(max_iteration),
        "--ploidy", str(ploidy),
        "--na_char", str(na_char),
        "--prefix", str(prefix),
        "--threads", str(threads),
    ]

    if matrix_binary:
        cmd += ["--matrix_binary", str(matrix_binary)]
    if matrix_rds:
        cmd += ["--matrix_rds", str(matrix_rds)]
    if label_fws:
        cmd += ["--label_fws", str(label_fws), "--fws_th", str(fws_th)]
    if remove_chr:
        cmd += ["--remove_chr", str(remove_chr)]
    if force_rebuild:
        cmd += ["--force_rebuild"]

    subprocess.run(cmd, check=True)
