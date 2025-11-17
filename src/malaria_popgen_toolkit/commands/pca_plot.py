# SPDX-License-Identifier: Apache-2.0
"""
Distance-based PCA/PCoA-style plots from either:
  - a binary matrix TSV (0 / 0.5 / 1, 'N' as missing), or
  - a multi-sample VCF (will be converted to 0/0.5/1 via bcftools %GT).

Method:
  1) Compute pairwise Manhattan SNP distances between samples:
       - Missing values are allowed.
       - For each pair (i,j), only loci where BOTH are non-missing are used.
       - The sum of |x_k - y_k| is scaled up proportionally to the total
         number of loci (matching the behaviour of R's amap::Dist). 
  2) Run classical MDS / PCoA (cmdscale-style) on the distance matrix.
  3) Plot requested PC pairs (default PC1–PC2 and PC1–PC3) colored by metadata
     group columns (e.g. region, country, year).

No imputation of missing genotypes is performed.
"""

from __future__ import annotations
import os
import sys
import subprocess
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------
# Utilities
# ------------------------
def _require_tool(name: str):
    if subprocess.run(["bash", "-lc", f"command -v {name} >/dev/null 2>&1"]).returncode != 0:
        sys.exit(f"ERROR: required tool '{name}' not found in PATH.")

def _run_cmd(cmd: List[str]) -> str | None:
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(f"[ERROR] {' '.join(cmd)}\n{p.stderr}")
        return None
    return p.stdout

def _vcf_samples(vcf_path: str) -> List[str]:
    out = _run_cmd(["bcftools", "query", "-l", vcf_path])
    if out is None:
        return []
    return [s.strip() for s in out.splitlines() if s.strip()]

def _intersect_ordered(a: List[str], b: List[str]) -> List[str]:
    bset = set(b)
    return [x for x in a if x in bset]

def _encode_gt_to_numeric(gt: str) -> float | np.nan:
    # '.', './.', '.|.' => NaN
    if gt is None or gt == "" or gt == "." or gt in (".|.", "./."):
        return np.nan
    al = gt.replace("|", "/").split("/")
    al = [x for x in al if x in ("0", "1")]
    if not al:
        return np.nan
    alt_count = sum(1 for x in al if x == "1")
    return alt_count / float(len(al))


# ------------------------
# Loaders
# ------------------------
def _load_matrix(matrix_path: str, metadata: pd.DataFrame, sample_col: str):
    """
    Expect a TSV with at least columns:
      - 'chr','pos','ref' (optional; will be dropped if present)
      - then per-sample columns
    Values: 0, 0.5, 1 or 'N' for missing.

    This version reads directly into float32 to reduce memory usage.
    """
    # First pass: just get columns
    header_df = pd.read_csv(matrix_path, sep="\t", nrows=0)
    cols = header_df.columns.tolist()

    start_idx = 0
    lower = [c.lower() for c in cols[:3]]
    if len(cols) >= 3 and lower == ["chr", "pos", "ref"]:
        start_idx = 3
    all_sample_cols = cols[start_idx:]

    meta_samples = metadata[sample_col].astype(str).tolist()
    sample_cols = _intersect_ordered(meta_samples, all_sample_cols)
    if not sample_cols:
        raise SystemExit("No overlapping samples between matrix and metadata.")

    # Build dtype dict: non-sample columns as 'category' or 'string',
    # sample columns as 'float32'. We read everything in one pass.
    dtype_map = {c: "string" for c in cols[:start_idx]}
    for c in sample_cols:
        dtype_map[c] = "float32"

    # Read full matrix; only columns we care about
    df = pd.read_csv(
        matrix_path,
        sep="\t",
        usecols=cols[:start_idx] + sample_cols,
        dtype=dtype_map,
        na_values=["N"],
    )

    # Extract sample matrix as float32, variants x samples
    M = df[sample_cols].to_numpy(dtype="float32")

    # Transpose to samples x variants
    X = M.T

    # Metadata subset in same sample order
    meta_sub = metadata.set_index(sample_col).loc[sample_cols].reset_index()

    return X, sample_cols, meta_sub



def _load_vcf_as_matrix(vcf_path: str, metadata: pd.DataFrame, sample_col: str):
    """
    Convert VCF to 0/0.5/1 matrix using bcftools query, no imputation.
    """
    _require_tool("bcftools")
    vcf_samples = _vcf_samples(vcf_path)
    if not vcf_samples:
        raise SystemExit("Could not read samples from VCF via bcftools.")
    meta_samples = metadata[sample_col].astype(str).tolist()
    samples = _intersect_ordered(meta_samples, vcf_samples)
    if not samples:
        raise SystemExit("No overlapping samples between VCF and metadata.")

    # Get VCF sample order
    hdr = _run_cmd(["bcftools", "query", "-H", "-f", "%CHROM\t%POS[\t%SAMPLE]\n", vcf_path])
    if hdr is None:
        raise SystemExit("bcftools query failed reading header.")
    all_samples = hdr.strip().splitlines()[-1].split("\t")[2:]

    # Query genotypes
    fmt = "%CHROM\t%POS[\t%GT]\n"
    body = _run_cmd(["bcftools", "query", "-f", fmt, vcf_path])
    if body is None or body.strip() == "":
        raise SystemExit("No variant rows returned by bcftools query.")

    lines = body.strip().splitlines()
    nvar, nsamp = len(lines), len(all_samples)
    M = np.full((nvar, nsamp), np.nan, dtype=float)
    for i, line in enumerate(lines):
        toks = line.split("\t")
        gts = toks[2:]
        for j, gt in enumerate(gts):
            M[i, j] = _encode_gt_to_numeric(gt)

    idx = [all_samples.index(s) for s in samples]
    M = M[:, idx]
    X = M.T  # samples x variants
    meta_sub = metadata.set_index(sample_col).loc[samples].reset_index()
    return X, samples, meta_sub


# ------------------------
# Distance + PCoA
# ------------------------
def _manhattan_dist_matrix(X: np.ndarray) -> np.ndarray:
    """
    Compute amap::Dist-style Manhattan distances between rows of X,
    allowing missing values (NaN) and scaling up when some columns
    are missing for a given pair.

    For each pair (i,j):
      - mask = finite entries in both i and j
      - if any shared loci:
          d_raw = sum( |x_ik - x_jk| over mask )
          d_scaled = d_raw * (n_total / n_used)
        else:
          d_scaled = 0 (no information; treat as 0 here)
    """
    n_samples, n_loci = X.shape
    D = np.zeros((n_samples, n_samples), dtype=float)

    for i in range(n_samples):
        xi = X[i, :]
        fi = np.isfinite(xi)
        for j in range(i + 1, n_samples):
            xj = X[j, :]
            fj = np.isfinite(xj)
            mask = fi & fj
            n_used = int(mask.sum())
            if n_used == 0:
                d = 0.0
            else:
                diff = np.abs(xi[mask] - xj[mask])
                d_raw = float(np.nansum(diff))
                d = d_raw * (n_loci / n_used)
            D[i, j] = D[j, i] = d

    return D


def _pcoa(distance_matrix: np.ndarray, n_components: int):
    """
    Classical MDS / PCoA (cmdscale-style) on a symmetric distance matrix.
    Returns:
      coords  : (n_samples, n_components)
      evr     : explained variance ratio per axis
    """
    D = np.asarray(distance_matrix, dtype=float)
    n = D.shape[0]

    # Double-centering of squared distances
    D2 = D ** 2
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J

    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eigh(B)
    # Sort descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # Keep positive components only
    eigvals = np.where(eigvals > 0, eigvals, 0.0)
    total_var = eigvals.sum()
    if total_var > 0:
        evr = eigvals / total_var
    else:
        evr = np.zeros_like(eigvals)

    k = min(n_components, (eigvals > 0).sum() if (eigvals > 0).any() else n)
    coords = eigvecs[:, :k] * np.sqrt(eigvals[:k])
    return coords, evr[:k]


# ------------------------
# Plotting
# ------------------------
def _scatter(scores: np.ndarray, evr: np.ndarray, meta: pd.DataFrame,
             group_col: str, pc_i: int, pc_j: int, out_pdf: str):
    """
    Plot PC{pc_i} vs PC{pc_j} (principal coordinates) colored by group_col.
    """
    if group_col not in meta.columns:
        print(f"[WARN] Skipping '{group_col}' (not in metadata).")
        return

    i, j = pc_i - 1, pc_j - 1
    if scores.shape[1] <= max(i, j):
        print(f"[WARN] Skipping PC{pc_i}-PC{pc_j}: not enough components.")
        return

    labels = meta[group_col].astype(str).fillna("NA").tolist()
    uniq = list(dict.fromkeys(labels))
    cmap = plt.cm.tab20.colors
    color_map = {g: cmap[k % len(cmap)] for k, g in enumerate(uniq)}
    colors = [color_map[g] for g in labels]

    def _axis_label(pc_idx: int) -> str:
        if pc_idx < len(evr):
            return f"PC{pc_idx+1} ({evr[pc_idx]*100:.2f}% var)"
        return f"PC{pc_idx+1}"

    xl = _axis_label(i)
    yl = _axis_label(j)

    plt.figure(figsize=(8.5, 6.2))
    ax = plt.gca()
    ax.scatter(scores[:, i], scores[:, j],
               s=32, c=colors, edgecolor="black", linewidth=0.3, alpha=0.9)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)

    handles = [
        plt.Line2D(
            [0], [0],
            marker="o", color="w",
            markerfacecolor=color_map[g],
            markeredgecolor="black",
            markersize=7,
            label=g,
        )
        for g in uniq
    ]
    ax.legend(handles=handles, loc="best", fontsize=8, frameon=True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)) or ".", exist_ok=True)
    plt.savefig(out_pdf)
    plt.close()
    print(f"[OK] Saved: {out_pdf}")


# ------------------------
# CLI entry
# ------------------------
def _parse_pcs_list(pcs_list: Optional[List[str]]) -> List[Tuple[int, int]]:
    """
    Parse ['1,2','1,3'] -> [(1,2),(1,3)].
    Defaults to [(1,2),(1,3)] if None given.
    """
    if not pcs_list:
        return [(1, 2), (1, 3)]
    out: List[Tuple[int, int]] = []
    for item in pcs_list:
        parts = item.split(",")
        if len(parts) != 2:
            sys.exit(f"Invalid --pcs '{item}'. Use format like --pcs 1,2")
        try:
            a, b = int(parts[0]), int(parts[1])
        except ValueError:
            sys.exit(f"Invalid --pcs '{item}'. Components must be integers.")
        if a < 1 or b < 1:
            sys.exit("--pcs values must be >= 1")
        out.append((a, b))
    return out


def run(
    matrix: Optional[str],
    vcf: Optional[str],
    metadata_path: str,
    outdir: str = "pca_plots",
    sample_col: str = "sample_id",
    group_by: Optional[List[str]] = None,
    max_sample_missing: Optional[float] = None,
    pcs: Optional[List[str]] = None,
):
    """
    Either --matrix or --vcf must be provided (mutually exclusive).

    Steps:
      - Load genotype matrix (0/0.5/1, 'N' -> NaN), samples x variants
      - Optionally drop samples with too much missing data (--max-sample-missing)
      - Compute amap-like Manhattan distance matrix with NA handling & scaling
      - Run classical MDS/PCoA
      - Plot requested PC pairs per group_by column
    """
    if (matrix is None) == (vcf is None):
        raise SystemExit("Provide exactly one of --matrix or --vcf.")

    os.makedirs(outdir, exist_ok=True)
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if sample_col not in meta.columns:
        raise SystemExit(f"Metadata missing required sample column '{sample_col}'.")

    # Load genotype data
    if matrix:
        X, samples, meta_sub = _load_matrix(matrix, meta, sample_col=sample_col)
    else:
        X, samples, meta_sub = _load_vcf_as_matrix(vcf, meta, sample_col=sample_col)

    # Optional sample missingness filter
    if max_sample_missing is not None:
        miss_prop = np.mean(np.isnan(X), axis=1)
        keep = miss_prop <= float(max_sample_missing)
        if keep.sum() < X.shape[0]:
            dropped = int((~keep).sum())
            print(f"[INFO] Dropping {dropped} samples for missingness > {max_sample_missing:.2f}")
            X = X[keep, :]
            meta_sub = meta_sub.iloc[keep].reset_index(drop=True)
            samples = [s for k, s in zip(keep, samples) if k]

    # Distance matrix (amap-like Manhattan)
    print("[INFO] Computing Manhattan distance matrix (amap-style, with NA handling)...")
    D = _manhattan_dist_matrix(X)

    # Decide PC pairs and number of components
    pc_pairs = _parse_pcs_list(pcs)
    max_pc = max(max(a, b) for a, b in pc_pairs)
    n_comp = max(max_pc, 3)

    # PCoA / cmdscale
    print("[INFO] Running classical MDS / PCoA...")
    coords, evr = _pcoa(D, n_components=n_comp)

    # Decide grouping columns
    if not group_by:
        candidates = ["region", "country", "year"]
        group_by = [c for c in candidates if c in meta_sub.columns]
        if not group_by:
            print("[WARN] No default group columns (region/country/year) found; will plot 'all' only.")
            group_by = []

    # Plot
    if group_by:
        for gcol in group_by:
            for (a, b) in pc_pairs:
                out_pdf = os.path.join(outdir, f"pca_{gcol}_PC{a}_PC{b}.pdf")
                _scatter(coords, evr, meta_sub, gcol, a, b, out_pdf)
    else:
        meta_all = meta_sub.copy()
        meta_all["all"] = "all"
        for (a, b) in pc_pairs:
            out_pdf = os.path.join(outdir, f"pca_all_PC{a}_PC{b}.pdf")
            _scatter(coords, evr, meta_all, "all", a, b, out_pdf)
