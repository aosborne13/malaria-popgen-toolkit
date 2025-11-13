# SPDX-License-Identifier: Apache-2.0
"""
Make jittered dot plots of Fws from a metadata TSV.

Features:
- one PDF per grouping column (e.g., region, country, year)
- jittered points (black edge, filled by group color)
- mean crossbar per group
- "N = ..." annotation above (at y ~ 1.02)
- numeric mean label to the right of the crossbar
- y-limits: 0..1.05
"""

from __future__ import annotations
import os
import math
from typing import List

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def _plot_one_group(df: pd.DataFrame, group_var: str, out_pdf: str,
                    width: float = 10.0, height: float = 6.0):
    """
    Draw one jittered dot plot of Fws by group_var to a PDF.
    Expects columns: 'fws' and group_var.
    """
    # Drop rows with missing fws or group
    d = df[[group_var, "fws"]].copy()
    d = d.dropna(subset=["fws", group_var])
    if d.empty:
        print(f"[WARN] No data to plot for {group_var}")
        return

    # Ensure Fws is numeric
    d["fws"] = pd.to_numeric(d["fws"], errors="coerce")
    d = d.dropna(subset=["fws"])
    if d.empty:
        print(f"[WARN] No numeric Fws values for {group_var}")
        return

    # Factor order is natural order of appearance in file
    categories = pd.Categorical(d[group_var], ordered=True)
    d["group_factor"] = categories

    # Summary: mean and N
    summary = (
        d.groupby("group_factor", observed=True)
         .agg(mean_fws=("fws", "mean"), N=("fws", "size"))
         .reset_index()
    )
    summary["group_x"] = summary.index.astype(float) + 1.0  # numeric positions 1..K

    # Build color map per group for fill
    groups = summary["group_factor"].astype(str).tolist()
    unique_groups = list(dict.fromkeys(groups))  # preserve order
    cmap = plt.cm.tab20.colors
    color_lookup = {g: cmap[i % len(cmap)] for i, g in enumerate(unique_groups)}

    # Prepare jittered x positions
    # Centered at 1..K with small uniform jitter
    def jitter(series_idx, width=0.25):
        return series_idx + np.random.uniform(-width, width, size=series_idx.size)

    # Numeric x per row in d
    d = d.merge(summary[["group_factor", "group_x"]], on="group_factor", how="left")
    np.random.seed(1)  # reproducible jitter
    d["xj"] = jitter(d["group_x"].values)

    # Figure
    plt.figure(figsize=(width, height))
    ax = plt.gca()

    # Scatter points: filled by group color with black edge
    facecolors = [color_lookup[str(g)] for g in d["group_factor"].astype(str)]
    ax.scatter(d["xj"], d["fws"], s=35, c=facecolors, edgecolor="black", linewidth=0.6, alpha=0.85, zorder=2)

    # Mean crossbars
    for _, row in summary.iterrows():
        x = row["group_x"]
        y = row["mean_fws"]
        ax.hlines(y, x - 0.25, x + 0.25, colors="0.25", linewidth=1.0, zorder=3)

    # N labels at top (~1.02)
    for _, row in summary.iterrows():
        x = row["group_x"]
        ax.text(x, 1.02, f"N = {row['N']}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    # Mean labels just to the right of crossbar
    for _, row in summary.iterrows():
        x = row["group_x"] + 0.28  # offset to the right
        y = row["mean_fws"]
        ax.text(x, y, f"{row['mean_fws']:.2f}", ha="left", va="center", fontsize=9, color="black")

    # Axes styling
    ax.set_ylim(0.0, 1.05)
    ax.set_xlim(0.5, len(summary) + 0.5)
    ax.set_ylabel("Fws")
    ax.set_xlabel(group_var)

    # X tick labels as group names
    ax.set_xticks(summary["group_x"].values)
    ax.set_xticklabels([str(g) for g in summary["group_factor"]], rotation=45, ha="right")

    # Minimal theme
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")

    # Tight layout and save
    os.makedirs(os.path.dirname(os.path.abspath(out_pdf)) or ".", exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_pdf)
    plt.close()
    print(f"[OK] Saved: {out_pdf}")


def run(metadata_path: str,
        outdir: str = "fws_plots",
        group_by: List[str] | None = None,
        width: float = 10.0,
        height: float = 6.0):
    """
    Entry for CLI.
    - metadata_path: TSV with at least column 'fws' and one grouping column
    - outdir: where PDFs are saved
    - group_by: list of metadata columns to group by (default: any of ['region','country','year'] present)
    """
    os.makedirs(outdir, exist_ok=True)

    # Read metadata
    df = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if "fws" not in df.columns:
        raise SystemExit("Metadata must contain a 'fws' column.")

    # Convert Fws to numeric once here
    df["fws"] = pd.to_numeric(df["fws"], errors="coerce")

    # Default group_by selection
    if not group_by:
        candidates = ["region", "country", "year"]
        group_by = [c for c in candidates if c in df.columns]
        if not group_by:
            raise SystemExit("No grouping columns found. Provide --group-by with at least one metadata column.")

    # Loop over requested grouping columns (ignore any that are missing)
    for gcol in group_by:
        if gcol not in df.columns:
            print(f"[WARN] Skipping '{gcol}' (column not in metadata).")
            continue
        out_pdf = os.path.join(outdir, f"fws_by_{gcol}_jittered_labels.pdf")
        _plot_one_group(df, gcol, out_pdf, width=width, height=height)
