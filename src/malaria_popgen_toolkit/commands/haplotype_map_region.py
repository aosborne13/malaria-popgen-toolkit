# SPDX-License-Identifier: Apache-2.0
"""
Regional haplotype maps for drug-resistance genes (VCF + DP filter).
Regions: 'africa', 'south_america', 'southeast_asia'

- Reads VCF with bcftools and pulls GT:DP for specific Pf3D7 positions.
- Applies per-sample min DP (default 5).
- Metadata TSV must have: sample_id (default), country (default).
- Outputs per-gene CSVs and a combined regional PNG.
"""

import os
import re
import subprocess
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Pf3D7 v3 codon reference per gene (positions in 1-based genome coords)
GENE_CODON_POSITIONS = {
    "CRT": {
        "chrom": "Pf3D7_07_v3",
        "positions": {
            403596: (72, 'C'), 403599: (73, 'V'), 403602: (74, 'M'),
            403605: (75, 'N'), 403625: (76, 'K')
        }
    },
    "MDR1": {
        "chrom": "Pf3D7_05_v3",
        "positions": {
            958145: (86, 'N'), 958440: (184, 'Y'), 961462: (1034, 'S'),
            961494: (1042, 'N'), 961625: (1246, 'D')
        }
    },
    "DHFR": {
        "chrom": "Pf3D7_04_v3",
        "positions": {
            748239: (51, 'N'), 748262: (59, 'C'), 748410: (108, 'S'), 748577: (164, 'I')
        }
    },
    "DHPS": {
        "chrom": "Pf3D7_08_v3",
        "positions": {
            549681: (436, 'S'), 549685: (437, 'G'), 549993: (540, 'K'), 550117: (581, 'A')
        }
    }
}

MUTANT_AA = {
    "CRT": {72: 'R', 73: 'F', 74: 'I', 75: 'E', 76: 'T'},
    "MDR1": {86: 'Y', 184: 'F', 1034: 'C', 1042: 'S', 1246: 'Y'},
    "DHFR": {51: 'I', 59: 'R', 108: 'N', 164: 'L'},
    "DHPS": {436: 'H', 437: 'A', 540: 'E', 581: 'G'}
}

WILDTYPE_HAPLO = {"DHFR": "NCSI", "DHPS": "SAKA", "CRT": "CVMNK", "MDR1": "NYSND"}

COUNTRY_ALIASES = {
    "DRC": "Democratic Republic of the Congo",
    "Congo (Kinshasa)": "Democratic Republic of the Congo",
    "Congo, Dem. Rep.": "Democratic Republic of the Congo",
    "Congo (Brazzaville)": "Republic of the Congo",
    "Congo, Rep.": "Republic of the Congo",
    "Ivory Coast": "Côte d’Ivoire",
    "Cote d'Ivoire": "Côte d’Ivoire",
    "Côte d'Ivoire": "Côte d’Ivoire",
    "Eswatini": "Swaziland",
    "The Gambia": "Gambia",
    "eSwatini": "Swaziland",
    "Laos": "Lao PDR",
    "Timor Leste": "Timor-Leste",
}
SE_ASIA_COUNTRIES = {
    "Brunei", "Cambodia", "Indonesia", "Lao PDR", "Laos", "Malaysia",
    "Myanmar", "Philippines", "Singapore", "Thailand", "Timor-Leste", "Vietnam"
}

def _safe(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')

def _normalize_country(name: str) -> str:
    return COUNTRY_ALIASES.get(name, name)

def _run(cmd: list[str]) -> str | None:
    print("Running:", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print("[ERROR]", r.stderr.strip())
        return None
    return r.stdout

def _query_gene_positions(vcf: str, samples: list[str], gene: str, min_dp: int) -> pd.DataFrame:
    """Query GT:DP for the gene's positions using bcftools; return tidy long DF."""
    info = GENE_CODON_POSITIONS[gene]
    chrom = info["chrom"]
    pos_list = sorted(info["positions"].keys())

    # Compose -r ranges as individual positions (chrom:pos-pos)
    regions = ",".join([f"{chrom}:{p}-{p}" for p in pos_list])

    # bcftools query: CHROM POS REF ALT  [  sample1(GT:DP)  sample2(GT:DP) ... ]
    fmt = "%CHROM\t%POS[\t%GT:%DP]\n"
    cmd = ["bcftools", "query", "-f", fmt, "-r", regions, "-s", ",".join(samples), vcf]
    out = _run(cmd)
    if out is None or out.strip() == "":
        return pd.DataFrame(columns=["chrom", "pos", "sample", "gt", "dp"])

    # Parse to long format
    rows = []
    for line in out.strip().splitlines():
        parts = line.split("\t")
        chrom_v, pos_v = parts[0], int(parts[1])
        sample_fields = parts[2:]
        for samp, tok in zip(samples, sample_fields):
            if ":" in tok:
                gt, dp = tok.split(":", 1)
            else:
                gt, dp = tok, "."
            # apply DP filter here by marking low-DP as missing
            try:
                dp_int = int(dp)
            except ValueError:
                dp_int = None
            if dp_int is not None and dp_int < min_dp:
                gt_use = "."
            else:
                gt_use = gt
            rows.append({"chrom": chrom_v, "pos": pos_v, "sample": samp, "gt": gt_use, "dp": dp_int})
    return pd.DataFrame(rows)

def _hap_table_from_vcf(df_long: pd.DataFrame, gene: str) -> pd.DataFrame:
    """
    Convert long GT table to per-sample AA calls per codon; drop samples with any missing codon.
    0/0 -> ref AA; genotypes containing '1' -> mutant AA; missing -> NaN.
    """
    info = GENE_CODON_POSITIONS[gene]
    pos2codon = info["positions"]  # pos -> (codon_num, ref_aa)

    # Map GT to AA
    df = df_long.copy()
    def gt_to_state(gt: str) -> float | None:
        if gt in (".", "./.", ".|."):
            return None
        # haploid/mixed/diploid encodings: treat any '1' as mutant
        if "1" in gt:
            return 1.0
        return 0.0

    df["state"] = df["gt"].map(gt_to_state)

    # Pivot to samples x positions (state)
    pivot = df.pivot_table(index="sample", columns="pos", values="state", aggfunc="first")

    # Ensure all positions exist
    for pos in pos2codon:
        if pos not in pivot.columns:
            pivot[pos] = np.nan
    pivot = pivot[sorted(pos2codon.keys())]

    # Build AA calls per codon number
    aa = pd.DataFrame(index=pivot.index)
    for pos in pivot.columns:
        codon_num, ref_aa = pos2codon[pos]
        mut_aa = MUTANT_AA.get(gene, {}).get(codon_num, ref_aa)
        col = []
        for v in pivot[pos]:
            if pd.isna(v):
                col.append(np.nan)
            elif v == 0.0:
                col.append(ref_aa)
            else:  # 1.0
                col.append(mut_aa)
        aa[codon_num] = col

    # drop samples with any missing codon
    aa = aa.dropna(axis=0, how="any")
    return aa

def _hap_string(row: pd.Series, gene: str) -> str:
    codon_order = [v[0] for v in GENE_CODON_POSITIONS[gene]["positions"].values()]
    return "".join(row.get(c, "") for c in codon_order)

def _haplotypes_by_country(hap_table: pd.DataFrame,
                           metadata: pd.DataFrame,
                           sample_col: str,
                           country_col: str,
                           gene: str) -> dict[str, list[str]]:
    meta = metadata[[sample_col, country_col]].dropna()
    meta[country_col] = meta[country_col].map(_normalize_country)
    s2c = dict(zip(meta[sample_col], meta[country_col]))

    haps: dict[str, list[str]] = {}
    for sample in hap_table.index:
        country = s2c.get(sample)
        if not country:
            continue
        hap = _hap_string(hap_table.loc[sample], gene)
        haps.setdefault(country, []).append(hap)
    return haps

def _export_haplotype_proportions(hap_by_country: dict, gene: str, out_dir: str):
    records = []
    for country, hap_list in hap_by_country.items():
        n = len(hap_list)
        if n == 0:
            continue
        counts = pd.Series(hap_list).value_counts(normalize=True)
        for hap, freq in counts.items():
            records.append({
                "Gene": gene,
                "Country": country,
                "Haplotype": hap,
                "Frequency": round(float(freq), 4),
                "SampleSize": n,
                "Is_Wildtype": hap == WILDTYPE_HAPLO.get(gene, "")
            })
    if not records:
        return
    df_out = pd.DataFrame(records)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{gene}_haplotype_frequencies.csv")
    df_out.to_csv(out_path, index=False)
    print(f"Saved CSV: {out_path}")

def _region_geodf(region: str) -> gpd.GeoDataFrame:
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres")).copy()
    if region == "africa":
        df = world[world["continent"] == "Africa"].copy()
    elif region == "south_america":
        df = world[world["continent"] == "South America"].copy()
    elif region == "southeast_asia":
        df = world[world["name"].isin(SE_ASIA_COUNTRIES)].copy()
    else:
        raise SystemExit(f"Unsupported region: {region}")
    df["name_norm"] = df["name"].apply(_normalize_country)
    return df

def _region_bounds(ax, region: str):
    if region == "africa":
        ax.set_xlim(-25, 55); ax.set_ylim(-35, 25)
    elif region == "south_america":
        ax.set_xlim(-85, -30); ax.set_ylim(-60, 15)
    elif region == "southeast_asia":
        ax.set_xlim(90, 150); ax.set_ylim(-15, 30)

def _plot_region_pies(all_haps: dict, region: str, out_png: str):
    region_df = _region_geodf(region)
    fig, axes = plt.subplots(2, 2, figsize=(22, 18))
    axes = axes.flatten()

    for ax, (gene, hap_by_country) in zip(axes, all_haps.items()):
        region_df.plot(ax=ax, color='whitesmoke', edgecolor='gray')

        # Unique hap strings for stable colors
        unique_haps = sorted(set().union(*[set(pd.Series(h).unique()) for h in hap_by_country.values()]) or [])
        cmap = plt.cm.tab20.colors
        hap_colors = {h: cmap[i % len(cmap)] for i, h in enumerate(unique_haps)}

        for country, hlist in hap_by_country.items():
            if not hlist:
                continue
            row = region_df[region_df["name_norm"] == _normalize_country(country)]
            if row.empty:
                continue
            pt = row.geometry.values[0].representative_point()
            x, y = pt.x, pt.y

            counts = pd.Series(hlist).value_counts(normalize=True)
            sizes = [counts.get(h, 0.0) for h in unique_haps]
            if sum(sizes) == 0:
                continue

            ax.pie(sizes, center=(x, y), radius=2.5,
                   colors=[hap_colors[h] for h in unique_haps],
                   wedgeprops=dict(edgecolor='black'))
            ax.text(x, y + 2.5, country, ha='center', fontsize=12, fontweight='bold')

        legend_handles = [
            mpatches.Patch(color=hap_colors[h],
                           label=f"{h}{' (WT)' if h == WILDTYPE_HAPLO.get(gene, '') else ''}")
            for h in unique_haps
        ]
        codons = [str(v[0]) for v in GENE_CODON_POSITIONS[gene]['positions'].values()]
        ax.legend(handles=legend_handles, loc='lower left', fontsize=11,
                  title=f"{gene} (codons: {', '.join(codons)})", title_fontsize=12)

        _region_bounds(ax, region)
        ax.axis("off")

    plt.tight_layout(pad=4)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f"Saved map: {out_png}")

def run(region: str,
        vcf: str,
        metadata_path: str,
        outdir: str,
        min_dp: int = 5,
        sample_col: str = "sample_id",
        country_col: str = "country"):
    """
    Region-aware CLI entrypoint (VCF + DP).
    region: 'africa' | 'south_america' | 'southeast_asia'
    """
    os.makedirs(outdir, exist_ok=True)
    plots_dir = os.path.join(outdir, "plots")
    csv_dir = os.path.join(outdir, "haplotype_freqs")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(csv_dir, exist_ok=True)

    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if sample_col not in meta.columns or country_col not in meta.columns:
        raise SystemExit(f"Metadata must contain '{sample_col}' and '{country_col}'.")
    samples = meta[sample_col].dropna().unique().tolist()

    all_haps = {}
    for gene in ["CRT", "MDR1", "DHFR", "DHPS"]:
        print(f"[Gene] {gene}")
        df_long = _query_gene_positions(vcf, samples=samples, gene=gene, min_dp=min_dp)
        if df_long.empty:
            print(f"[WARN] No data for {gene}")
            all_haps[gene] = {}
            continue
        hap_table = _hap_table_from_vcf(df_long, gene)
        hap_by_country = _haplotypes_by_country(
            hap_table, meta,
            sample_col=sample_col, country_col=country_col, gene=gene
        )
        all_haps[gene] = hap_by_country
        _export_haplotype_proportions(hap_by_country, gene, csv_dir)

    out_png = os.path.join(plots_dir, f"haplotype_map_{region}.png")
    _plot_region_pies(all_haps, region, out_png)
