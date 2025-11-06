# SPDX-License-Identifier: Apache-2.0
"""
Haplotype map for Africa (drug-resistance genes).
- Works with any African countries present in metadata (no hard-coded list).
- Expects a binary matrix with columns: chr, pos, ref, <sample1> <sample2> ...
  Values per sample: 0=ref, 0.5=mixed, 1=alt, N/NA=missing
- Expects metadata TSV with columns: sample_id, country (configurable)
- Outputs:
  - CSV per gene with haplotype frequencies by country
  - A combined PNG map with per-country pie charts
"""

import os
import re
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from shapely.geometry import Point

# Reference amino acids and expected chromosome per gene (Pf3D7 v3)
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

WILDTYPE_HAPLO = {
    "DHFR": "NCSI",
    "DHPS": "SAKA",
    "CRT": "CVMNK",
    "MDR1": "NYSND"
}

# Common country alias mapping for Natural Earth names (minimal set)
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
}

def _safe(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')

def _normalize_country(name: str) -> str:
    return COUNTRY_ALIASES.get(name, name)

def _subset_matrix_for_gene(matrix_path: str, gene_key: str) -> pd.DataFrame:
    """Return a DataFrame containing only rows for the gene's chrom+positions."""
    info = GENE_CODON_POSITIONS[gene_key]
    chrom = info["chrom"]
    positions = set(info["positions"].keys())

    # Read base to find indices for desired rows
    base = pd.read_csv(matrix_path, sep='\t', usecols=['chr', 'pos', 'ref'])
    idx = base[(base['chr'] == chrom) & (base['pos'].isin(positions))].index

    # Read only header + desired rows (skip all others)
    df = pd.read_csv(
        matrix_path, sep='\t',
        skiprows=lambda x: x != 0 and (x - 1) not in set(idx)
    )
    # Keep sanity
    df = df[df['chr'] == chrom]
    return df

def _haplotype_table(df_gene: pd.DataFrame, gene_key: str) -> pd.DataFrame:
    """
    Build a per-sample table of AA calls at specified codons and return a DataFrame:
    index = sample, columns = codon numbers (e.g., 72, 73, ...), values = AA letter
    """
    info = GENE_CODON_POSITIONS[gene_key]
    codons_map = info["positions"]  # pos -> (codon_number, ref_aa)

    # Set index to position; drop non-sample cols
    df = df_gene.set_index('pos').drop(columns=['chr', 'ref'], errors='ignore')
    sample_cols = [c for c in df.columns if c not in ('chr', 'ref')]

    # Ensure all codons exist as rows (fill missing with NaN)
    for pos in codons_map:
        if pos not in df.index:
            df.loc[pos, sample_cols] = np.nan

    # Reorder by the codon order
    order = [p for p in codons_map]
    df = df.loc[order]
    df = df.sort_index()

    # Transpose to samples x positions
    dfT = df.T  # rows=samples, cols=genome positions

    # Build AA calls per codon number
    out = pd.DataFrame(index=dfT.index)
    for pos, (codon_num, ref_aa) in codons_map.items():
        vals = dfT[pos]
        aa_calls = []
        mut_aa = MUTANT_AA.get(gene_key, {}).get(codon_num, ref_aa)
        for v in vals:
            try:
                if float(v) == 0:
                    aa_calls.append(ref_aa)
                elif float(v) in (0.5, 1.0):
                    aa_calls.append(mut_aa)
                else:
                    aa_calls.append(ref_aa)
            except Exception:
                aa_calls.append(ref_aa)
        out[codon_num] = aa_calls
    return out

def _hap_string(row: pd.Series, gene_key: str) -> str:
    # Preserve defined codon order for the gene
    codon_order = [v[0] for v in GENE_CODON_POSITIONS[gene_key]["positions"].values()]
    parts = []
    for c in codon_order:
        parts.append(row.get(c, ''))
    return ''.join(parts)

def _haplotypes_by_country(hap_table: pd.DataFrame,
                           metadata: pd.DataFrame,
                           sample_col: str,
                           country_col: str,
                           gene_key: str) -> dict:
    """
    Returns dict: {country_name: [hap1, hap2, ...]}
    Only countries present in metadata are included.
    """
    # Map sample -> country
    meta = metadata[[sample_col, country_col]].dropna()
    meta[country_col] = meta[country_col].map(_normalize_country)
    sample_to_country = dict(zip(meta[sample_col], meta[country_col]))

    haps = {}
    for sample in hap_table.index:
        country = sample_to_country.get(sample)
        if not country:
            continue
        hap = _hap_string(hap_table.loc[sample], gene_key)
        haps.setdefault(country, []).append(hap)
    return haps

def _export_haplotype_proportions(hap_by_country: dict, gene_key: str, out_dir: str):
    records = []
    for country, hap_list in hap_by_country.items():
        total = len(hap_list)
        if total == 0:
            continue
        counts = pd.Series(hap_list).value_counts(normalize=True)
        for hap, freq in counts.items():
            records.append({
                "Gene": gene_key,
                "Country": country,
                "Haplotype": hap,
                "Frequency": round(float(freq), 4),
                "SampleSize": total,
                "Is_Wildtype": hap == WILDTYPE_HAPLO.get(gene_key, "")
            })
    if not records:
        return None
    df_out = pd.DataFrame(records)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{gene_key}_haplotype_frequencies.csv")
    df_out.to_csv(out_path, index=False)
    print(f"Saved CSV: {out_path}")
    return df_out

def _plot_africa_pies(all_haps: dict, out_png: str):
    """
    Plot per-country pies for each gene on an Africa basemap.
    all_haps: {gene: {country: [hap, hap, ...]}}
    """
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    africa = world[world["continent"] == "Africa"].copy()
    africa["name_norm"] = africa["name"].apply(_normalize_country)

    fig, axes = plt.subplots(2, 2, figsize=(22, 18))
    axes = axes.flatten()

    for ax, (gene, hap_by_country) in zip(axes, all_haps.items()):
        africa.plot(ax=ax, color='whitesmoke', edgecolor='gray')

        # Build frequencies per country
        # First pass: collect all unique hap strings for a stable color map
        unique_haps = set()
        for hlist in hap_by_country.values():
            unique_haps |= set(pd.Series(hlist).unique().tolist())
        unique_haps = sorted(unique_haps)
        cmap = plt.cm.tab20.colors
        hap_colors = {h: cmap[i % len(cmap)] for i, h in enumerate(unique_haps)}

        summary_rows = []
        for country, hlist in hap_by_country.items():
            if not hlist:
                continue
            # Find country geometry
            name_norm = _normalize_country(country)
            row = africa[africa["name_norm"] == name_norm]
            if row.empty:
                print(f"[WARN] Country not found on map: {country}")
                continue
            # representative_point is safer than centroid for multipolygons
            pt = row.geometry.values[0].representative_point()
            x, y = pt.x, pt.y

            counts = pd.Series(hlist).value_counts(normalize=True)
            sizes = [counts.get(h, 0.0) for h in unique_haps]
            if sum(sizes) == 0:
                continue

            ax.pie(
                sizes,
                center=(x, y),
                radius=2.5,
                colors=[hap_colors[h] for h in unique_haps],
                wedgeprops=dict(edgecolor='black')
            )
            ax.text(x, y + 3.0, country, ha='center', fontsize=13, fontweight='bold')

            summary_rows.append({
                "Country": country,
                "Longitude": x,
                "Latitude": y,
                "SampleSize": len(hlist),
            })

        # Legend
        legend_handles = [
            mpatches.Patch(color=hap_colors[h],
                           label=f"{h}{' (WT)' if h == WILDTYPE_HAPLO.get(gene, '') else ''}")
            for h in unique_haps
        ]
        codon_numbers = [str(v[0]) for v in GENE_CODON_POSITIONS[gene]['positions'].values()]
        ax.legend(handles=legend_handles, loc='lower left', fontsize=11,
                  title=f"{gene}  (codons: {', '.join(codon_numbers)})", title_fontsize=12)

        ax.set_xlim(-25, 55)
        ax.set_ylim(-35, 25)
        ax.axis("off")

    plt.tight_layout(pad=4)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f"Saved map: {out_png}")

def run(matrix_path: str,
        metadata_path: str,
        outdir: str,
        sample_col: str = "sample_id",
        country_col: str = "country"):
    """
    CLI entrypoint for Africa haplotype map.
    """
    os.makedirs(outdir, exist_ok=True)
    plots_dir = os.path.join(outdir, "plots")
    csv_dir = os.path.join(outdir, "haplotype_freqs")
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(csv_dir, exist_ok=True)

    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if sample_col not in metadata.columns or country_col not in metadata.columns:
        raise SystemExit(f"Metadata must contain '{sample_col}' and '{country_col}'.")

    # Build haplotypes per gene
    all_haps = {}
    for gene in ["CRT", "MDR1", "DHFR", "DHPS"]:
        print(f"\n[Gene] {gene}")
        df_gene = _subset_matrix_for_gene(matrix_path, gene)
        hap_table = _haplotype_table(df_gene, gene)
        hap_by_country = _haplotypes_by_country(
            hap_table, metadata,
            sample_col=sample_col, country_col=country_col,
            gene_key=gene
        )
        all_haps[gene] = hap_by_country
        # CSV per gene
        _export_haplotype_proportions(hap_by_country, gene, csv_dir)

    # Combined map
    out_png = os.path.join(plots_dir, "haplotype_map_africa.png")
    _plot_africa_pies(all_haps, out_png)
