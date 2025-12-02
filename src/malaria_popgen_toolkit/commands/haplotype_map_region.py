# SPDX-License-Identifier: Apache-2.0
"""
Haplotype maps for drug-resistance genes by region (Africa/S. America/SE Asia).
- Reads multi-sample VCF + metadata TSV
- Intersects metadata sample IDs with VCF samples (avoids bcftools failure)
- Applies per-sample DP filtering (DP >= min_dp) when deciding haplotypes
- Assembles per-gene haplotypes from fixed codon sets
- Exports per-gene haplotype frequency CSVs
- Plots regional pie charts with GeoPandas (geodatasets or Natural Earth fallback)

Note: If a codon site is absent in the VCF (invariant), treat as REF.
"""

from __future__ import annotations
import os
import re
import subprocess
from collections import Counter
from difflib import get_close_matches
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.patches import Wedge
import matplotlib.patches as mpatches

# Prefer geodatasets if available; otherwise fetch Natural Earth zip directly.
try:
    from geodatasets import get_path as _gd_get_path
    HAS_GEODATASETS = True
except Exception:
    HAS_GEODATASETS = False


# ---------------------------------------------------------------------
# Gene definitions (Pf3D7 v3). POS -> (codon_no, REF_AA, MUT_AA)
# ---------------------------------------------------------------------
GENE_CODON_POSITIONS: Dict[str, Dict[str, Dict[int, Tuple[int, str, str]]]] = {
    "CRT": {
        "chrom": "Pf3D7_07_v3",
        "positions": {
            403596: (72, 'C', 'R'),
            403599: (73, 'V', 'F'),
            403602: (74, 'M', 'I'),
            403618: (74, 'MNK', 'IET'),
            403605: (75, 'N', 'E'),
            403625: (76, 'K', 'T'),
        },
    },
    "MDR1": {
        "chrom": "Pf3D7_05_v3",
        "positions": {
            958145: (86, 'N', 'Y'),
            958440: (184, 'Y', 'F'),
            961462: (1034, 'S', 'C'),
            961494: (1042, 'N', 'S'),
            961625: (1246, 'D', 'Y'),
        },
    },
    "DHFR": {
        "chrom": "Pf3D7_04_v3",
        "positions": {
            748239: (51,  'N', 'I'),
            748262: (59,  'C', 'R'),
            748410: (108, 'S', 'N'),
            748577: (164, 'I', 'L'),
        },
    },
    "DHPS": {
        "chrom": "Pf3D7_08_v3",
        "positions": {
            549681: (436, 'S', 'H'),
            549685: (437, 'G', 'A'),
            549993: (540, 'K', 'E'),
            550117: (581, 'A', 'G'),
        },
    },
}

WILDTYPE_HAPLO = {
    "DHFR": "NCSI",
    "DHPS": "SGKA",
    "CRT":  "CVMNK",
    "MDR1": "NYSND",
}


# ---------------------------------------------------------------------
# Country-name canonicalization (aliases -> Natural Earth names)
# ---------------------------------------------------------------------
COUNTRY_ALIASES = {
    # Explicit
    "drc": "Democratic Republic of the Congo",
    "gambia": "The Gambia",
    # Common variants
    "democratic republic of congo": "Democratic Republic of the Congo",
    "congo (kinshasa)": "Democratic Republic of the Congo",
    "congo, the democratic republic of the": "Democratic Republic of the Congo",
    "the gambia": "The Gambia",
    # Extras
    "ivory coast": "Côte d’Ivoire",
    "cote d'ivoire": "Côte d’Ivoire",
    "swaziland": "Eswatini",
    "cape verde": "Cabo Verde",
}

def _norm(s: str) -> str:
    s = s.lower().strip()
    s = re.sub(r"\s+", " ", s)
    return s

def _canonical_country_name(user_name: str, shape_names: list[str]) -> str | None:
    """Return shapefile country name that best matches user_name."""
    if user_name in shape_names:
        return user_name
    alias = COUNTRY_ALIASES.get(_norm(user_name))
    if alias and alias in shape_names:
        return alias
    norm_map = {_norm(n): n for n in shape_names}
    if _norm(user_name) in norm_map:
        return norm_map[_norm(user_name)]
    cand = get_close_matches(user_name, shape_names, n=1, cutoff=0.88)
    return cand[0] if cand else None


# ---------------------------------------------------------------------
# bcftools helpers
# ---------------------------------------------------------------------
def _run_cmd(cmd: List[str]) -> str | None:
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(f"[ERROR] {' '.join(cmd)}\n{p.stderr}")
        return None
    return p.stdout

def _vcf_samples(vcf_path: str) -> set:
    out = _run_cmd(['bcftools', 'query', '-l', vcf_path])
    if out is None:
        return set()
    return {s.strip() for s in out.splitlines() if s.strip()}

def _query_gene_region_with_gt_dp(vcf_path: str, chrom: str, start: int, end: int) -> str | None:
    fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP]\n'
    return _run_cmd(['bcftools', 'query', '-f', fmt, '-r', f'{chrom}:{start}-{end}', vcf_path])

def _sample_order_for_query(vcf_path: str, chrom: str, start: int, end: int) -> List[str]:
    hdr = _run_cmd([
        'bcftools', 'query', '-H',
        '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE]\n',
        '-r', f'{chrom}:{start}-{end}', vcf_path
    ])
    if hdr:
        cols = hdr.strip().split('\n')[-1].split('\t')
        return cols[4:]  # after CHROM,POS,REF,ALT
    return sorted(_vcf_samples(vcf_path))


# ---------------------------------------------------------------------
# Geo loader (robust across versions)
# ---------------------------------------------------------------------
def _load_world_countries(cache_dir: str | None = None) -> gpd.GeoDataFrame:
    if HAS_GEODATASETS:
        for key in [
            "naturalearth.cultural.admin_0_countries",
            "naturalearth.cultural.110m_admin_0_countries",
            "naturalearth.cultural.ne_110m_admin_0_countries",
            "naturalearth.cultural.v5_0.110m_admin_0_countries",
        ]:
            try:
                path = _gd_get_path(key)
                gdf = gpd.read_file(path).copy()
                if not gdf.empty:
                    return gdf
            except Exception:
                pass

    # Fallback: download Natural Earth 110m admin_0
    import urllib.request
    url = "https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip"
    if cache_dir is None:
        cache_dir = os.getcwd()
    os.makedirs(cache_dir, exist_ok=True)
    local_zip = os.path.join(cache_dir, "ne_110m_admin_0_countries.zip")
    if not os.path.exists(local_zip):
        print(f"[INFO] Downloading Natural Earth to {local_zip} ...")
        urllib.request.urlretrieve(url, local_zip)
    return gpd.read_file(f"zip://{local_zip}").copy()

def _region_geodf(region: str, cache_dir: str | None = None) -> gpd.GeoDataFrame:
    world = _load_world_countries(cache_dir)
    cols = {c.lower(): c for c in world.columns}
    cont_col = cols.get("continent") or cols.get("CONTINENT".lower())
    if cont_col is None:
        raise RuntimeError("No continent column in country layer.")
    world = world.to_crs(4326)
    cont = world[cont_col].astype(str).str.lower()
    region = region.lower()
    if region == "africa":
        return world[cont == "africa"].copy()
    elif region in ("south_america", "samerica"):
        return world[cont == "south america"].copy()
    elif region in ("southeast_asia", "seasia"):
        return world[cont == "asia"].copy()  # crude; refine with a country list if desired
    return world

def _country_name_column(gdf: gpd.GeoDataFrame) -> str:
    for cand in ("NAME_EN", "ADMIN", "NAME", "name_en", "admin", "name"):
        if cand in gdf.columns:
            return cand
    raise RuntimeError("Could not find a country name column in shapes.")

def _country_centroids(gdf: gpd.GeoDataFrame) -> Dict[str, Tuple[float, float]]:
    name_col = _country_name_column(gdf)
    pts = {}
    for _, row in gdf.iterrows():
        if row.geometry is None or row.geometry.is_empty:
            continue
        c = row.geometry.representative_point()
        pts[str(row[name_col])] = (float(c.x), float(c.y))
    return pts


# ---------------------------------------------------------------------
# Haplotype building
# ---------------------------------------------------------------------
def _parse_gt_dp_token(tok: str, min_dp: int) -> Tuple[int, int] | None:
    if ':' not in tok:
        return None
    gt_str, dp_str = tok.split(':', 1)
    if dp_str in ('.', ''):
        return None
    try:
        dp = int(dp_str)
    except ValueError:
        return None
    if dp < min_dp:
        return None
    if gt_str in ('.', './.', '.|.'):
        return None
    alleles = [a for a in gt_str.replace('|', '/').split('/') if a in ('0', '1')]
    if not alleles:
        return None
    alt_ct = sum(1 for a in alleles if a == '1')
    return (alt_ct, len(alleles))

def _haplotypes_for_gene(
    vcf_path: str,
    samples: List[str],
    gene: str,
    chrom: str,
    pos_map: Dict[int, Tuple[int, str, str]],
    min_dp: int
) -> Dict[str, str]:
    start, end = min(pos_map.keys()), max(pos_map.keys())
    raw = _query_gene_region_with_gt_dp(vcf_path, chrom, start, end)
    if raw is None or not raw.strip():
        return {}

    sample_order = _sample_order_for_query(vcf_path, chrom, start, end)
    sidx = {s: i for i, s in enumerate(sample_order)}
    per_sample_tokens: Dict[str, Dict[int, str]] = {s: {} for s in samples}

    for line in raw.strip().splitlines():
        f = line.split('\t')
        if len(f) < 5:
            continue
        pos = int(f[1])
        if pos not in pos_map:
            continue
        toks = f[4:]
        for s in samples:
            i = sidx.get(s)
            if i is None or i >= len(toks):
                continue
            per_sample_tokens[s][pos] = toks[i]

    ordered_positions = sorted(pos_map.keys(), key=lambda p: pos_map[p][0])
    hap_by_sample: Dict[str, str] = {}
    for s in samples:
        aaseq = []
        for p in ordered_positions:
            codon_no, refAA, mutAA = pos_map[p]
            tok = per_sample_tokens[s].get(p)
            if tok is None:
                aaseq.append(refAA)  # assume invariant/REF if absent
                continue
            parsed = _parse_gt_dp_token(tok, min_dp)
            if parsed is None:
                aaseq.append(refAA)  # fails DP or missing -> treat as REF
                continue
            alt_ct, _ = parsed
            aaseq.append(mutAA if alt_ct > 0 else refAA)
        hap_by_sample[s] = ''.join(aaseq)

    return hap_by_sample


# ---------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------
def _draw_pie(ax, center_xy, frac_list, color_list, radius_deg):
    x, y = center_xy
    start_ang = 90.0
    for f, c in zip(frac_list, color_list):
        if f <= 0:
            continue
        theta = 360.0 * f
        w = Wedge(center=(x, y), r=radius_deg,
                  theta1=start_ang, theta2=start_ang + theta,
                  facecolor=c, edgecolor="black", linewidth=0.5)
        ax.add_patch(w)
        start_ang += theta


# ---------------------------------------------------------------------
# Public entrypoint
# ---------------------------------------------------------------------
def run(region: str,
        vcf: str,
        metadata_path: str,
        outdir: str,
        min_dp: int = 5,
        sample_col: str = "sample_id",
        country_col: str = "country"):
    """
    Build and plot haplotype pies for a region; export per-gene haplotype frequencies.
    """
    os.makedirs(outdir, exist_ok=True)

    # Load metadata and sanitize sample IDs
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if sample_col not in meta.columns or country_col not in meta.columns:
        raise SystemExit(f"Metadata must contain '{sample_col}' and '{country_col}'.")
    meta[sample_col] = meta[sample_col].astype(str).str.strip()
    meta = meta[meta[sample_col].ne('') & meta[sample_col].notna()]

    # Intersect with VCF samples
    vcf_samps = _vcf_samples(vcf)
    if not vcf_samps:
        raise SystemExit("Could not read samples from VCF via 'bcftools query -l'.")
    groups: Dict[str, List[str]] = {}
    for country, dfc in meta.groupby(country_col):
        wanted = dfc[sample_col].tolist()
        present = [s for s in wanted if s in vcf_samps]
        missing = sorted(set(wanted) - set(present))
        if missing:
            print(f"[WARN] {country}: {len(missing)} metadata samples not in VCF. "
                  f"Examples: {', '.join(missing[:10])}")
        if present:
            groups[country] = present
    if not groups:
        raise SystemExit("No usable groups: none of the metadata samples matched the VCF.")

    # Build haplotypes per gene and aggregate per country (using user country names for now)
    all_hap_freqs: Dict[str, Dict[str, Counter]] = {}      # gene -> user_country -> Counter(hap)
    all_records: List[Dict[str, str | int | float]] = []   # for CSV export

    for gene, ginfo in GENE_CODON_POSITIONS.items():
        chrom = ginfo["chrom"]
        pos_map = ginfo["positions"]
        gene_country_counts: Dict[str, Counter] = {}

        for country, present_samples in groups.items():
            hap_by_sample = _haplotypes_for_gene(
                vcf_path=vcf,
                samples=present_samples,
                gene=gene,
                chrom=chrom,
                pos_map=pos_map,
                min_dp=min_dp
            )
            if not hap_by_sample:
                continue

            counts = Counter(hap_by_sample.values())
            total = sum(counts.values())
            gene_country_counts[country] = counts

            for hap, c in counts.items():
                all_records.append({
                    "Region": region,
                    "Country": country,
                    "Gene": gene,
                    "Haplotype": hap,
                    "Frequency": round(c / total, 6) if total > 0 else 0.0,
                    "SampleSize": total,
                    "Is_Wildtype": hap == WILDTYPE_HAPLO.get(gene, "")
                })

        all_hap_freqs[gene] = gene_country_counts

        # Export per-gene CSV (keeps user country names as provided)
        gene_out = os.path.join(outdir, f"{gene}_haplotype_frequencies_{region}.csv")
        gene_recs = [r for r in all_records if r["Gene"] == gene]
        if gene_recs:
            pd.DataFrame(gene_recs).to_csv(gene_out, index=False)
            print(f"[OK] Saved CSV: {gene_out}")
        else:
            print(f"[WARN] No data for {gene}")

    # --- Remap country keys to shapefile names before plotting ---
    cache_dir = os.path.abspath(outdir)
    region_gdf = _region_geodf(region, cache_dir=cache_dir)
    name_col = _country_name_column(region_gdf)
    shape_names = region_gdf[name_col].astype(str).tolist()

    remapped: Dict[str, Dict[str, Counter]] = {}
    for gene, country_counts in all_hap_freqs.items():
        cc2: Dict[str, Counter] = {}
        for user_country, counter in country_counts.items():
            canon = _canonical_country_name(user_country, shape_names)
            if canon is None:
                print(f"[WARN] No map match for '{user_country}' – skipping pies for this name.")
                continue
            cc2.setdefault(canon, Counter()).update(counter)  # merge if multiple map to same
        remapped[gene] = cc2

    out_png = os.path.join(outdir, f"haplotype_map_{region}.png")
    _plot_region_pies(remapped, region_gdf, out_png)


# ---------------------------------------------------------------------
# Plotting (2x2 panels)
# ---------------------------------------------------------------------
def _plot_region_pies(hap_dict: Dict[str, Dict[str, Counter]],
                      region_gdf: gpd.GeoDataFrame,
                      out_png: str):
    name_col = _country_name_column(region_gdf)
    centroids = _country_centroids(region_gdf)

    gene_order = [g for g in ["CRT", "MDR1", "DHFR", "DHPS"] if g in hap_dict]
    if not gene_order:
        print("[WARN] No data to plot.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    axes = axes.flatten()

    # Lock map extent & aspect, draw basemap
    xmin, ymin, xmax, ymax = region_gdf.total_bounds
    padx, pady = (xmax - xmin) * 0.04, (ymax - ymin) * 0.04
    for ax in axes:
        ax.axis("off")
        region_gdf.plot(ax=ax, color="white", edgecolor="0.6", linewidth=0.8)
        ax.set_xlim(xmin - padx, xmax + padx)
        ax.set_ylim(ymin - pady, ymax + pady)
        ax.set_aspect("equal", adjustable="box")

    # Panel letters
    panel_labels = ["A", "B", "C", "D"]

    for i, gene in enumerate(gene_order):
        ax = axes[i]
        ax.text(0.01, 0.98, panel_labels[i], transform=ax.transAxes,
                va="top", ha="left", fontsize=16, fontweight="bold")
        ax.set_title(gene, fontsize=14, pad=6)

        country_counts = hap_dict[gene]
        all_haps = sorted({h for c in country_counts.values() for h in c.keys()})
        colors = plt.cm.tab20.colors
        hap_colors = {h: colors[j % len(colors)] for j, h in enumerate(all_haps)}

        # Legend with codon numbers
        codon_numbers = [str(v[0]) for v in GENE_CODON_POSITIONS[gene]["positions"].values()]
        legend_handles = [
            mpatches.Patch(color=hap_colors[h],
                           label=f"{h}{' (WT)' if WILDTYPE_HAPLO.get(gene) == h else ''}")
            for h in all_haps
        ]
        ax.legend(
            handles=legend_handles,
            loc="lower left",
            fontsize=8,
            title=f"Codons: {', '.join(codon_numbers)}",
            title_fontsize=9,
            frameon=True,
        )

        # Pies
        map_w = xmax - xmin
        radius_deg = max(0.6, 0.036 * map_w)  # slightly bigger pies; adjust factor to taste
        for _, row in region_gdf.iterrows():
            country = str(row[name_col])
            if country not in country_counts:
                continue
            lon, lat = centroids.get(country, (None, None))
            if lon is None:
                continue
            cnt = country_counts[country]
            total = sum(cnt.values())
            if total == 0:
                continue
            fracs = [cnt.get(h, 0) / total for h in all_haps]
            _draw_pie(ax, (lon, lat), fracs, [hap_colors[h] for h in all_haps], radius_deg)

            # --- Wrapped country labels (2 lines if >2 words) + optional extra offset ---
            words = country.split()
            if len(words) > 2:
                mid = 2 if len(words) > 4 else len(words) // 2
                wrapped_name = "\n".join([" ".join(words[:mid]), " ".join(words[mid:])])
            else:
                wrapped_name = country

            ax.text(
                lon, lat + radius_deg * 1.5,  # optional extra offset applied
                wrapped_name,
                ha="center", va="bottom",
                fontsize=8, fontweight="bold",
                linespacing=1.1,
            )

    # Hide any unused axes
    for j in range(len(gene_order), 4):
        axes[j].axis("off")

    plt.tight_layout(pad=2.0)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved map: {out_png}")



