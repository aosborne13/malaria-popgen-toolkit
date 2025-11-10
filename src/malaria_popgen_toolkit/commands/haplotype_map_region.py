# SPDX-License-Identifier: Apache-2.0
"""
Haplotype maps for drug-resistance genes by region.
- Reads a multi-sample VCF and metadata TSV
- Intersects metadata sample IDs with VCF samples (avoids bcftools failure)
- Applies per-sample DP filtering at variant sites (DP >= min_dp)
- Assembles per-gene haplotypes from a fixed set of codon positions
- Exports per-gene haplotype frequency CSVs per region
- Plots regional pie charts with GeoPandas (via geodatasets or Natural Earth fallback)

NOTE:
If a codon position is absent from the VCF, we treat those samples as REF at
that position (VCF typically omits invariant sites). DP filtering therefore
applies when the site is variant in at least one sample. For strict DP on
reference at invariant sites, a gVCF or coverage track would be needed.
"""

import os
import subprocess
from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.patches as mpatches

# geodatasets is preferred; if not available or the key changes, we fall back to Natural Earth zip
try:
    from geodatasets import get_path as _gd_get_path
    HAS_GEODATASETS = True
except Exception:
    HAS_GEODATASETS = False

# ----------------------------
# Gene definitions (Pf3D7 v3)
# ----------------------------
# positions: genomic POS -> (codon_number, ref_aa, mut_aa)
GENE_CODON_POSITIONS: Dict[str, Dict[str, Dict[int, Tuple[int, str, str]]]] = {
    "CRT": {
        "chrom": "Pf3D7_07_v3",
        "positions": {
            403596: (72, 'C', 'R'),
            403599: (73, 'V', 'F'),
            403602: (74, 'M', 'I'),
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

# ----------------------------
# bcftools helpers
# ----------------------------

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
    """
    Query VCF for GT:DP per sample over a window; later we filter to exact positions.
    """
    fmt = '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP]\n'
    return _run_cmd(['bcftools', 'query', '-f', fmt, '-r', f'{chrom}:{start}-{end}', vcf_path])

# ----------------------------
# Geo loader (robust across versions)
# ----------------------------

def _load_world_countries(cache_dir: str | None = None) -> gpd.GeoDataFrame:
    """
    Try geodatasets (several common keys), else download Natural Earth 110m admin_0.
    """
    # Try geodatasets keys (varies across versions)
    if HAS_GEODATASETS:
        candidate_keys = [
            # common/older keys
            "naturalearth.cultural.admin_0_countries",
            "naturalearth.cultural.110m_admin_0_countries",
            "naturalearth.cultural.ne_110m_admin_0_countries",
            # versioned keys (examples)
            "naturalearth.cultural.v1_1.50m_admin_0_countries",
            "naturalearth.cultural.v4_1.110m_admin_0_countries",
            "naturalearth.cultural.v5_0.110m_admin_0_countries",
        ]
        for key in candidate_keys:
            try:
                path = _gd_get_path(key)
                world = gpd.read_file(path).copy()
                if not world.empty:
                    return world
            except Exception:
                continue

    # Fallback: fetch Natural Earth zip
    import urllib.request
    url = "https://naturalearth.s3.amazonaws.com/110m_cultural/ne_110m_admin_0_countries.zip"
    if cache_dir is None:
        cache_dir = os.getcwd()
    os.makedirs(cache_dir, exist_ok=True)
    local_zip = os.path.join(cache_dir, "ne_110m_admin_0_countries.zip")
    if not os.path.exists(local_zip):
        print(f"[INFO] Downloading Natural Earth admin_0 countries to {local_zip} ...")
        urllib.request.urlretrieve(url, local_zip)
    world = gpd.read_file(f"zip://{local_zip}").copy()
    return world

def _region_geodf(region: str, cache_dir: str | None = None) -> gpd.GeoDataFrame:
    """
    Load countries and filter by region.
    """
    world = _load_world_countries(cache_dir=cache_dir)
    # Normalize schema
    cols = {c.lower(): c for c in world.columns}
    cont_col = cols.get("continent", None) or cols.get("CONTINENT".lower(), None)
    name_col = cols.get("name_en", None) or cols.get("admin", None) or cols.get("name", None)
    if cont_col is None:
        raise RuntimeError("No continent column found in country layer.")
    if name_col is None:
        # best-effort: use the first string-like column
        for c in world.columns:
            if world[c].dtype == object:
                name_col = c
                break
    world = world.to_crs(4326)

    region = region.lower()
    cont_series = world[cont_col].astype(str).str.lower()
    if region == "africa":
        return world[cont_series == "africa"].copy()
    elif region in ("south_america", "samerica"):
        return world[cont_series == "south america"].copy()
    elif region in ("southeast_asia", "seasia"):
        # crude: entire Asia (optionally restrict to a list if desired later)
        return world[cont_series == "asia"].copy()
    else:
        return world

def _country_centroids(gdf: gpd.GeoDataFrame) -> Dict[str, Tuple[float, float]]:
    """
    Compute (lon, lat) representative points for placing pies, keyed by country name.
    """
    name_col = None
    for cand in ("NAME_EN", "ADMIN", "NAME", "name_en", "admin", "name"):
        if cand in gdf.columns:
            name_col = cand
            break
    if name_col is None:
        raise RuntimeError("Could not find a name column for countries.")

    centroids = {}
    for _, row in gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        c = geom.representative_point()  # better for multipolygons
        centroids[str(row[name_col])] = (float(c.x), float(c.y))
    return centroids

# ----------------------------
# Haplotype building
# ----------------------------

def _parse_gt_dp_token(tok: str, min_dp: int) -> Tuple[int, int] | None:
    """
    Parse 'GT:DP' -> (alt_alleles_in_GT, total_alleles) if DP passes; else None.
    """
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
    """
    Build per-sample haplotype strings for a single gene.
    """
    # Query once over window
    start, end = min(pos_map.keys()), max(pos_map.keys())
    raw = _query_gene_region_with_gt_dp(vcf_path, chrom, start, end)
    if raw is None or not raw.strip():
        return {}  # no variants in window

    # Establish sample order returned by bcftools query
    header_line = _run_cmd([
        'bcftools', 'query', '-H',
        '-f', '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE]\n',
        '-r', f'{chrom}:{start}-{end}', vcf_path
    ])
    if header_line:
        cols = header_line.strip().split('\n')[-1].split('\t')
        sample_order = cols[4:]  # after CHROM,POS,REF,ALT
    else:
        sample_order = sorted(_vcf_samples(vcf_path))

    sample_idx = {s: i for i, s in enumerate(sample_order)}

    # Collect tokens per sample at our codon positions
    per_sample_tokens: Dict[str, Dict[int, str]] = {s: {} for s in samples}

    for line in raw.strip().splitlines():
        f = line.split('\t')
        if len(f) < 5:
            continue
        pos = int(f[1])
        if pos not in pos_map:
            continue  # not one of our codons
        tokens = f[4:]
        for s in samples:
            idx = sample_idx.get(s)
            if idx is None or idx >= len(tokens):
                continue
            per_sample_tokens[s][pos] = tokens[idx]

    # Assemble hap strings
    hap_by_sample: Dict[str, str] = {}
    ordered_positions = sorted(pos_map.keys(), key=lambda p: pos_map[p][0])

    for s in samples:
        aaseq = []
        for p in ordered_positions:
            codon_no, refAA, mutAA = pos_map[p]
            tok = per_sample_tokens[s].get(p)
            if tok is None:
                # site invariant / not present -> assume REF
                aaseq.append(refAA)
                continue
            parsed = _parse_gt_dp_token(tok, min_dp)
            if parsed is None:
                # fails DP or missing -> treat as REF for hap display
                aaseq.append(refAA)
                continue
            alt_ct, _ = parsed
            aaseq.append(mutAA if alt_ct > 0 else refAA)
        hap_by_sample[s] = ''.join(aaseq)

    return hap_by_sample

# ----------------------------
# Plotting
# ----------------------------

def _plot_region_pies(hap_dict: Dict[str, Dict[str, Counter]], region: str, out_png: str):
    """
    hap_dict: dict[gene][country] -> Counter({hap: count})
    """
    # Cache Natural Earth zip next to output (works offline later)
    cache_dir = os.path.dirname(os.path.abspath(out_png)) or os.getcwd()
    region_gdf = _region_geodf(region, cache_dir=cache_dir)
    centroids = _country_centroids(region_gdf)

    # Prepare figure: one subplot per gene (up to 4 genes -> 2x2)
    n_genes = len(hap_dict)
    rows = 2 if n_genes > 2 else 1
    cols = 2 if n_genes > 1 else 1
    fig, axes = plt.subplots(rows, cols, figsize=(20, 12))
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    axes = axes.flatten()

    for ax in axes:
        ax.axis("off")
        region_gdf.plot(ax=ax, color='whitesmoke', edgecolor='gray')

    for i, (gene, country_counts) in enumerate(hap_dict.items()):
        ax = axes[i]
        # Build palette for unique haplotypes across all countries for this gene
        all_haps = sorted({h for c in country_counts.values() for h in c.keys()})
        colors = plt.cm.tab20.colors
        hap_colors = {h: colors[j % len(colors)] for j, h in enumerate(all_haps)}

        # legend entries
        legend = [
            mpatches.Patch(color=hap_colors[h],
                           label=f"{h}{' (WT)' if WILDTYPE_HAPLO.get(gene) == h else ''}")
            for h in all_haps
        ]

        for country, counter in country_counts.items():
            if country not in centroids:
                continue
            lon, lat = centroids[country]
            total = sum(counter.values())
            if total == 0:
                continue
            sizes = [counter.get(h, 0) / total for h in all_haps]
            ax.pie(
                sizes,
                center=(lon, lat),
                radius=2.0,
                colors=[hap_colors[h] for h in all_haps],
                wedgeprops=dict(edgecolor='black', linewidth=0.5)
            )
            ax.text(lon, lat + 3.0, country, ha='center', fontsize=10, fontweight='bold')

        ax.set_title(f"{gene}", fontsize=14)
        ax.legend(handles=legend, loc='lower left', fontsize=9)

    plt.tight_layout(pad=3)
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"[OK] Saved map: {out_png}")

# ----------------------------
# Public entrypoint
# ----------------------------

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
        raise SystemExit(f"Metadata must contain '{sample_col}' and '{country_col}' columns.")
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

    # Build haplotypes per gene and aggregate per country
    all_hap_freqs: Dict[str, Dict[str, Counter]] = {}      # gene -> country -> Counter(hap)
    all_records: List[Dict[str, str | int | float]] = []   # for CSV export

    for gene, ginfo in GENE_CODON_POSITIONS.items():
        chrom = ginfo["chrom"]
        pos_map = ginfo["positions"]
        gene_country_counts: Dict[str, Counter] = {}

        for country, present_samples in groups.items():
            # per-sample hap strings for this gene
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

            # records for CSV
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

        # Export per-gene CSV
        gene_out = os.path.join(outdir, f"{gene}_haplotype_frequencies_{region}.csv")
        gene_recs = [r for r in all_records if r["Gene"] == gene]
        if gene_recs:
            pd.DataFrame(gene_recs).to_csv(gene_out, index=False)
            print(f"[OK] Saved CSV: {gene_out}")
        else:
            print(f"[WARN] No data for {gene}")

    # Plot combined map (dataset cached alongside outputs for reproducibility)
    out_png = os.path.join(outdir, f"haplotype_map_{region}.png")
    _plot_region_pies(all_hap_freqs, region, out_png)


