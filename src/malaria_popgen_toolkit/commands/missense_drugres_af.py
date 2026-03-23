# SPDX-License-Identifier: Apache-2.0
"""
Compute allele frequencies for missense variants in drug-resistance genes,
filtering out sample genotypes with depth DP < min_dp (default 5).

Supports hierarchical grouping by one or more metadata columns, e.g.
--group-by country
--group-by country year

For multiple grouping columns, outputs are written for:
  - level 1: country
  - level 2: country/year
"""

import subprocess
import pandas as pd
import os
import re
from typing import List

genes = {
    "CRT":  "Pf3D7_07_v3:403000-406000",
    "K13":  "Pf3D7_13_v3:1724817-1728000",
    "MDR1": "Pf3D7_05_v3:957500-962000",
    "DHFR": "Pf3D7_04_v3:746000-751000",
    "DHPS": "Pf3D7_08_v3:548000-552000",
    "PX1":  "Pf3D7_07_v3:889800-899213",
    "UBP1": "Pf3D7_01_v3:188400-201400",
    "AP2MU": "Pf3D7_12_v3:716700-720700",
    "AAT1": "Pf3D7_06_v3:1213100-1217350",
}

gene_aliases = {
    "CRT":   ["CRT"],
    "K13":   ["K13"],
    "MDR1":  ["MDR1"],
    "DHFR":  ["DHFR", "DHFR-TS"],
    "DHPS":  ["DHPS", "PPPK-DHPS"],
    "PX1":   ["PX1", "PF3D7_0720700"],
    "UBP1":  ["UBP1"],
    "AP2MU": ["AP2-MU", "AP2MU"],
    "AAT1":  ["AAT1", "PF3D7_0629500"],
}

def run_cmd(cmd):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running {' '.join(cmd)}:\n{result.stderr}")
        return None
    return result.stdout

def annotate_vcf(vcf, ref_fasta, gff3, out_annotated_vcf):
    cmd = ['bcftools', 'csq', '-p', 'a', '-f', ref_fasta, '-g', gff3,
           '-o', out_annotated_vcf, '-O', 'z', vcf]
    if run_cmd(cmd) is None:
        return False
    return run_cmd(['bcftools', 'index', out_annotated_vcf]) is not None

def query_with_gt_dp(vcf, region):
    fmt = '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%DP]\n'
    cmd = ['bcftools', 'query', '-f', fmt, '-r', region, vcf]
    return run_cmd(cmd)

def get_vcf_samples(vcf):
    out = run_cmd(['bcftools', 'query', '-l', vcf])
    if out is None:
        return set()
    return {x.strip() for x in out.strip().splitlines() if x.strip()}

def parse_output_and_recalculate_af(output, gene_name, min_dp=5):
    rows = []
    if not output:
        return pd.DataFrame(rows)

    alias_set = {a.upper() for a in gene_aliases.get(gene_name, [gene_name])}

    for line in output.strip().split('\n'):
        fields = line.split('\t')
        if len(fields) < 7:
            continue
        chrom, pos, ref, alt, csq = fields[0], fields[1], fields[2], fields[3], fields[4]
        sample_tokens = fields[5:]

        alt_count = 0
        total_alleles = 0
        contrib_samples = 0

        for tok in sample_tokens:
            if ':' not in tok:
                continue
            gt_str, dp_str = tok.split(':', 1)

            if dp_str in ('.', ''):
                continue
            try:
                dp = int(dp_str)
            except ValueError:
                continue
            if dp < min_dp:
                continue

            if gt_str in ('.', './.', '.|.'):
                continue

            alleles = [a for a in gt_str.replace('|', '/').split('/') if a in ('0', '1')]
            if not alleles:
                continue

            contrib_samples += 1
            alt_count += sum(1 for a in alleles if a == '1')
            total_alleles += len(alleles)

        if total_alleles == 0:
            continue

        af = alt_count / total_alleles

        for eff in csq.split(','):
            parts = eff.split('|')
            if len(parts) < 6:
                continue
            consequence = parts[0].lower()
            gene_tok1 = parts[1].upper() if len(parts) > 1 else ""
            gene_tok2 = parts[2].upper() if len(parts) > 2 else ""
            aa_change = parts[5] if len(parts) > 5 else f"{ref}>{alt}"

            if "missense" in consequence and ({gene_tok1, gene_tok2} & alias_set):
                rows.append({
                    "gene": gene_name,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "aa_change": aa_change,
                    "AF": af,
                    "N": contrib_samples,
                })

    return pd.DataFrame(rows)

def _safe(name: str) -> str:
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')

def _norm_group_by(group_by):
    if isinstance(group_by, str):
        return [group_by]
    return list(group_by)

def _make_group_label(group_cols, group_vals):
    return " | ".join(f"{c}={v}" for c, v in zip(group_cols, group_vals))

def _make_group_dir(base_out, group_vals):
    path = base_out
    for val in group_vals:
        path = os.path.join(path, _safe(val))
    return path

def _build_groups(metadata, group_cols):
    """
    Returns dict:
      key   = tuple of grouping values
      value = list of sample_ids
    """
    grouped = metadata.groupby(group_cols, dropna=False)['sample_id'].apply(list)

    groups = {}
    for key, sample_ids in grouped.items():
        if not isinstance(key, tuple):
            key = (key,)
        groups[key] = sample_ids
    return groups

def _run_one_grouping_level(vcf, ref_fasta, gff3, metadata, outdir, vcf_samples,
                            min_dp, group_cols, write_combined=True):
    """
    Run one grouping level, e.g.:
      ['country']
      ['country', 'year']
    """
    level_name = "_".join(group_cols)
    base_out = os.path.join(outdir, level_name)
    os.makedirs(base_out, exist_ok=True)

    groups = _build_groups(metadata, group_cols)
    all_data = []

    for group_vals, sample_ids in groups.items():
        present = [s for s in sample_ids if s in vcf_samples]
        missing = sorted(set(sample_ids) - set(present))

        group_desc = _make_group_label(group_cols, group_vals)

        if missing:
            print(f"[WARN] {group_desc}: {len(missing)} metadata samples not found in VCF "
                  f"(showing up to 10): {', '.join(missing[:10])}")

        if not present:
            print(f"[WARN] {group_desc}: no matching samples in VCF; skipping group.")
            continue

        print(f"\nProcessing group: {group_desc} (n_present={len(present)}, n_missing={len(missing)})")

        group_dir = _make_group_dir(base_out, group_vals)
        os.makedirs(group_dir, exist_ok=True)

        group_frames = []

        for gene, region in genes.items():
            print(f"  Gene {gene} at region {region} ...")

            subset_vcf_file = os.path.join(group_dir, f"{gene}_subset.vcf.gz")
            annotated_vcf_file = os.path.join(group_dir, f"{gene}_subset_csq.vcf.gz")

            cmd = [
                'bcftools', 'view',
                '-r', region,
                '-s', ','.join(present),
                '-Oz', '-o', subset_vcf_file,
                vcf
            ]
            if run_cmd(cmd) is None or run_cmd(['bcftools', 'index', subset_vcf_file]) is None:
                continue

            if not annotate_vcf(subset_vcf_file, ref_fasta, gff3, annotated_vcf_file):
                continue

            output = query_with_gt_dp(annotated_vcf_file, region)
            if output is None:
                continue

            df = parse_output_and_recalculate_af(output, gene_name=gene, min_dp=min_dp)
            if df.empty:
                continue

            for col, val in zip(group_cols, group_vals):
                df[col] = str(val)

            df["grouping_level"] = level_name
            df["group_label"] = _make_group_label(group_cols, group_vals)

            group_frames.append(df)

        if group_frames:
            group_df = pd.concat(group_frames, ignore_index=True)
            out_csv = os.path.join(group_dir, "missense_AF_per_gene.csv")
            group_df.to_csv(out_csv, index=False)
            print(f"  Saved CSV for {group_desc} to {out_csv}")
            all_data.append(group_df)
        else:
            print(f"[INFO] {group_desc}: no variants retained after filtering.")

    if all_data and write_combined:
        combined = pd.concat(all_data, ignore_index=True)

        sort_cols = [c for c in group_cols if c in combined.columns] + ["gene", "pos"]
        combined = combined.sort_values(by=sort_cols)

        combined_out = os.path.join(base_out, f"missense_AF_per_gene_by_{level_name}.csv")
        combined.to_csv(combined_out, index=False)
        print(f"Saved combined CSV to {combined_out}")

def run(vcf, ref_fasta, gff3, metadata_path, outdir,
        min_dp=5, group_by=None, write_combined=True):
    """
    CLI-facing entrypoint.

    group_by can be:
      - "country"
      - ["country"]
      - ["country", "year"]

    If multiple columns are provided, hierarchical outputs are written for:
      level 1: first column
      level 2: first+second
      ...
    """
    os.makedirs(outdir, exist_ok=True)

    group_by = _norm_group_by(group_by or ["country"])

    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if 'sample_id' not in metadata.columns:
        raise SystemExit("Metadata must have a 'sample_id' column.")

    for col in group_by:
        if col not in metadata.columns:
            raise SystemExit(f"Metadata does not contain the requested grouping column '{col}'.")

    metadata = metadata.dropna(subset=['sample_id']).copy()

    vcf_samples = get_vcf_samples(vcf)
    if not vcf_samples:
        raise SystemExit("Failed to read sample names from VCF (bcftools query -l).")

    # hierarchical levels:
    # ['country']
    # ['country', 'year']
    for i in range(1, len(group_by) + 1):
        level_cols = group_by[:i]
        _run_one_grouping_level(
            vcf=vcf,
            ref_fasta=ref_fasta,
            gff3=gff3,
            metadata=metadata,
            outdir=outdir,
            vcf_samples=vcf_samples,
            min_dp=min_dp,
            group_cols=level_cols,
            write_combined=write_combined,
        )


