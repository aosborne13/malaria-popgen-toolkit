# SPDX-License-Identifier: Apache-2.0
"""
Compute allele frequencies for missense variants in drug-resistance genes,
filtering out sample genotypes with depth DP < min_dp (default 5).
Groups by an arbitrary metadata column (e.g., country/region/site/year).
Writes per-group CSVs and one combined CSV.
"""

import subprocess
import pandas as pd
import os
import re

# Define gene regions to subset (update as needed)
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

# Define gene aliases to match in BCSQ field
# NOTE: Add PF3D7 IDs where needed; AAT1 commonly appears as PF3D7_0629500 in BCSQ.
gene_aliases = {
    "CRT":   ["CRT"],
    "K13":   ["K13"],
    "MDR1":  ["MDR1"],
    "DHFR":  ["DHFR", "DHFR-TS"],
    "DHPS":  ["DHPS", "PPPK-DHPS"],
    "PX1":   ["PX1"],
    "UBP1":  ["UBP1"],
    "AP2MU": ["AP2-MU", "AP2MU"],
    "AAT1":  ["AAT1", "PF3D7_0629500"],  # <-- important for your case
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

def parse_output_and_recalculate_af(output, gene_name, min_dp=5):
    """
    Recalculate AF using ONLY genotypes with DP >= min_dp.
    Accept BCSQ gene tokens in either field [1] or [2] (ID vs symbol), case-insensitive.
    """
    rows = []
    if not output:
        return pd.DataFrame(rows)

    # Normalize alias list to uppercase for robust matching
    alias_set = {a.upper() for a in gene_aliases.get(gene_name, [gene_name])}

    for line in output.strip().split('\n'):
        fields = line.split('\t')
        if len(fields) < 7:
            continue
        chrom, pos, ref, alt, csq = fields[0], fields[1], fields[2], fields[3], fields[4]
        sample_tokens = fields[5:]

        alt_count = 0
        total_alleles = 0

        for tok in sample_tokens:
            # expect GT:DP
            if ':' not in tok:
                continue
            gt_str, dp_str = tok.split(':', 1)

            # DP filter
            if dp_str in ('.', ''):
                continue
            try:
                dp = int(dp_str)
            except ValueError:
                continue
            if dp < min_dp:
                continue

            # skip missing GT
            if gt_str in ('.', './.', '.|.'):
                continue

            alleles = [a for a in gt_str.replace('|', '/').split('/') if a in ('0', '1')]
            if not alleles:
                continue

            alt_count += sum(1 for a in alleles if a == '1')
            total_alleles += len(alleles)

        if total_alleles == 0:
            continue

        af = alt_count / total_alleles

        # parse BCSQ; match missense + gene alias against either token
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
                    "AF": af
                })

    return pd.DataFrame(rows)

def _safe(name: str) -> str:
    """Filesystem-safe folder/file token."""
    return re.sub(r'[^A-Za-z0-9._-]+', '_', str(name)).strip('_')

def run(vcf, ref_fasta, gff3, metadata_path, outdir,
        min_dp=5, group_by="country", write_combined=True):
    """
    CLI-facing entrypoint.
    - group_by: metadata column to group samples by (e.g., country/region/site/year)
    """
    os.makedirs(outdir, exist_ok=True)

    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if 'sample_id' not in metadata.columns:
        raise SystemExit("Metadata must have a 'sample_id' column.")
    if group_by not in metadata.columns:
        raise SystemExit(f"Metadata does not contain the requested grouping column '{group_by}'.")

    groups = metadata.groupby(group_by)['sample_id'].apply(list).to_dict()
    all_data = []

    # put outputs under outdir/<group_by>/<group_value>/
    base_out = os.path.join(outdir, group_by)
    os.makedirs(base_out, exist_ok=True)

    for group_value, sample_ids in groups.items():
        print(f"\nProcessing group: {group_by} = {group_value} ...")
        group_dir = os.path.join(base_out, _safe(group_value))
        os.makedirs(group_dir, exist_ok=True)
        group_frames = []

        for gene, region in genes.items():
            print(f"  Gene {gene} at region {region} ...")

            subset_vcf_file = os.path.join(group_dir, f"{gene}_subset.vcf.gz")
            annotated_vcf_file = os.path.join(group_dir, f"{gene}_subset_csq.vcf.gz")

            cmd = [
                'bcftools', 'view',
                '-r', region,
                '-s', ','.join(sample_ids),
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

            # add group label columns
            df[group_by] = str(group_value)
            df["group_label"] = f"{group_value} (n={len(sample_ids)})"
            group_frames.append(df)

        if group_frames:
            group_df = pd.concat(group_frames, ignore_index=True)
            out_csv = os.path.join(group_dir, "missense_AF_per_gene.csv")
            group_df.to_csv(out_csv, index=False)
            print(f"  Saved CSV for {group_by}={group_value} to {out_csv}")
            all_data.append(group_df)

    if not all_data:
        print("No data generated (after DP filtering).")
        return

    if write_combined:
        combined = pd.concat(all_data, ignore_index=True).sort_values(
            by=[group_by, "gene", "pos"]
        )
        combined_out = os.path.join(base_out, f"missense_AF_per_gene_by_{group_by}.csv")
        combined.to_csv(combined_out, index=False)
        print(f"Saved combined CSV to {combined_out}")


