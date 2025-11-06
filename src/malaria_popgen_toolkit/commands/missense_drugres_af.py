# SPDX-License-Identifier: Apache-2.0
"""
Compute allele frequencies for missense variants in drug-resistance genes,
filtering out sample genotypes with depth DP < min_dp (default 5).
Writes per-country CSVs and a combined CSV. No plotting.
"""

import subprocess
import pandas as pd
import os

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
gene_aliases = {
    "CRT": ["CRT"],
    "K13": ["K13"],
    "MDR1": ["MDR1"],
    "DHFR": ["DHFR", "DHFR-TS"],
    "DHPS": ["DHPS", "PPPK-DHPS"],
    "PX1": ["PX1"],
    "UBP1": ["UBP1"],
    "AP2MU": ["AP2-MU", "AP2MU"],
    "AAT1": ["AAT1"],
}

def run_cmd(cmd):
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running {' '.join(cmd)}:\n{result.stderr}")
        return None
    return result.stdout

def annotate_vcf(vcf, ref_fasta, gff3, out_annotated_vcf):
    """
    Annotate with bcftools csq and index.
    """
    cmd = ['bcftools', 'csq', '-p', 'a', '-f', ref_fasta, '-g', gff3,
           '-o', out_annotated_vcf, '-O', 'z', vcf]
    if run_cmd(cmd) is None:
        return False
    return run_cmd(['bcftools', 'index', out_annotated_vcf]) is not None

def query_with_gt_dp(vcf, region):
    """
    Query GT:DP per sample.
    """
    fmt = '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%DP]\n'
    cmd = ['bcftools', 'query', '-f', fmt, '-r', region, vcf]
    return run_cmd(cmd)

def parse_output_and_recalculate_af(output, gene_name, min_dp=5):
    """
    Recalculate AF using ONLY genotypes with DP >= min_dp.
    """
    rows = []
    if not output:
        return pd.DataFrame(rows)

    aliases = gene_aliases.get(gene_name, [gene_name])

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

        # now filter to missense + correct gene
        for eff in csq.split(','):
            parts = eff.split('|')
            if len(parts) < 6:
                continue
            consequence = parts[0].lower()
            gene = parts[1]
            aa_change = parts[5] if len(parts) > 5 else f"{ref}>{alt}"
            if "missense" in consequence and gene in aliases:
                rows.append({
                    "gene": gene_name,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "aa_change": aa_change,
                    "AF": af
                })

    return pd.DataFrame(rows)

def run(vcf, ref_fasta, gff3, metadata_path, outdir,
        min_dp=5, write_combined=True):
    """
    CLI-facing entrypoint.
    """
    os.makedirs(outdir, exist_ok=True)

    metadata = pd.read_csv(metadata_path, sep="\t")
    if 'sample_id' not in metadata.columns or 'country' not in metadata.columns:
        raise SystemExit("Metadata must have 'sample_id' and 'country' columns.")

    country_to_samples = metadata.groupby('country')['sample_id'].apply(list).to_dict()
    all_data = []

    for country, sample_ids in country_to_samples.items():
        print(f"\nProcessing samples from country: {country} ...")
        country_dir = os.path.join(outdir, f"{country}")
        os.makedirs(country_dir, exist_ok=True)
        country_data = []

        for gene, region in genes.items():
            print(f"  Gene {gene} at region {region} ...")

            subset_vcf_file = os.path.join(country_dir, f"{gene}_subset.vcf.gz")
            annotated_vcf_file = os.path.join(country_dir, f"{gene}_subset_csq.vcf.gz")

            # subset
            cmd = [
                'bcftools', 'view',
                '-r', region,
                '-s', ','.join(sample_ids),
                '-Oz', '-o', subset_vcf_file,
                vcf
            ]
            if run_cmd(cmd) is None or run_cmd(['bcftools', 'index', subset_vcf_file]) is None:
                continue

            # annotate
            if not annotate_vcf(subset_vcf_file, ref_fasta, gff3, annotated_vcf_file):
                continue

            # query & calc
            output = query_with_gt_dp(annotated_vcf_file, region)
            if output is None:
                continue

            df = parse_output_and_recalculate_af(output, gene_name=gene, min_dp=min_dp)
            if df.empty:
                continue

            df["country"] = f"{country} (n={len(sample_ids)})"
            country_data.append(df)

        if country_data:
            country_df = pd.concat(country_data, ignore_index=True)
            out_csv = os.path.join(country_dir, "missense_AF_per_gene.csv")
            country_df.to_csv(out_csv, index=False)
            print(f"  Saved CSV for country {country} to {out_csv}")
            all_data.append(country_df)

    if not all_data:
        print("No data generated (after DP filtering).")
        return

    if write_combined:
        combined = pd.concat(all_data, ignore_index=True).sort_values(
            by=["country", "gene", "pos"]
        )
        combined_out = os.path.join(outdir, "missense_AF_per_gene_all_countries.csv")
        combined.to_csv(combined_out, index=False)
        print(f"Saved combined CSV to {combined_out}")
