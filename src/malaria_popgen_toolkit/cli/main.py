# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.

Program name:
  malaria-pipeline

Core commands:
  - missense-drugres-af   Allele frequencies of missense drug-resistance variants
  - hapmap-africa         Haplotype map for Africa
  - hapmap-samerica       Haplotype map for South America
  - hapmap-seasia         Haplotype map for Southeast Asia
  - fws-dotplot           Fws jittered dot plots by metadata groups
  - pca                   Distance-based PCA / PCoA from binary matrix
  - dataset-stats         Quick sample/variant counts for VCF or matrix

IBD workflow (hmmIBD):
  - hmmibd-matrix         Prepare and run hmmIBD from a binary SNP matrix
  - hmmibd-summary        Summarise hmmIBD output into sliding windows + annotate genes
  - hmmibd-ibdplots       Plot IBD summaries (boxplot, genome-wide, chromosome painting)

Selection workflow (rehh):
  - ihs-selection         iHS pipeline from binary SNP matrix (optionally by subgroup)
  - xpehh-selection       XP-EHH comparisons from scanned_haplotypes_<pop>.tsv
"""

import argparse
import sys
import shutil

from malaria_popgen_toolkit.commands import (
    missense_drugres_af,
    haplotype_map_region,
    fws_dotplot,
    pca_plot,
    hmmibd_matrix,
    hmmibd_summary,
    hmmibd_ibdplots,
    ihs_selection,
    xpehh_selection,
    dataset_stats,
)


def require_tool(name: str):
    if shutil.which(name) is None:
        sys.exit(f"ERROR: Required tool '{name}' not found in PATH.")


def main():
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ------------------------------------------------------------------
    # missense-drugres-af
    # ------------------------------------------------------------------
    p1 = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes",
    )
    p1.add_argument("--vcf", required=True)
    p1.add_argument("--ref", required=True)
    p1.add_argument("--gff3", required=True)
    p1.add_argument("--metadata", required=True)
    p1.add_argument("--outdir", default="missense_af")
    p1.add_argument("--min-dp", type=int, default=5)
    p1.add_argument("--group-by", default="country")

    # ------------------------------------------------------------------
    # Haplotype maps
    # ------------------------------------------------------------------
    for cmd, region in [
        ("hapmap-africa", "africa"),
        ("hapmap-samerica", "south_america"),
        ("hapmap-seasia", "southeast_asia"),
    ]:
        p = sub.add_parser(cmd, help=f"Haplotype map for {region}")
        p.add_argument("--vcf", required=True)
        p.add_argument("--metadata", required=True)
        p.add_argument("--outdir", required=True)
        p.add_argument("--min-dp", type=int, default=5)
        p.add_argument("--sample-col", default="sample_id")
        p.add_argument("--country-col", default="country")

    # ------------------------------------------------------------------
    # Fws dotplot
    # ------------------------------------------------------------------
    p5 = sub.add_parser("fws-dotplot", help="Fws jittered dot plots")
    p5.add_argument("--metadata", required=True)
    p5.add_argument("--outdir", default="fws_plots")
    p5.add_argument(
        "--group-by",
        action="append",
        help="Metadata column(s) to group by; repeat for multiple. "
        "If omitted, tries region, country, year.",
    )
    p5.add_argument("--width", type=float, default=10)
    p5.add_argument("--height", type=float, default=6)

    # ------------------------------------------------------------------
    # PCA / PCoA
    # ------------------------------------------------------------------
    p6 = sub.add_parser("pca", help="Distance-based PCA / PCoA")
    p6.add_argument(
        "--matrix",
        required=True,
        help="Binary genotype matrix (.tsv) with samples as columns",
    )
    p6.add_argument("--metadata", required=True)
    p6.add_argument("--outdir", default="pca_plots")
    p6.add_argument("--sample-col", default="sample_id")
    p6.add_argument(
        "--group-by",
        nargs="+",
        help="Metadata columns to color by (e.g. country region year)",
    )
    p6.add_argument(
        "--pcs",
        nargs="+",
        help="PC pairs to plot, e.g. --pcs 1,2 1,3 (default: 1,2 and 1,3)",
    )
    p6.add_argument(
        "--max-sample-missing",
        type=float,
        help="Drop samples with >X fraction missing (e.g. 0.3)",
    )

    # ------------------------------------------------------------------
    # hmmIBD (matrix -> input + run)
    # ------------------------------------------------------------------
    p7 = sub.add_parser(
        "hmmibd-matrix",
        help="Prepare input and run hmmIBD from a binary SNP matrix",
    )
    p7.add_argument(
        "--matrix",
        required=True,
        help="Binary SNP matrix with columns: chr,pos,ref,<samples>",
    )
    p7.add_argument(
        "--metadata",
        required=True,
        help="Metadata TSV with sample IDs, category, and Fws column",
    )
    p7.add_argument(
        "--outdir",
        required=True,
        help="Working/output directory for hmmIBD inputs/outputs",
    )
    p7.add_argument("--category-col", default="country")
    p7.add_argument(
        "--category",
        default=None,
        help="Category name or comma-separated list (e.g. Ethiopia,Kenya). "
        "If omitted, run all categories.",
    )
    p7.add_argument(
        "--subgroup-col",
        default=None,
        help="Optional second-level grouping (e.g. year) to split within category",
    )
    p7.add_argument("--sample-col", default="sample_id")
    p7.add_argument("--fws-col", default="fws")
    p7.add_argument("--fws-th", type=float, default=0.95)
    p7.add_argument("--maf", type=float, default=0.01)
    p7.add_argument("--threads", type=int, default=4)
    p7.add_argument("--hmmibd-bin", default="hmmIBD")
    p7.add_argument("--skip-hmmibd", action="store_true")
    p7.add_argument(
        "--na-char",
        default="N",
        help="Character used for missing genotypes in matrix "
        "(also '.' treated as missing)",
    )
    p7.add_argument(
        "--exclude-chr",
        default="Pf3D7_API_v3,Pf3D7_MIT_v3",
        help="Comma-separated list of chromosomes to drop "
        "(e.g. Pf3D7_API_v3,Pf3D7_MIT_v3)",
    )
    p7.add_argument(
        "--regex-chr",
        default="(.*?)_(.+)_(.*)",
        help="Regex to parse chromosome names, default matches Pf3D7_01_v3",
    )
    p7.add_argument(
        "--regex-group",
        type=int,
        default=3,
        help="Capture group index in regex_chr that contains numeric chromosome",
    )

    # ------------------------------------------------------------------
    # hmmIBD summary: windows + annotation
    # ------------------------------------------------------------------
    p_hmm_sum = sub.add_parser(
        "hmmibd-summary",
        help="Summarise hmmIBD output into sliding windows and annotate genes",
    )
    p_hmm_sum.add_argument(
        "--workdir",
        required=True,
        help="Directory with hmmIBD outputs and ibd_matrix_hap_leg.tsv",
    )
    p_hmm_sum.add_argument(
        "--ref_index",
        required=True,
        help="Reference FASTA .fai index (Pf3D7)",
    )
    p_hmm_sum.add_argument(
        "--gene_product",
        required=True,
        help="Pf3D7 genome product TSV annotation (pf_genome_product_v3.tsv)",
    )
    p_hmm_sum.add_argument(
        "--suffix",
        required=True,
        help="Prefix for output files (e.g. 10_12_2025)",
    )
    p_hmm_sum.add_argument("--window_size", type=int, default=50000)
    p_hmm_sum.add_argument("--maf", type=float, default=0.01)
    p_hmm_sum.add_argument("--quantile_cutoff", type=float, default=0.95)
    p_hmm_sum.add_argument(
        "--remove_chr",
        default="Pf3D7_API_v3,Pf3D7_MIT_v3",
        help="Comma-separated list of chromosomes to drop",
    )
    p_hmm_sum.add_argument(
        "--regex_chr",
        default="(.*?)_(.+)_(.*)",
        help="Regex for extracting numeric chromosome from chr names",
    )
    p_hmm_sum.add_argument(
        "--regex_groupid",
        type=int,
        default=3,
        help="Capture group index in regex_chr giving numeric chromosome",
    )

    # ------------------------------------------------------------------
    # hmmIBD plots (boxplot, genome-wide, chromosome painting)
    # ------------------------------------------------------------------
    p_hmm_plot = sub.add_parser(
        "hmmibd-ibdplots",
        help="Generate IBD plots from summarised hmmIBD windows",
    )
    # Back-compat alias (old name)
    p_hmm_plot_alias = sub.add_parser(
        "hmmibd-ibdplot",
        help=argparse.SUPPRESS,
    )

    def _add_hmmibd_plot_args(p):
        p.add_argument(
            "--workdir",
            required=True,
            help="Directory containing <suffix>_hmmIBD_ibd_winXkb.tsv and <suffix>_hmmIBD_fraction.tsv",
        )
        p.add_argument(
            "--ref_index",
            required=True,
            help="Reference FASTA .fai index (Pf3D7)",
        )
        p.add_argument(
            "--gene_product",
            required=True,
            help="Pf3D7 genome product TSV annotation (pf_genome_product_v3.tsv)",
        )
        p.add_argument(
            "--suffix",
            required=True,
            help="Prefix used in hmmIBD summarisation, e.g. '10_12_2025'",
        )
        p.add_argument("--window_size", type=int, default=50000)
        p.add_argument("--quantile_cutoff", type=float, default=0.95)
        p.add_argument(
            "--remove_chr",
            default="Pf3D7_API_v3,Pf3D7_MIT_v3",
            help="Comma-separated chromosomes to remove [default: Pf3D7_API_v3,Pf3D7_MIT_v3]",
        )
        p.add_argument(
            "--regex_chr",
            default="(.*?)_(.+)_(.*)",
            help="Regex for extracting numeric chromosome from chr names",
        )
        p.add_argument(
            "--regex_groupid",
            type=int,
            default=3,
            help="Capture group index in regex_chr giving numeric chromosome",
        )
        p.add_argument(
            "--outdir",
            default=None,
            help="Optional explicit output directory for plots (default: <workdir>/win_<window_kb>kb)",
        )

    _add_hmmibd_plot_args(p_hmm_plot)
    _add_hmmibd_plot_args(p_hmm_plot_alias)

    # ------------------------------------------------------------------
    # dataset-stats
    # ------------------------------------------------------------------
    p8 = sub.add_parser(
        "dataset-stats",
        help="Quickly report #samples and #variants from a VCF or binary matrix",
    )
    p8.add_argument("--vcf", help="VCF/BCF (bgzipped) file")
    p8.add_argument("--matrix", help="Binary matrix (.tsv) with samples as columns")

    # ------------------------------------------------------------------
    # iHS selection
    # ------------------------------------------------------------------
    p10 = sub.add_parser(
        "ihs-selection",
        help="Run iHS scans from a binary SNP matrix (per category/subgroup) and generate plots/TSVs",
    )
    p10.add_argument("--workdir", required=True)
    p10.add_argument("--matrix_binary", required=True)
    p10.add_argument("--metadata", required=True)
    p10.add_argument("--annotation", required=True)
    p10.add_argument("--genome-file", required=True)
    p10.add_argument("--label_category", default="country")
    p10.add_argument("--subgroup_col", default=None)
    p10.add_argument("--label_id", default="sample_id")
    p10.add_argument("--label_fws", default="fws")
    p10.add_argument("--category", default=None)
    p10.add_argument("--focus-pop", default=None)
    p10.add_argument("--fws_th", type=float, default=0.95)
    p10.add_argument("--maf", type=float, default=0.01)
    p10.add_argument("--rehh_min_perc_hap", type=float, default=80.0)
    p10.add_argument("--rehh_min_perc_mrk", type=float, default=70.0)
    p10.add_argument("--na_char", default="NA")
    p10.add_argument("--forced_recode", action="store_true")
    p10.add_argument("--forced_mixed", action="store_true")
    p10.add_argument("--remove_chr", default=None)
    p10.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p10.add_argument("--regex_groupid", type=int, default=3)
    p10.add_argument("--threads", type=int, default=4)
    p10.add_argument("--min-maf-ihs", type=float, default=0.0)
    p10.add_argument("--freqbin", type=float, default=0.05)
    p10.add_argument("--ihs-thresh", type=float, default=2.0)
    p10.add_argument("--logp-thresh", type=float, default=5.0)

    # ------------------------------------------------------------------
    # XP-EHH selection
    # ------------------------------------------------------------------
    p11 = sub.add_parser(
        "xpehh-selection",
        help="Run XP-EHH comparisons from scan_hh outputs (scanned_haplotypes_<pop>.tsv)",
    )
    p11.add_argument("--workdir", required=True)
    p11.add_argument("--genome-file", required=True)
    p11.add_argument("--focus-pop", required=True)
    p11.add_argument("--min-abs-xpehh", type=float, default=2.0)
    p11.add_argument("--min-logp", type=float, default=1.3)
    p11.add_argument("--remove_chr", default=None)
    p11.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p11.add_argument("--regex_groupid", type=int, default=3)
    p11.add_argument(
        "--panel-groups",
        default=None,
        help="Optional comma-separated list of comparisons for a stacked multi-panel plot, "
        'e.g. "Ethiopia_vs_Ghana,Ethiopia_vs_Malawi,Ethiopia_vs_Sudan"',
    )

    args = parser.parse_args()

    # ==========================
    # Dispatch
    # ==========================

    if args.command == "missense-drugres-af":
        require_tool("bcftools")
        missense_drugres_af.run(
            vcf=args.vcf,
            ref_fasta=args.ref,
            gff3=args.gff3,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            group_by=args.group_by,
        )

    elif args.command in ("hapmap-africa", "hapmap-samerica", "hapmap-seasia"):
        require_tool("bcftools")
        region = {
            "hapmap-africa": "africa",
            "hapmap-samerica": "south_america",
            "hapmap-seasia": "southeast_asia",
        }[args.command]
        haplotype_map_region.run(
            region=region,
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            sample_col=args.sample_col,
            country_col=args.country_col,
        )

    elif args.command == "fws-dotplot":
        fws_dotplot.run(
            metadata_path=args.metadata,
            outdir=args.outdir,
            group_by=args.group_by,
            width=args.width,
            height=args.height,
        )

    elif args.command == "pca":
        pca_plot.run(
            matrix=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            group_by=args.group_by,
            pcs=args.pcs,
            max_sample_missing=args.max_sample_missing,
        )

    elif args.command == "hmmibd-matrix":
        require_tool("Rscript")
        require_tool(args.hmmibd_bin)
        hmmibd_matrix.run(
            matrix=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            category_col=args.category_col,
            category=args.category,
            subgroup_col=args.subgroup_col,
            sample_col=args.sample_col,
            fws_col=args.fws_col,
            fws_th=args.fws_th,
            maf=args.maf,
            na_char=args.na_char,
            threads=args.threads,
            exclude_chr=args.exclude_chr,
            regex_chr=args.regex_chr,
            regex_group=args.regex_group,
            hmmibd_bin=args.hmmibd_bin,
            skip_hmmibd=args.skip_hmmibd,
        )

    elif args.command == "hmmibd-summary":
        require_tool("Rscript")
        hmmibd_summary.run(
            workdir=args.workdir,
            ref_index=args.ref_index,
            gene_product=args.gene_product,
            suffix=args.suffix,
            window_size=args.window_size,
            maf=args.maf,
            quantile_cutoff=args.quantile_cutoff,
            remove_chr=args.remove_chr,
            regex_chr=args.regex_chr,
            regex_groupid=args.regex_groupid,
        )

    elif args.command in ("hmmibd-ibdplots", "hmmibd-ibdplot"):
        require_tool("Rscript")
        hmmibd_ibdplots.run(
            workdir=args.workdir,
            ref_index=args.ref_index,
            gene_product=args.gene_product,
            suffix=args.suffix,
            window_size=args.window_size,
            quantile_cutoff=args.quantile_cutoff,
            remove_chr=args.remove_chr,
            regex_chr=args.regex_chr,
            regex_groupid=args.regex_groupid,
            outdir=args.outdir,
        )

    elif args.command == "dataset-stats":
        dataset_stats.run(
            vcf=args.vcf,
            matrix=args.matrix,
        )

    elif args.command == "ihs-selection":
        require_tool("Rscript")
        ihs_selection.run(
            workdir=args.workdir,
            matrix_binary=args.matrix_binary,
            metadata_path=args.metadata,
            annotation_path=args.annotation,
            genome_file=args.genome_file,
            label_category=args.label_category,
            subgroup_col=args.subgroup_col,
            label_id=args.label_id,
            label_fws=args.label_fws,
            category=args.category,
            focus_pop=args.focus_pop,
            fws_th=args.fws_th,
            maf=args.maf,
            rehh_min_perc_hap=args.rehh_min_perc_hap,
            rehh_min_perc_mrk=args.rehh_min_perc_mrk,
            na_char=args.na_char,
            forced_recode=args.forced_recode,
            forced_mixed=args.forced_mixed,
            remove_chr=args.remove_chr,
            regex_chr=args.regex_chr,
            regex_groupid=args.regex_groupid,
            threads=args.threads,
            min_maf_ihs=args.min_maf_ihs,
            freqbin=args.freqbin,
            ihs_thresh=args.ihs_thresh,
            logp_thresh=args.logp_thresh,
        )

    elif args.command == "xpehh-selection":
        require_tool("Rscript")
        xpehh_selection.run(
            workdir=args.workdir,
            genome_file=args.genome_file,
            focus_pop=args.focus_pop,
            min_abs_xpehh=args.min_abs_xpehh,
            min_logp=args.min_logp,
            remove_chr=args.remove_chr,
            regex_chr=args.regex_chr,
            regex_groupid=args.regex_groupid,
            panel_groups=args.panel_groups,
        )


if __name__ == "__main__":
    main()
