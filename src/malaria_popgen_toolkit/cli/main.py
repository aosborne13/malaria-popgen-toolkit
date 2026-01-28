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
    (alias: hmmibd-ibdplot)

Selection workflow (rehh):
  - ihs-selection         iHS pipeline from binary SNP matrix (optionally by subgroup)
  - xpehh-selection       XP-EHH comparisons from scanned_haplotypes_<pop>.tsv
"""

from __future__ import annotations

import argparse
import shutil
import sys

from malaria_popgen_toolkit.commands import (
    dataset_stats,
    fws_dotplot,
    haplotype_map_region,
    hmmibd_ibdplots,
    hmmibd_matrix,
    hmmibd_summary,
    ihs_selection,
    missense_drugres_af,
    pca_plot,
    xpehh_selection,
)


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        sys.exit(f"ERROR: Required tool '{name}' not found in PATH.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="malaria-pipeline",
        description="Malaria population genomics command-line toolkit",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ------------------------------------------------------------------
    # missense-drugres-af
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "missense-drugres-af",
        help="Compute missense allele frequencies in drug-resistance genes",
    )
    p.add_argument("--vcf", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--gff3", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", default="missense_af")
    p.add_argument("--min-dp", type=int, default=5)
    p.add_argument("--group-by", default="country")

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
    p = sub.add_parser("fws-dotplot", help="Fws jittered dot plots")
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", default="fws_plots")
    p.add_argument(
        "--group-by",
        action="append",
        help="Metadata column(s) to group by; repeat for multiple. "
        "If omitted, tries region, country, year.",
    )
    p.add_argument("--width", type=float, default=10)
    p.add_argument("--height", type=float, default=6)

    # ------------------------------------------------------------------
    # PCA / PCoA
    # ------------------------------------------------------------------
    p = sub.add_parser("pca", help="Distance-based PCA / PCoA")
    p.add_argument(
        "--matrix",
        required=True,
        help="Binary genotype matrix (.tsv) with samples as columns",
    )
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", default="pca_plots")
    p.add_argument("--sample-col", default="sample_id")
    p.add_argument(
        "--group-by",
        nargs="+",
        help="Metadata columns to color by (e.g. country region year)",
    )
    p.add_argument(
        "--pcs",
        nargs="+",
        help="PC pairs to plot, e.g. --pcs 1,2 1,3 (default: 1,2 and 1,3)",
    )
    p.add_argument(
        "--max-sample-missing",
        type=float,
        help="Drop samples with >X fraction missing (e.g. 0.3)",
    )

    # ------------------------------------------------------------------
    # dataset-stats
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "dataset-stats",
        help="Quickly report #samples and #variants from a VCF or binary matrix",
    )
    p.add_argument("--vcf", help="VCF/BCF (bgzipped) file")
    p.add_argument("--matrix", help="Binary matrix (.tsv) with samples as columns")

    # ------------------------------------------------------------------
    # hmmIBD (matrix -> input + run)
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "hmmibd-matrix",
        help="Prepare input and run hmmIBD from a binary SNP matrix",
    )
    p.add_argument(
        "--matrix",
        required=True,
        help="Binary SNP matrix with columns: chr,pos,ref,<samples>",
    )
    p.add_argument(
        "--metadata",
        required=True,
        help="Metadata TSV with sample IDs, category, and Fws column",
    )
    p.add_argument(
        "--outdir",
        required=True,
        help="Working/output directory for hmmIBD inputs/outputs",
    )
    p.add_argument("--category-col", default="country")
    p.add_argument(
        "--category",
        default=None,
        help="Category value(s) (comma-separated) within --category-col. "
        "If omitted/ALL, run all categories.",
    )
    p.add_argument(
        "--subgroup-col",
        default=None,
        help="Optional second-level grouping (e.g. year) to split within category",
    )
    p.add_argument("--sample-col", default="sample_id")
    p.add_argument("--fws-col", default="fws")
    p.add_argument("--fws-th", type=float, default=0.95)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--hmmibd-bin", default="hmmIBD")
    p.add_argument("--skip-hmmibd", action="store_true")
    p.add_argument(
        "--na-char",
        default="N",
        help="Character used for missing genotypes in matrix (also '.' treated as missing)",
    )
    p.add_argument(
        "--exclude-chr",
        default="Pf3D7_API_v3,Pf3D7_MIT_v3",
        help="Comma-separated list of chromosomes to drop",
    )
    p.add_argument(
        "--regex-chr",
        default="(.*?)_(.+)_(.*)",
        help="Regex to parse chromosome names, default matches Pf3D7_01_v3",
    )
    p.add_argument(
        "--regex-group",
        type=int,
        default=3,
        help="Capture group index in regex_chr that contains numeric chromosome",
    )

    # ------------------------------------------------------------------
    # hmmIBD summary: windows + annotation
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "hmmibd-summary",
        help="Summarise hmmIBD output into sliding windows and annotate genes",
    )
    p.add_argument(
        "--workdir",
        required=True,
        help="Directory with hmmIBD outputs and ibd_matrix_hap_leg.tsv",
    )
    p.add_argument("--ref_index", required=True, help="Reference FASTA .fai index (Pf3D7)")
    p.add_argument(
        "--gene_product",
        required=True,
        help="Pf3D7 genome product TSV annotation (pf_genome_product_v3.tsv)",
    )
    p.add_argument("--suffix", required=True, help="Prefix for output files (e.g. 10_12_2025)")
    p.add_argument("--window_size", type=int, default=50000)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--quantile_cutoff", type=float, default=0.95)
    p.add_argument(
        "--remove_chr",
        default="Pf3D7_API_v3,Pf3D7_MIT_v3",
        help="Comma-separated list of chromosomes to drop",
    )
    p.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p.add_argument("--regex_groupid", type=int, default=3)

    # ------------------------------------------------------------------
    # hmmIBD plots (boxplot, genome-wide, chromosome painting)
    # ------------------------------------------------------------------
    p_plot = sub.add_parser(
        "hmmibd-ibdplots",
        help="Generate IBD plots from summarised hmmIBD windows",
    )
    p_alias = sub.add_parser("hmmibd-ibdplot", help=argparse.SUPPRESS)

    def _add_hmmibd_plot_args(pp: argparse.ArgumentParser) -> None:
        pp.add_argument("--workdir", required=True)
        pp.add_argument("--ref_index", required=True)
        pp.add_argument("--gene_product", required=True)
        pp.add_argument("--suffix", required=True)
        pp.add_argument("--window_size", type=int, default=50000)
        pp.add_argument("--quantile_cutoff", type=float, default=0.95)
        pp.add_argument(
            "--remove_chr",
            default="Pf3D7_API_v3,Pf3D7_MIT_v3",
            help="Comma-separated chromosomes to remove",
        )
        pp.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
        pp.add_argument("--regex_groupid", type=int, default=3)
        pp.add_argument(
            "--outdir",
            default=None,
            help="Optional explicit output directory for plots "
            "(default: <workdir>/win_<window_kb>kb)",
        )

    _add_hmmibd_plot_args(p_plot)
    _add_hmmibd_plot_args(p_alias)

    # ------------------------------------------------------------------
    # iHS selection
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "ihs-selection",
        help="Run iHS scans from a binary SNP matrix (per category/subgroup) and generate plots/TSVs",
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--matrix_binary", required=True)
    p.add_argument("--metadata", required=True)
    # NOTE: --annotation has been removed; R script now builds marker map from the matrix.
    p.add_argument(
        "--genome-file",
        required=True,
        help="Genome product TSV with columns: chr,pos_start,pos_end,gene_id,gene_product,gene_name",
    )
    p.add_argument("--label_category", default="country")
    p.add_argument("--subgroup_col", default=None)
    p.add_argument("--label_id", default="sample_id")
    p.add_argument("--label_fws", default="fws")
    p.add_argument(
        "--category",
        default=None,
        help="Category value(s) (comma-separated) within --label_category. If omitted/ALL, run all.",
    )
    p.add_argument("--focus-pop", default=None)
    p.add_argument("--fws_th", type=float, default=0.95)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--rehh_min_perc_hap", type=float, default=80.0)
    p.add_argument("--rehh_min_perc_mrk", type=float, default=70.0)
    p.add_argument("--na_char", default="NA")
    p.add_argument("--forced_recode", action="store_true")
    p.add_argument("--forced_mixed", action="store_true")
    p.add_argument("--remove_chr", default=None)
    p.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p.add_argument("--regex_groupid", type=int, default=3)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-maf-ihs", type=float, default=0.0)
    p.add_argument("--freqbin", type=float, default=0.05)
    p.add_argument("--ihs-thresh", type=float, default=2.0)
    p.add_argument("--logp-thresh", type=float, default=5.0)

    # ------------------------------------------------------------------
    # XP-EHH selection
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "xpehh-selection",
        help="Run XP-EHH comparisons from scan_hh outputs (scanned_haplotypes_<pop>.tsv)",
    )
    p.add_argument("--workdir", required=True)
    p.add_argument("--genome-file", required=True)
    p.add_argument("--focus-pop", required=True)
    p.add_argument("--min-abs-xpehh", type=float, default=2.0)
    p.add_argument("--min-logp", type=float, default=1.3)
    p.add_argument("--remove_chr", default=None)
    p.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p.add_argument("--regex_groupid", type=int, default=3)
    p.add_argument(
        "--panel-groups",
        default=None,
        help="Optional comma-separated list of comparisons for a stacked multi-panel plot, "
        'e.g. "Ethiopia_vs_Ghana,Ethiopia_vs_Malawi,Ethiopia_vs_Sudan"',
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

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
        return

    if args.command in ("hapmap-africa", "hapmap-samerica", "hapmap-seasia"):
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
        return

    if args.command == "fws-dotplot":
        fws_dotplot.run(
            metadata_path=args.metadata,
            outdir=args.outdir,
            group_by=args.group_by,
            width=args.width,
            height=args.height,
        )
        return

    if args.command == "pca":
        pca_plot.run(
            matrix=args.matrix,
            metadata_path=args.metadata,
            outdir=args.outdir,
            sample_col=args.sample_col,
            group_by=args.group_by,
            pcs=args.pcs,
            max_sample_missing=args.max_sample_missing,
        )
        return

    if args.command == "dataset-stats":
        dataset_stats.run(vcf=args.vcf, matrix=args.matrix)
        return

    if args.command == "hmmibd-matrix":
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
        return

    if args.command == "hmmibd-summary":
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
        return

    if args.command in ("hmmibd-ibdplots", "hmmibd-ibdplot"):
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
        return

    if args.command == "ihs-selection":
        require_tool("Rscript")
        ihs_selection.run(
            workdir=args.workdir,
            matrix_binary=args.matrix_binary,
            metadata_path=args.metadata,
            # NOTE: annotation_path removed
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
        return

    if args.command == "xpehh-selection":
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
        return

    # Should never get here due to required=True on subparsers
    parser.error(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()


