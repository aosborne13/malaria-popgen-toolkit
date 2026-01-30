# SPDX-License-Identifier: Apache-2.0
"""
Main CLI entry point for the malaria-popgen-toolkit.

Program name:
  malaria-pipeline
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

from malaria_popgen_toolkit.resources.resolver import resolve_species


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
    p.add_argument("--vcf", required=True, help="Input VCF/BCF (bgzipped) with BCSQ or annotatable by bcftools csq")
    p.add_argument(
        "--species",
        default="Pf3D7",
        help="Species/reference bundle ID to fetch (default: Pf3D7)",
    )
    # Optional overrides (if provided, they win over resolver)
    p.add_argument("--ref", default=None, help="Override reference FASTA path (otherwise resolved from --species)")
    p.add_argument("--gff3", default=None, help="Override GFF3 annotation path (otherwise resolved from --species)")
    p.add_argument("--metadata", required=True, help="Metadata TSV with 'sample_id' and grouping column")
    p.add_argument("--outdir", default="missense_af", help="Output directory")
    p.add_argument("--min-dp", type=int, default=5, help="Minimum DP to include a sample genotype at a site")
    p.add_argument("--group-by", default="country", help="Metadata column to group by (e.g. country/region/year)")

    # ------------------------------------------------------------------
    # Haplotype map by region (your existing command module)
    # ------------------------------------------------------------------
    p = sub.add_parser("haplotype-map-region", help="Haplotype map by region/category")
    p.add_argument("--vcf", required=True, help="Input VCF/BCF (bgzipped)")
    p.add_argument("--metadata", required=True, help="Metadata TSV")
    p.add_argument("--outdir", default="haplotype_map", help="Output directory")
    p.add_argument(
        "--region",
        required=True,
        choices=["africa", "samerica", "south_america", "seasia", "southeast_asia"],
        help="Region preset (allowed): africa, samerica/south_america, seasia/southeast_asia",
    )
    p.add_argument("--min-dp", type=int, default=5, help="Minimum DP to include a sample genotype at a site")
    p.add_argument("--sample-col", default="sample_id")
    p.add_argument("--country-col", default="country")

    # Back-compat aliases (old command names)
    p_af = sub.add_parser("hapmap-africa", help=argparse.SUPPRESS)
    p_sa = sub.add_parser("hapmap-samerica", help=argparse.SUPPRESS)
    p_se = sub.add_parser("hapmap-seasia", help=argparse.SUPPRESS)

    def _add_hapmap_args(pp: argparse.ArgumentParser) -> None:
        pp.add_argument("--vcf", required=True, help="Input VCF/BCF (bgzipped)")
        pp.add_argument("--metadata", required=True, help="Metadata TSV")
        pp.add_argument("--outdir", default="haplotype_map", help="Output directory")
        pp.add_argument("--min-dp", type=int, default=5, help="Minimum DP to include a sample genotype at a site")
        pp.add_argument("--sample-col", default="sample_id")
        pp.add_argument("--country-col", default="country")

    _add_hapmap_args(p_af)
    _add_hapmap_args(p_sa)
    _add_hapmap_args(p_se)

    # ------------------------------------------------------------------
    # Fws dotplot
    # ------------------------------------------------------------------
    p = sub.add_parser("fws-dotplot", help="Fws jittered dot plots")
    p.add_argument("--metadata", required=True, help="Metadata TSV")
    p.add_argument("--outdir", default="fws_plots")
    p.add_argument(
        "--group-by",
        action="append",
        help="Metadata column(s) to group by; repeat for multiple. If omitted, tries region, country, year.",
    )
    p.add_argument("--width", type=float, default=10)
    p.add_argument("--height", type=float, default=6)

    # ------------------------------------------------------------------
    # PCA / PCoA
    # ------------------------------------------------------------------
    p = sub.add_parser("pca-plot", help="Distance-based PCA / PCoA from a matrix or VCF")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--matrix", help="Binary genotype matrix (.tsv) with samples as columns")
    g.add_argument("--vcf", help="Multi-sample VCF/BCF (bgzipped)")
    p.add_argument("--metadata", required=True, help="Metadata TSV")
    p.add_argument("--outdir", default="pca_plots")
    p.add_argument("--sample-col", default="sample_id")
    p.add_argument("--group-by", nargs="+", help="Metadata columns to color by (e.g. country region year)")
    p.add_argument("--pcs", nargs="+", help="PC pairs to plot, e.g. --pcs 1,2 1,3 (default: 1,2 and 1,3)")
    p.add_argument(
        "--max-sample-missing",
        type=float,
        default=None,
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
    p.add_argument("--matrix", required=True, help="Binary SNP matrix: chr,pos,ref,<samples>")
    p.add_argument("--metadata", required=True, help="Metadata TSV with sample IDs, category, and Fws column")
    p.add_argument("--outdir", required=True, help="Working/output directory for hmmIBD inputs/outputs")
    p.add_argument("--category-col", default="country")
    p.add_argument(
        "--category",
        default=None,
        help="Category value(s) (comma-separated) within --category-col. If omitted/ALL, run all categories.",
    )
    p.add_argument("--subgroup-col", default=None, help="Optional second-level grouping (e.g. year)")
    p.add_argument("--sample-col", default="sample_id")
    p.add_argument("--fws-col", default="fws")
    p.add_argument("--fws-th", type=float, default=0.95)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--hmmibd-bin", default="hmmIBD")
    p.add_argument("--skip-hmmibd", action="store_true")
    p.add_argument("--na-char", default="N", help="Character used for missing genotypes in matrix (also '.' treated missing)")
    p.add_argument("--exclude-chr", default="Pf3D7_API_v3,Pf3D7_MIT_v3", help="Comma-separated chromosomes to drop")
    p.add_argument("--regex-chr", default="(.*?)_(.+)_(.*)", help="Regex to parse chromosome names")
    p.add_argument("--regex-group", type=int, default=3, help="Capture group index in regex-chr that contains numeric chromosome")

    # ------------------------------------------------------------------
    # hmmIBD summary: windows + annotation
    # ------------------------------------------------------------------
    p = sub.add_parser(
        "hmmibd-summary",
        help="Summarise hmmIBD output into sliding windows and annotate genes",
    )
    p.add_argument("--workdir", required=True, help="Directory with hmmIBD outputs and ibd_matrix_hap_leg.tsv")
    p.add_argument(
        "--species",
        default="Pf3D7",
        help="Species/reference bundle ID to fetch (default: Pf3D7)",
    )
    # Optional overrides (if provided, they win)
    p.add_argument("--ref_index", default=None, help="Override FASTA .fai path (otherwise resolved from --species)")
    p.add_argument("--gene_product", default=None, help="Override genome product TSV path (otherwise resolved from --species)")
    p.add_argument("--suffix", required=True, help="Prefix for output files (e.g. 10_12_2025)")
    p.add_argument("--window_size", type=int, default=50000)
    p.add_argument("--maf", type=float, default=0.01)
    p.add_argument("--quantile_cutoff", type=float, default=0.95)
    p.add_argument("--remove_chr", default="Pf3D7_API_v3,Pf3D7_MIT_v3", help="Comma-separated chromosomes to drop")
    p.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p.add_argument("--regex_groupid", type=int, default=3)

    # ------------------------------------------------------------------
    # hmmIBD plots (boxplot, genome-wide, chromosome painting)
    # ------------------------------------------------------------------
    p_plot = sub.add_parser("hmmibd-ibdplots", help="Generate IBD plots from summarised hmmIBD windows")
    p_alias = sub.add_parser("hmmibd-ibdplot", help=argparse.SUPPRESS)

    def _add_hmmibd_plot_args(pp: argparse.ArgumentParser) -> None:
        pp.add_argument("--workdir", required=True)
        pp.add_argument("--species", default="Pf3D7", help="Species/reference bundle ID to fetch (default: Pf3D7)")
        pp.add_argument("--ref_index", default=None, help="Override FASTA .fai path (otherwise resolved from --species)")
        pp.add_argument("--gene_product", default=None, help="Override gene product TSV path (otherwise resolved from --species)")
        pp.add_argument("--suffix", required=True)
        pp.add_argument("--window_size", type=int, default=50000)
        pp.add_argument("--quantile_cutoff", type=float, default=0.95)
        pp.add_argument("--remove_chr", default="Pf3D7_API_v3,Pf3D7_MIT_v3", help="Comma-separated chromosomes to remove")
        pp.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
        pp.add_argument("--regex_groupid", type=int, default=3)
        pp.add_argument(
            "--outdir",
            default=None,
            help="Optional explicit output directory for plots (default: <workdir>/win_<window_kb>kb)",
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
    p.add_argument(
        "--species",
        default="Pf3D7",
        help="Species/reference bundle ID to fetch (default: Pf3D7)",
    )
    p.add_argument("--matrix_binary", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument(
        "--genome-file",
        default=None,
        help="Override genome product TSV path (otherwise resolved from --species)",
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
    p.add_argument(
        "--species",
        default="Pf3D7",
        help="Species/reference bundle ID to fetch (default: Pf3D7)",
    )
    p.add_argument(
        "--genome-file",
        default=None,
        help="Override genome product TSV path (otherwise resolved from --species)",
    )
    p.add_argument("--focus-pop", required=True)
    p.add_argument("--min-abs-xpehh", type=float, default=2.0)
    p.add_argument("--min-logp", type=float, default=1.3)
    p.add_argument("--remove_chr", default=None)
    p.add_argument("--regex_chr", default="(.*?)_(.+)_(.*)")
    p.add_argument("--regex_groupid", type=int, default=3)
    p.add_argument(
        "--panel-groups",
        default=None,
        help='Optional comma-separated list of comparisons, e.g. "Ethiopia_vs_Ghana,Ethiopia_vs_Malawi"',
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
        refs = resolve_species(args.species)
        ref_fasta = args.ref or refs["fasta"]
        gff3 = args.gff3 or refs["gff3"]

        missense_drugres_af.run(
            vcf=args.vcf,
            ref_fasta=ref_fasta,
            gff3=gff3,
            metadata_path=args.metadata,
            outdir=args.outdir,
            min_dp=args.min_dp,
            group_by=args.group_by,
        )
        return

    if args.command in ("haplotype-map-region", "hapmap-africa", "hapmap-samerica", "hapmap-seasia"):
        require_tool("bcftools")

        if args.command == "hapmap-africa":
            region = "africa"
        elif args.command == "hapmap-samerica":
            region = "samerica"
        elif args.command == "hapmap-seasia":
            region = "seasia"
        else:
            region = args.region

        haplotype_map_region.run(
            vcf=args.vcf,
            metadata_path=args.metadata,
            outdir=args.outdir,
            region=region,
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

    if args.command == "pca-plot":
        pca_plot.run(
          matrix=args.matrix,
          vcf=args.vcf,
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
            threads=args.threads,
            hmmibd_bin=args.hmmibd_bin,
            skip_hmmibd=args.skip_hmmibd,
            na_char=args.na_char,
            exclude_chr=args.exclude_chr,
            regex_chr=args.regex_chr,
            regex_group=args.regex_group,
        )
        return

    if args.command == "hmmibd-summary":
        require_tool("Rscript")
        refs = resolve_species(args.species)
        ref_index = args.ref_index or refs["fai"]
        gene_product = args.gene_product or refs["gene_product"]

        hmmibd_summary.run(
            workdir=args.workdir,
            ref_index=ref_index,
            gene_product=gene_product,
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
        refs = resolve_species(args.species)
        ref_index = args.ref_index or refs["fai"]
        gene_product = args.gene_product or refs["gene_product"]

        hmmibd_ibdplots.run(
            workdir=args.workdir,
            ref_index=ref_index,
            gene_product=gene_product,
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
        refs = resolve_species(args.species)
        genome_file = args.genome_file or refs["gene_product"]

        ihs_selection.run(
            workdir=args.workdir,
            matrix_binary=args.matrix_binary,
            metadata=args.metadata,
            genome_file=genome_file,
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
        refs = resolve_species(args.species)
        genome_file = args.genome_file or refs["gene_product"]

        xpehh_selection.run(
            workdir=args.workdir,
            genome_file=genome_file,
            focus_pop=args.focus_pop,
            min_abs_xpehh=args.min_abs_xpehh,
            min_logp=args.min_logp,
            remove_chr=args.remove_chr,
            regex_chr=args.regex_chr,
            regex_groupid=args.regex_groupid,
            panel_groups=args.panel_groups,
        )
        return

    parser.error(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()

