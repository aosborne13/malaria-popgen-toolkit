# malaria-popgen-toolkit
A collection of reproducible population genomics workflows for *Plasmodium* species, designed to support transparent, accessible, and scalable malaria genomic analysis.

## Installation and Prerequisites
This toolkit has been formatted to interpret and analyse whole genome sequencing (WGS) data that has been processed using an open access *ngs-pipeline* toolkit designed by **VivaxGEN**, available at: https://github.com/vivaxgen/ngs-pipeline

Courtesy of Hidayat Trimarsanto (Anto) and Ludwig Kian Soon Hoon.

Install VivaxGEN NGS pipeline prerequisites (if you’re using those outputs):
```
"${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/ngs-pipeline/main/install.sh)

```
Filtered WGS data following mapping and QC using alternative pipelines can also be used with this toolkit. This includes outputs such as those from the *fastq2matrix* pipeline, courtesy of the **LSHTMPathogenSeqLab**, available at: https://github.com/LSHTMPathogenSeqLab/fastq2matrix

For vcf files, the csq field is prefered which can be generated using *bcftools*, available at: https://github.com/samtools/bcftools.

### Install **malaria-popgen-toolkit**.
```
pip install malaria-popgen-toolkit
# or from source while developing:
# pip install -e .

```

## Usage and Documentation (work in progress)
Further information and tutorials coming soon.

### Quick Start
How to use **malaria-popgen-toolkit**.

```
malaria-pipeline run \
  --vcf path/to/variants.vcf.gz \
  --ref Pf3D7.fasta \
  --gff3 Pf3D7.gff3 \
  --metadata samples.tsv \
  --outdir out \
  --min-dp 10

```
