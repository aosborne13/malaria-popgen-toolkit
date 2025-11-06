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

COMING SOON...

## Drug Resistance - *P. falciparum*
Candidates for *P. vivax* coming soon...

### Missense variants on validated and candidate drug resistance genes
Compute allele frequencies for missense variants in drug-resistance genes, grouping by any metadata column (e.g. country, region, site, or year).

The metadata file **must have a sample_id column** *and* **chosen grouping column**.

#### Included markers...

Validated: K13, CRT, MDR1, DHFR, DHPS

Candidate: AAT1, PX1, UBP1, AP2MU

```
malaria-pipeline missense-drugres-af \
  --vcf path/to/filtered.vcf.gz \
  --ref path/to/reference.fasta \
  --gff3 path/to/annotation.gff3 \
  --metadata path/to/metadata.tsv \
  --outdir results_missense_AF \
  --min-dp 5 \
  --group-by country

```
##### Output example
Coming soon...

### Haplotype analysis for CRT, MDR1, DHFR, and DHPS

```
malaria-pipeline hapmap-africa \
  --matrix path/to/matrix.tsv \
  --metadata path/to/metadata.tsv \
  --outdir hapmap_africa_output \
  --sample-col sample_id \
  --country-col country
```

#### Example output: African haplotype map

Below is an example output from the `hapmap-africa` command:

![Example African haplotype map](docs/images/hapmap_africa_example.png)
