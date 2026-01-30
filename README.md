# malaria-popgen-toolkit
A collection of reproducible population genomics workflows for *Plasmodium* species, designed to support transparent, accessible, and scalable malaria genomic analysis.

*Tools and documentation were made for WGS data. Amplicon sequencing is technically not supported at the moment but many of these tools should be compatible with it.*

*Undergoing renovations...* 
- Adding a ```--species``` tag and including genome reference files into the installation (e.g. Pf3D7 v3 for now, PvP01 v2 to come).

## Installation and Prerequisites
This toolkit has been formatted to interpret and analyse whole genome sequencing (WGS) data that has been processed using an open access *ngs-pipeline* toolkit designed by **VivaxGEN**, available at: 

https://github.com/vivaxgen/ngs-pipeline

Courtesy of Hidayat Trimarsanto (Anto) and Ludwig Kian Soon Hoon.

Downstream processing can also be carried out using the *fastq2matrix* and *malaria-hub* pipeline, courtesy of the **LSHTMPathogenSeqLab**, available at: 

https://github.com/LSHTMPathogenSeqLab/

For vcf files:
- Files should be filtered to only include the **core genome** - sub-telmomeric and hypervariable regions should be removed (a core genome bed file will be provided in a future release of this package - available through *MalariaGEN*)
- The csq field is prefered which can be generated using *bcftools*, available at: https://github.com/samtools/bcftools.

### Install **malaria-popgen-toolkit**.
Step 1: Install micromamba (if you don't already have it).
```
curl -L micro.mamba.pm/install.sh | bash
exec bash  # restart your shell so micromamba is available
```
Step 2: Clone this repository.
```
git clone https://github.com/aosborne13/malaria-popgen-toolkit.git
cd malaria-popgen-toolkit
```
Step 3: Create the environment.
This installs all prerequisites (Python, R packages, bcftools, hmmIBD, TESS3r, etc.) in one command. Requires a stead internet connection.
```
bash scripts/setup_env.sh
```
Step 4: Activate the new environment once the installation is complete.
```
micromamba activate malaria-popgen
```
Step 5: Verify installation.
```
malaria-pipeline --help
```
To update, activate the *malaria-popgen* enironment and then run the following:
```
git pull
python -m pip install -e .
```

# Usage and Documentation (work in progress)
## Formatting your input files

### Input *metadata.tsv* file:
The input metadata.tsv file should be a tab-deliminated file formatted like the example below. The column "sample_id" is required for all tools. A group, such as "country", "region", or "year" is required for many plotting tools to calculate or plot population-specific information. The geographical mapping tools require "country" information, "region" is not suffient (for now).

*The "**fws**" column is optional for some tools but required for others - documentation will be made clear when it is required. Information on how to calculate Fws is provided below.*

**COLUMN LABELS ARE CASE-SENSITIVE - "sample_id" can *NOT* be "SAMPLE_ID"**

Columns unique to your dataset (such as a sampling site) can be included in your *metadata.tsv* file.

*Remember not to include spaces in your column names - replace any spaces with an underscore "_" or a dash "-" to prevent errors. Example: Use "sample_site" or "sample-site", **NOT** "sample site"*

![Example Metadata](docs/images/metadata_example.png)

### Example *input.vcf.gz* file:

The "sample_id" column in your *metadata.tsv* file, must match the sample names, or IDs, in your VCF file **exactly**.

![Example VCF](docs/images/vcf_example.png)

# Drug Resistance - *P. falciparum*
Candidates for *P. vivax* coming soon...

## Missense variants on validated and candidate drug resistance genes
Compute allele frequencies for missense variants in drug-resistance genes, grouping by any metadata column (e.g. country, region, site, or year).

### Included markers...

Validated: K13, CRT, MDR1, DHFR, DHPS

Candidate: AAT1, PX1, UBP1, AP2MU

```
malaria-pipeline missense-drugres-af \
  --vcf path/to/filtered.vcf.gz \
  --species Pf3D7 \
  --metadata path/to/metadata.tsv \
  --outdir results_missense_AF \
  --min-dp 5 \
  --group-by country

```
*For users that want to supply their own reference and gff files, usage:*
```
malaria-pipeline missense-drugres-af \
  --vcf path/to/filtered.vcf.gz \
  --ref path/to/reference_genome.fasta \
  --gff3 path/to/gff_annotation_file.gff3 \
  --metadata path/to/metadata.tsv \
  --outdir results_missense_AF \
  --min-dp 5 \
  --group-by country

```

#### Output example:
A global dataset run using the "country" group label and a minimum read depth (DP) of 5.

![Example Drug Resistance](docs/images/drug_res_example.png)

## Haplotype analysis for CRT, MDR1, DHFR, and DHPS
Support currently includes maps of Africa, South America, and Southeast Asia. As of now, K13 has not been included in these maps as a wild-type haplotype has not been established.

Geographical region can be specified using the argument ```--region``` and currently accepts:
- Africa = ```africa```
- Southeast Asia = ```seasia```,```southeast_asia```
- South America = ```samerica```,```south_america```

Here is an example of how to run the commands to generate haplotype maps depending on the region you are working in:

```
malaria-pipeline haplotype-map-region \
  --vcf path/to/filtered.vcf.gz \
  --metadata path/to/metadata.tsv \
  --outdir hapmap_africa_output \
  --region africa \
  --min-dp 5 \
  --sample-col sample_id \
  --country-col country

```

### Example output of the African haplotype map

Below is an example output from the `hapmap-africa` command:

![Example African haplotype map](docs/images/hapmap_africa_example.png)

# Population Genomics - *PopGen* Tools
## Complexity of Infection (COI) - using *Fws* estimates
*Fws*, is a required metric for some downstream WGS processing. This can be calculated using *moimix* R-based package, available at: https://github.com/bahlolab/moimix

User-friendly R-script for running filtered VCF files through *moimix* available from the *malaria-hub* pipeline.

https://github.com/LSHTMPathogenSeqLab/malaria-hub/tree/master/moi

### Plot *Fws* by Population
Plotting *Fws* by population can be an intuitive way to visualise diversity between populations, either between different regions or the same region over different time points. Changes in diversity can give insight into transmission intensity. An *Fws* estimate > 0.95 generally correlates to a monoclonal infection and is more commonly identified in low transmission regions. An isolate with an *Fws* < 0.95 likely consists of multiple clones, with the diversity or number of clones increasing as he *Fws* value decreases.

The default plotting script will automatically search for columns labeled "country", "region", and "year" in your *metadata.tsv* file. A column with *Fws* values is **required** for this script - the column should be named **"fws"**.
```
malaria-pipeline fws-dotplot \
  --metadata /path/to/metadata.tsv \
  --outdir fws_plots
```
To select specific columns for plotting, or a group unique to your dataset, specify the group, or groups, using this option:

```
malaria-pipeline fws-dotplot \
  --metadata /path/to/metadata.tsv \
  --group-by clinical_trial --group-by sample_site \
  --outdir fws_plots_trial
```

To adjust the size and dimesnsions of your plot, use this option:
```
malaria-pipeline fws-dotplot \
  --metadata /path/to/metadata.tsv \
  --group-by country \
  --width 12 --height 6 \
  --outdir fws_plots
```
#### Output:
![Example fws](docs/images/fws_plot_example.png)

## Quick Stats - How many samples and variants are in my filtered input file?
Use this command to quickly report the number of samples and variants. You can provide either:

```--matrix``` — a binary matrix with values 0, 0.5, 1, or N

or

```--vcf``` — a multi-sample VCF (which is internally converted into numeric genotypes)

Example - running this tool using a binary matrix:
```
malaria-pipeline dataset-stats \
  --matrix popgen_africa.mat.bin
```
Running using a filtered VCF file:
```
malaria-pipeline dataset-stats \
  --vcf popgen_africa.vcf.gz
```
*Note - VCF runtime can be a long for large datasets, matrix will be faster*

## Distance-based PCA/PCoA
This command performs pairwise SNP-difference distances using a Manhattan metric:
- Missing genotypes (N or .) are allowed
- Loci missing in one sample are ignored for that pair
- Distances are scaled following the behavior of the R-package *amap* ```amap::Dist```

The PCoA / classical multidimensional-scaling (MDS) is run on the distance matrix, similar to R’s ```cmdscale()```

### Input:
You can provide either a binary matrix or a multi-sample VCF. Metadata should be formatted as described above.

### Running:
Run time will vary depending on the sizer of the dataset - memory limits have been placed to promote running on local machines (vs. HPCs) when necessary.

Example - running this tool using a binary matrix:
```
malaria-pipeline pca-plot \
  --matrix popgen_africa.mat.bin \
  --metadata metadata.tsv \
  --group-by country region year \
  --pcs 1,2 1,3 \
  --max-sample-missing 0.2 \
  --outdir pca_plots
```
Running using a filtered VCF file:

*Note: ```--max-sample-missing``` currently only works if you supply a ```--matrix```, not a VCF file*
```
malaria-pipeline pca-plot \
  --vcf popgen_africa.vcf.gz \
  --metadata metadata.tsv \
  --group-by country \
  --outdir pca_plots
```

### Output:
The command produces:
- PCA/PCoA plots (in PDF format) for each grouping variable (e.g. year, country, region, or unique column value)
- One plot per principal component (PC) pair (default: PC1–PC2, PC1–PC3), specific component pairs can be requested


## Trees

## Identity-by-Descent (IBD)
This toolkit supports an end-to-end IBD workflow using hmmIBD, starting from a binary SNP matrix, producing pairwise IBD fractions, windowed IBD summaries, gene annotation, and article-style plots.

#### Inputs
Input files should include:
1. A binary SNP matrix
- rows = SNPs
- columns = samples
- values of 0 (Reference), 1 (Alternate), 0.5 (mixed-call), N (missing, NA)

2. Metadata (.tsv) - tab-deliminated and must contain a "sample_id" column, a category (e.g. country, year) column, and an "fws" column
3. Reference genome index - the .fai file alongside the Pf3D7 reference genome
4. Gene product annotation - the "pf_genome_product_v3.tsv" file can be downloaded from https://plasmodb.org/plasmo/app/downloads but will also be available alongside this GitHub repo installation (Work-in-Progress)

Chromosome names are expected to be Pf3D7-style (e.g. Pf3D7_01_v3). By default, the toolkit removes:
- Pf3D7_API_v3
- Pf3D7_MIT_v3

### IBD Network plots

### Selection using IBD - Step 1: Run HmmIBD from a binary matrix
Runs hmmIBD per category (e.g. per country), and optionally per subgroup within a category (e.g. Country split by year or sampling site).

```
malaria-pipeline hmmibd-matrix \
  --matrix snps.mat.bin \
  --metadata metadata.tsv \
  --category-col country \
  --category Ethiopia \
  --subgroup-col year \
  --outdir hmmibd_ethiopia_by_year \
  --na-char N \
  --exclude-chr Pf3D7_API_v3,Pf3D7_MIT_v3
```
Key Options:
- ```--category-col```: metadata column defining your primary groups (e.g. country, region)
- ```--category```: run a subset (comma-separated); leave blank (or do not use) to run all categories (e.g. all countries)
- ```--subgroup-col```: optional secondary grouping within a category (e.g. year)
- ```--fws-th```: *Fws* threshold required for samples to be included; Default set to 0.95 to only include "monoclonal" samples
- ```--maf```: Minor allele frequency minimum reuirement, Default set to 0.01
- ```--hmmibd-bin```: path/name of hmmIBD executable if not on PATH
- ```--skip-hmmibd```: only write inputs, do not run hmmIBD

#### Outputs
Output files will be written into the ```--outdir```

One subdirectory per category/subgroup containing hmmIBD raw outputs:
- hmmIBD_<group>_maf0.01.txt
- hmmIBD_<group>_maf0.01_out.hmm.txt
- hmmIBD_<group>_maf0.01_out.hmm_fract.txt

A shared legend file:
- ibd_matrix_hap_leg.tsv

### Selection - Step 2: Summarise IBD into windows + gene annotation
Aggregates hmmIBD segments into genome-wide IBD fraction windows (default 50 kb - some users might find 10kb more appropriate) per group, and annotates high-IBD windows with known genes.

```
malaria-pipeline hmmibd-summary \
  --workdir hmmibd_ethiopia_by_year \
  --species Pf3D7 \
  --suffix 10_12_2025 \
  --window_size 50000 \
  --maf 0.01 \
  --quantile_cutoff 0.90 \
  --remove_chr Pf3D7_API_v3,Pf3D7_MIT_v3
```
Key Options:
- ```--suffix```: prefix for output files (use dates in DD_MM_YYYY format, like "10_12_2025" - **USE UNDERSCORES - NO SPACES**)
- ```--quantile_cutoff```: defines “high-IBD windows” for annotation (e.g. 0.90 or 0.95)

#### Outputs
Output files will be written into the ```--workdir```
- <suffix>_hmmIBD_ibd_win50kb.tsv
- <suffix>_hmmIBD_fraction.tsv
- <suffix>_hmmIBD_ibd_win50kb_annotated_q0.90.tsv
- <suffix>_hmmIBD_ibd_win50kb_annotated_q0.90_wide.tsv

### Selection - Step 1: Plot IBD summaries
Generates:
- pairwise fraction IBD boxplots
- genome-wide fraction plots (Manhattan-style per group)
- chromosome painting plots (high-IBD windows coloured by group, drug-resistance overlays in red)

```
malaria-pipeline hmmibd-ibdplots \
  --workdir hmmibd_ethiopia_by_year \
  --species Pf3D7 \
  --suffix 10_12_2025 \
  --window_size 50000 \
  --quantile_cutoff 0.90 \
  --remove_chr Pf3D7_API_v3,Pf3D7_MIT_v3
```

#### Outputs
Written, by default, into ```/workdir/win_50kb/```
- ibd_fraction_boxplot.png
- ibd_genomewide_fraction_cleaned.png
- ibd_chromosome_painting.png

## Selection calculating iHS and XP-EHH
The toolkit provides two complementary selection scan workflows using the *rehh* R-based package:
- *iHS* (within-population scan) directly from the binary SNP matrix
- *XP-EHH* (between-population comparisons) from scan_hh outputs (```scanned_haplotypes_<pop>.tsv```)
Both support category + subgroup style grouping (e.g. country and year), consistent with the IBD workflow.

### *iHS* Workflow
Runs a full iHS pipeline:
1. subset samples by category (and optional subgroup)
2. build haplotype objects per chromosome using ```rehh::data2haplohh```
3. compute ```scan_hh```
4. compute iHS using ```ihh2ihs```
5. output per-group plots + extreme-site tables with gene and drug-resistance highlighting

Example - Ethiopia split by year:
```
malaria-pipeline ihs-selection \
  --workdir selection_ihs_ethiopia_by_year \
  --matrix_binary snps.mat.bin \
  --metadata metadata.tsv \
  --species Pf3D7 \
  --label_category country \
  --subgroup_col year \
  --category Ethiopia \
  --focus-pop Ethiopia_2013 \
  --remove_chr Pf3D7_API_v3,Pf3D7_MIT_v3
```
### *XP-EHH* Workflow
*Under construction...*

### Fst
*Under construction...*

## Structural variant detection
*Under construction...*

## Admixture ancestry analysis
*Under construction...*
