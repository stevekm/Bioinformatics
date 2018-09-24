# TCGA Tumor Analysis


# Setup & Run

Clone this repository

```
git clone ...
cd ...
```

Install Nextflow in the current directory

```
make install
```

Install `conda` in the current directory

```
make conda
make conda-install
```
- NOTE: `conda` installation may fail if current directory path length is >120 characters

Download the TCGA data required from the [GDC](https://portal.gdc.cancer.gov/repository)

```
make gdc
```

Run the pipeline

```
make run
```
- NOTE: by default the pipeline uses the `phoenix` Nextflow profile, configured for an HPC cluster running SGE

# Methods

Masked somatic mutation data (.maf) from whole exome experiments was retrieved from the [GDC repository](https://portal.gdc.cancer.gov/repository) for 10028 TCGA samples representing 33 tumor types (TCGA Projects). Variants were filtered, keeping only entries that matched the following criteria:

- Not listed as 'non_coding'
- Not lacking an HGVSp entry
- Lack COSMIC entries
- Lack 'dbSNP_RS' entires or were listed as 'novel'
- Lack ExAC_AF entries ('NA')

The sum of all remaining variants was calculated per sample. This value was divided by 38, representing the estimated number of Megabases covered across TCGA whole exome sequencing capture kits. [1] The resulting value represents the estimated tumor mutation burden for each sample.

## TCGA Projects Used

- ACC
- BLCA
- BRCA
- CESC
- CHOL
- COAD
- DLBC
- ESCA
- GBM
- HNSC
- KICH
- KIRC
- KIRP
- LAML
- LGG
- LIHC
- LUAD
- LUSC
- MESO
- OV
- PAAD
- PCPG
- PRAD
- READ
- SARC
- SKCM
- STAD
- TGCT
- THCA
- THYM
- UCEC
- UCS
- UVM

# Software Requirements

- Java 8 (Nextflow)

- Linux or macOS

# Resources & References

- 1. [Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0424-2). (https://doi.org/10.1186/s13073-017-0424-2)

- 2. [Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients](https://www.nature.com/articles/nm.4333). (http://dx.doi.org/10.1038/nm.4333)

- 3. [Mutational heterogeneity in cancer and the search for new cancer genes](https://www.nature.com/articles/nature12213). (http://dx.doi.org/10.1038/nature12213)
