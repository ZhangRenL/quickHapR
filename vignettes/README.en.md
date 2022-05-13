---
Title:  quickHapR
Author: ZhangRenL
Date:   2022.05.14
---

# quickHapR

## 1 Introduction

The gene haplotype were generated through vcf file.

Haplotype analysis was performed based on the genotype information of a large number of populations obtained by next-generation sequencing. Combined with phenotypic data to screen out excellent haplotypes for follow-up research quickHapR

## 2. Basic logic

1.  Import user's datasets (vcf data, GFF annotations, phenotype data)
2.  Calculate haplotype through vcf file and generate haplotype data
3.  The haplotype data were then combined with the annotation information to mark the variant sites on the gene diagram; and show the evolutionary relationship between the haplotypes
4.  Compare phenotypic differences between different haplotypes and screen for dominant haplotypes

## 3. Installation Tutorial

**Notice:** To avoid unnecessary installation problems:

-   It is recommended to use the latest version of R(\>= 4.2.1)and Rtools4.2 for Windows 10 and above

-   It is recommended to use R(= 4.1.3)and Rtools4.0 for Windows7 system(due to system differences, the dependent Packages may not be installed normally)

### 3.1 Installation preparation

-   [ ] Install Rtools software
-   [ ] Install git software
-   [ ] Install R packages: `devtools`, `BiocManager`

### 3.2 quickHapR automatic installation method

``` r
# auto install command
# 1. Choose to install from gitee(domestic code hosting platform)
# Note: This method needs to install Rtools and Git first
devtools :: install_git("https://gitee.com/zhangrenl/quickhapr")

# 2. Choose to install from github
devtools :: install_github("zhangrenl/quickhapr")
```

### 3.3 quickHapR manual installation

``` r
# If the installation of the above command fails, you can go to Gitee to download the precompiled R package and select local installation
# Note the latest version was not promised by this way.
# Before installing quickHapR locally, you need to manually install the dependent R packages:
install.packages("BiocManager")
library(BiocManager)
install(c("ggpubr", "vcfR", "tidyverse", "stringr", "reshape2", "randomcoloR",
          "rtracklayer", "trackViewer", "GenomicRanges", "IRanges"))
```

## 4. How to use

### 4.1 Package Testing

``` r
library(quickHapR)
data("quickHap_test")
hap <- get_hap(vcf, hyb_remove = TRUE, na.drop = TRUE)
hapVsPheno(hap = hap, pheno = pheno, phenoName = "GrainWeight.2021", minAcc = 3)
results <- hapVsPheno(hap = hap,
                      pheno = pheno,
                      phenoName = "GrainWeight.2021",
                      minAcc = 3,
                      mergeFigs = TRUE)
plot(results$figs)
```

### 4.2 Software usage

``` r
# Load quickhapR
library(quickHapR)
data("quickHap_test")# Load test data, you don't have to execute this line when processing your own data

# set working directory
setwd("/your/working/directory")

# Import DataSets
vcf = import_vcf("YourVcfFile")
gff = import_gff("YourGffFile")
phenos = import_pheno("YourPhenoFile")

# Calculate and output haplotype results
# hap, data.frame: The first column and the last column are fixed as HAP and Accession respectively, and the middle column is the position and the corresponding genotype
# The first four lines of comment information are: CHR, POS, ALLELE, INFO
hap = get_hap(vcf,                 # import_vcf()imported vcfR
              filter_Chr = FALSE,  # filter chromosome options
              Chr = "scaffold_1",  # Screen the vcf information by chromosome
              filter_POS = FALSE,  # filter by position
              startPOS = 136756,   # Numeric, starting position, filter vcf information by position
              endPOS = 144094)     # Numeric, end position, filter vcf information by position

# hapResult, data.frame: The first column is fixed as HAP, the last two columns are fixed as Accession and freq respectively, the middle column is the position and the corresponding genotype
# The first four lines of comment information are: CHR, POS, ALLELE, INFO
hapResult = hap_result(hap,          # hap result
                       out  = FALSE, # Whether to output the file, if TRUE, the output path file must be specified
                       file = "results/Seita.1G001600_hapResult.txt")  # output file path(tab separated table)


# Visualize haplotype results
plotGeneStructure(gff,                # gff annotation information
                  hapResult,          # haplotype result
                  Chr = "scaffold_1", # the chromosome where the gene is located
                  startPOS = 136756,  # the starting position of the schematic diagram of the gene structure
                  endPOS = 144094,    # the end position of the schematic diagram of the gene structure
                  type = "pin",       # SNP type
                  cex = 1,            # circle size
                  CDS_h = 0.05,       # height of different gene structures
                  fiveUTR_h = 0.02,
                  threeUTR_h = 0.01)
plotHapTable(hapResult,               # haplotype result
             hapPrefix = "H",         # Haplotype prefix(letters before Arabic numerals)
             geneID = "",             # gene ID, as chart Title
             title.color = "grey90")  # header background color

# Association analysis of haplotype and phenotype
res = hapVsPheno(hap,        # data.frame: The first column and the last column are fixed as HAP and Accession respectively, and the middle column is the position and the corresponding genotype
                 phenos,     # data.frame: The first column is fixed as Accession, then each column is phenotype data, phenoName is used as colnames
                 phenoName = "yourPhenoName", # phenotype name used in this analysis
                 hapPrefix = "H",             # prefix of haplotype number
                 geneID = "Seita.1G000000",   # Gene ID, as header information
                 mergeFigs = TRUE,    # Whether to merge the two images
                 minAcc = 5)          # The minimum amount of data contained in the haplotype to be analyzed

# plot(res$fig_pvalue)
# plot(res$fig_Violin)
plot(res$figs)
```
