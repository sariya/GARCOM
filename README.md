
# GARCOM

<!-- badges: start -->
<!-- badges: end -->

The goal of GARCOM is to provide mutation counts per individual within genetic boundaries (genes). It accepts different data formats with input file from [plink](https://www.cog-genomics.org/plink/1.9/index) (.raw), 
gene boundaries, SNP location. 


## Installation

You can install the released version of GARCOM from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("GARCOM")
```

## Example

This is a small example which shows you how to use GARCOM:

``` r
library(GARCOM)
## basic example code
## sample data provided with library: genecoord, snpgene, snppos and genecoord

## Input data requires output from PLINK --recode flag. plink --bfile input --recode A --out sample_output 

#input data: .raw formatted and SNP-gene (two columns)
gene_annot_counts(recodedgen,snpgene) 

#input data: .raw formatted, SNP location (two columns) and Gene boundaries (three columns)
gene_pos_counts(recodedgen, snppos, genecoord) 

##For more examples refer manual
```
## Citation

``` r
citation("GARCOM")
```

## Dependencies
```
data.table
```

## suggests
```
testthat
```

## Issues and suggestions
```
GARCOM welcomes suggestions and improvements. For issues please open on the github
```

## Origin
```
GARCOM is derived from French word garçom (/ɡaʁ.sɔ̃/); 
here garcom is ready to serve to obtain desired results for the genetics data 
```