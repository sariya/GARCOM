
# GARCOM

<!-- badges: start -->
<!-- badges: end -->

The goal of GARCOM is to provide mutation counts per individual within genetic boundaries (genes). It accepts different data formats with input file from plink (.raw), 
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

gene_annot_counts(recodedgen,snpgene) #input data: .raw formatted and SNP-gene (two columns)
```
## Citation

``` r
citation("GARCOM")
```

## Dependancies
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
