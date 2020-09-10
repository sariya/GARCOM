
# GARCOM

<!-- badges: start -->
<!-- badges: end -->

The goal of GARCOM is to provide mutation counts per individual within genetic boundaries (genes). It accepts different data formats with input file from [plink](https://www.cog-genomics.org/plink/1.9/index) (.raw), gene boundaries, SNP location. 
It also accepts [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file format. The vcf file is assumed to be pVCF, that is where all samples are merged/combined.  


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

#read VCF file vcf_data <- vcfR::read.vcfR("CHRXX.vcf.gz", verbose=TRUE)
vcf_counts_annot(vcf_data,df_snpgene) # pass vcf data read and data frame with SNP-gene annotation

#read VCF file vcf_data <- vcfR::read.vcfR("CHRXX.vcf.gz", verbose=TRUE)
vcf_counts_SNP_genecoords(vcf_data,df_snppos,df_genecoords) # pass vcf data read and data frame SNP position and third with gene coordinates

#subset individuals 
vcf_counts_SNP_genecoords(vcf_data,df_snptestpos, df_genecoordstestpos,keep_indiv=c("IID1","IID2"))

ind_select<-c("IID1","IID2") ## store in a vector
vcf_counts_SNP_genecoords(vcf_data,df_snptestpos, df_genecoordstestpos,keep_indiv=ind_select)

##For more examples refer manual
```
## Citation

``` r
citation("GARCOM")
```

## Dependencies (Imports)
```
data.table(>=v1.12.8)
stats
vcfR(>=v1.12.0)
```

## suggests
```
testthat
```

## Issues and suggestions
```
GARCOM welcomes suggestions and improvements. Please open issues on the github for bugs/suggestions.
```

## Origin
```
GARCOM is derived from French word garçom (/ɡaʁ.sɔ̃/); 
here garcom is ready to serve to obtain desired results for the genetics data 
```
## Version
```
Currently v1.1.0
```