#' @title recodedgen
#'
#' @description sample genetic data
#'
#' @details recodedgen is sample genetic data provided with the GARCOM package. It has with 10 rows and 16 columns.
#' Where the first 6 columns are FID, IID, PAT, MAT, SEX and PHENOTYPE which are inherent from the 
#' PLINK (recode) output. Columns followed by PHENOTYPE are SNP names which are suffixed with _A, or _C
#' or _T or _G. Each SNP column may have 0, 1, 2 or NA value. Where NA represents missingness.
#'

recodedgen<-data.frame("FID"=c(paste("FID",seq(1:10),sep="")),
"IID"=c(paste("IID_sample",seq=(1:10),sep="")),
"PAT"=c(rep(0,10)),"MAT"=c(rep(0,10)),
"SEX"=c( c( c(rep(1,5)),rep(0,5))),
"PHENOTYPE"=c( c(rep("NA",2),rep("1",4),
rep("0",4))),"SNP1_A"=c(rep(1,1),rep(0,9)),
"SNP2_T"=c(rep(1,2),rep(0,7),rep(NA,1)),
"SNP3_G"=c(rep(0,2),rep(1,8)),
"SNP4_C"=c(rep(NA,2),rep(0,7),rep(1,1)),
"SNP5_C"=c(rep(NA,1),rep(0,8),rep(1,1)),
"SNP6_G"=c(rep(1,10)),"SNP7_G"=c(rep(0,10)),
"SNP8_C"=c(rep(0,10)), "SNP9_T"=c(rep(0,10)), "SNP10_A"=c(rep(0,10)))

#' @title snpgene
#' @description sample SNP-gene annotation data
#'
#' @details snpgene is sample SNP-Gene annotation data provided with the GARCOM package. It has 10 rows and 2 columns. 
#' Column names are GENE and SNP, where GENE column contains GENE names and SNP column
#' contains SNP name that is annotated with the GENE
#'

snpgene<-data.frame("GENE"=
c("GENE1","GENE2","GENE3","GENE4",
"GENE5","GENE1","GENE1","GENE3","GENE2","GENE4"),
"SNP"=c("SNP1","SNP2","SNP3",
"SNP4","SNP5","SNP6","SNP7","SNP8","SNP9","SNP1"))

#' @title snppos
#'
#' @description sample data for SNP coordinates
#'
#' @details snppos is sample SNP-BP data provided with the GARCOM package. It has 2 columns and 10 rows. Column names are SNP and BP. 
#' where SNP columns contains SNP names BP column contains position of the SNP

snppos<-data.frame("SNP"=c(paste("SNP",seq(1:10),sep="")), "BP"=c(1100, 89200, 2500, 33000, 5500, 69500, 12000,8800, 23200, 27000))

#' @title genecoord 
#'
#' @description sample data for gene base pair boundaries
#'
#' @details genecoord is example data provided with the GARCOM package. It has 3 columns and 5 rows. 
#' Column names are GENE, START and END where GENE column contains gene name, START and END indicate 
#' start BP and end BP respectively.

genecoord<-data.frame("GENE"=c("GENE1","GENE2","GENE3","GENE4","GENE5"), 
"START"=c(1000,2100,5000,40000,23000), 
"END"=c(2000,3000,9000,45000,30000))
