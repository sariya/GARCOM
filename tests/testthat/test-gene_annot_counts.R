context("check gene_annot_counts function")

testexp_snpgene<-data.table::data.table("GENE"=c("GENE1"),"SNP"=c("SNP10")) ## 0 values in each person for this SNP
testexp_snpgene_noSNP<-data.table::data.table("GENE"=c("GENE1"),"SNP"=c("SNP11")) ## no SNP in the raw data
testexp_snpgene_pass<-data.table::data.table("GENE"=c("GENE1"),"SNP"=c("SNP4")) ## we'll return on row

test_gen_withnoSNP<-data.frame("FID"=c(paste("FID",seq(1:10),sep="")),"IID"=c(paste("IID_sample",seq=(1:10),sep="")),
"PAT"=c(rep(0,10)),"MAT"=c(rep(0,10)),"SEX"=c( c( c(rep(1,5)),rep(0,5))),"PHENOTYPE"=c( c(rep("NA",2),rep("1",4),
rep("0",4))), "SNP12_A"=c(rep(0,10))) ## this has no SNP in the annotation

test_gen_with_oneSNP<-data.frame("FID"=c(paste("FID",seq(1:10),sep="")),"IID"=c(paste("IID_sample",seq=(1:10),sep="")),
"PAT"=c(rep(0,10)),"MAT"=c(rep(0,10)),"SEX"=c( c( c(rep(1,5)),rep(0,5))),"PHENOTYPE"=c( c(rep("NA",2),rep("1",4),
rep("0",4))), "SNP1_A"=c(rep(1,1),rep(0,9)) ,"SNP12_A"=c(rep(0,10))) ## this has one SNP in the annotation

test_that("function testing",{
    expect_error(gene_annot_counts(-4, 1))
    expect_is(gene_annot_counts(recodedgen,snpgene),'data.table') ## recieve data.table 
    expect_null(gene_annot_counts(recodedgen,testexp_snpgene)) #testexp_snpgene NULL returned for this data set
    expect_null(gene_annot_counts(recodedgen,testexp_snpgene_noSNP)) 
    
    expect_equal(nrow(gene_annot_counts(recodedgen,testexp_snpgene_pass)), 1) #check nrow
    
    expect_null(gene_annot_counts(test_gen_withnoSNP,snpgene)) #genetic data have no SNP in the SNP-gene set
    expect_equal(nrow(gene_annot_counts(test_gen_with_oneSNP,snpgene)), 2)
})

test_that("function testing for subsetting",{

    expect_error(gene_annot_counts(recodedgen,snpgene,extract_SNP=c(""))) ##empty SNP names
    expect_error(gene_annot_counts(recodedgen,snpgene,extract_SNP=c("snp1"))) ##empty SNP names

    expect_error(gene_annot_counts(recodedgen,snpgene,keep_indiv=c(""))) ##empty iid
    expect_error(gene_annot_counts(recodedgen,snpgene,keep_indiv=c("snp1"))) ##wrong iid
    expect_error(gene_annot_counts(recodedgen,snpgene,filter_gene=c("snp1"))) ##wrong gene name

    expect_is(gene_annot_counts(recodedgen,snpgene,filter_gene=c("GENE1")),'data.table')
    expect_is(gene_annot_counts(recodedgen,snpgene,keep_indiv=c("IID_sample8")),'data.table') ##just one sample id
    
    expect_null(gene_annot_counts(recodedgen,snpgene,keep_indiv=c("IID_sample8"),filter_gene=c("GENE5"))) ## we provide IIDs, gene and subset. Return is null for this dataset
})
##testing ends

test_that("function testing for impute",{

    expect_is(gene_annot_counts(recodedgen,snpgene,impute_missing=TRUE,impute_method="median"),'data.table') ##test for median
    expect_equal( nrow(gene_annot_counts(recodedgen,snpgene,impute_missing=TRUE,filter_gene=c("GENE2","GENE1"),impute_method="median")),2) ##subset genes and impute using median. check number of rows as 2
 
})
##testing ends
