context("check gene_pos_counts function")

test_noSNP<-data.frame("BP"=c(12111),"SNP"=c("SNP11")) ## no SNP in the raw data - throw an error as the bounadries are not there for the SNP

test_noSNPinraw<-data.frame("BP"=c(1500),"SNP"=c("SNP11")) ## no SNP in the raw data but present in the boundaries

test_SNP_nosum<-data.frame("BP"=c(27000),"SNP"=c("SNP10"))  ##return null as no person carries allele of this SNP

##test with raw data that has SNP with 0 allele count
##test with raw data that has SNP but not in the boundaries

test_gen_withnoSNP<-data.frame("FID"=c(paste("FID",seq(1:10),sep="")),"IID"=c(paste("IID_sample",seq=(1:10),sep="")),
"PAT"=c(rep(0,10)),"MAT"=c(rep(0,10)),"SEX"=c( c( c(rep(1,5)),rep(0,5))),"PHENOTYPE"=c( c(rep("NA",2),rep("1",4),
rep("0",4))), "SNP12_A"=c(rep(0,10))) ## this has no SNP in the annotation

test_gen_with_oneSNP<-data.frame("FID"=c(paste("FID",seq(1:10),sep="")),"IID"=c(paste("IID_sample",seq=(1:10),sep="")),
"PAT"=c(rep(0,10)),"MAT"=c(rep(0,10)),"SEX"=c( c( c(rep(1,5)),rep(0,5))),"PHENOTYPE"=c( c(rep("NA",2),rep("1",4),
rep("0",4))), "SNP1_A"=c(rep(1,1),rep(0,9)) ,"SNP12_A"=c(rep(0,10))) ## this has one SNP in the annotation

test_genecoord_noSNP<-data.frame("GENE"=c("GENEX"), "START"=c(180000), "END"=c(190000)) #No snp is present in these boundaries

test_that("function testing",{
    expect_error(gene_pos_counts(-4, 1,1)) #empty and trash data to test around
    expect_is(gene_pos_counts(recodedgen,snppos,genecoord),'data.table') #returns data.table dataset. Check it
    expect_error(gene_pos_counts(recodedgen,test_noSNP,genecoord)) 
    expect_null(gene_pos_counts(recodedgen,test_SNP_nosum,genecoord))
    
    expect_error(gene_pos_counts(test_gen_withnoSNP,snppos,genecoord))
    expect_equal(nrow(gene_pos_counts(test_gen_with_oneSNP, snppos,genecoord)), 1) ##returns one row 
    expect_error(gene_pos_counts(recodedgen, snppos,test_genecoord_noSNP)) ## there are no SNPs test_genecoord_noSNP in gene boundaries for this    
})
##
test_that("function testing for subsetting",{

    expect_error(gene_pos_counts(recodedgen, snppos,genecoord, keep_indiv=c("IID121"))) ## subset IIDs
    expect_error(gene_pos_counts(recodedgen, snppos,genecoord, filter_gene=c("gene111"))) ## subset genes
    expect_is(gene_pos_counts(recodedgen, snppos,genecoord, filter_gene=c("gene111","GENE1")),'data.table')  ## subset genes, with one valid gene
    expect_error(gene_pos_counts(recodedgen, snppos,genecoord, extract_SNP=c("snps111"))) ## subset snps
    expect_null(gene_pos_counts(recodedgen, snppos,genecoord,extract_SNP=c("snps111","SNP10"))) ## subset snps. SNP10 has zero counts
    expect_is(gene_pos_counts(recodedgen, snppos,genecoord,extract_SNP=c("snps111","SNP1")),'data.table') ## subset snps. SNP1 has various counts
    
    expect_is(gene_pos_counts(recodedgen, snppos,genecoord,keep_indiv=c("IID121","IID_sample1")),'data.table')  ## subset IIDs, with one valid IIDs

    expect_is(gene_pos_counts(recodedgen, snppos,genecoord,keep_indiv=c("IID_sample1","IID_sample4"),filter_gene=c("GENE1")),'data.table') ##subset and further subset: IIDs and gene filters
    expect_equal(nrow(gene_pos_counts(recodedgen, snppos,genecoord,keep_indiv=c("IID_sample1","IID_sample4"),filter_gene=c("GENE1"))),1) ## we check row counts

    expect_null(gene_pos_counts(recodedgen, snppos,genecoord,extract_SNP=c("SNP1","SNP3","SNP5","SNP6"),keep_indiv=c("IID_sample2"))) ## sample2 has zero allele for these SNPs
})

test_that("function testing for impute",{

    expect_is(gene_pos_counts(recodedgen,snppos,genecoord,impute_missing=TRUE,impute_method="median"),'data.table') ##test for median
    expect_error(gene_pos_counts(recodedgen,snppos,genecoord,impute_missing=TRUE,impute_method="Mean")) ##make sure that spellings are OK
})
##testing ends


