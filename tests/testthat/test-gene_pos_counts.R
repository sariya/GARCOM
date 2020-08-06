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
              expect_error(gene_pos_counts(-4, 1,1))
	expect_is(gene_pos_counts(recodedgen,  snppos, genecoord),'data.table')
	expect_error(gene_pos_counts(recodedgen, test_noSNP,genecoord))
expect_null(gene_pos_counts(recodedgen, test_SNP_nosum,genecoord))

expect_error(gene_pos_counts(test_gen_withnoSNP, snppos,genecoord))

expect_equal(nrow(gene_pos_counts(test_gen_with_oneSNP, snppos,genecoord)), 1) ##returns one row 
expect_error(gene_pos_counts(recodedgen, snppos,test_genecoord_noSNP))

})


