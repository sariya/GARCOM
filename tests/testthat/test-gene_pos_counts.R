context("check gene_pos_counts function")

test_noSNP<-data.frame("BP"=c(12111),"SNP"=c("SNP11")) ## no SNP in the raw data - throw an error as the bounadries are not there for the SNP

test_noSNPinraw<-data.frame("BP"=c(1500),"SNP"=c("SNP11")) ## no SNP in the raw data but present in the boundaries

test_SNP_nosum<-data.frame("BP"=c(27000),"SNP"=c("SNP10"))  ##return null as no person carries allele of this SNP

##test with raw data that has SNP with 0 allele count
##test with raw data that has SNP but not in the boundaries

test_that("function testing",{
              expect_error(gene_pos_counts(-4, 1,1))
	expect_is(gene_pos_counts(recodedgen,  snppos, genecoord),'data.table')
	expect_error(gene_pos_counts(recodedgen, test_noSNP,genecoord))
expect_null(gene_pos_counts(recodedgen, test_SNP_nosum,genecoord))
})


