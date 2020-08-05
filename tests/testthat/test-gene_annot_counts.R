context("check gene_annot_counts function")

testexp_snpgene<-data.table("GENE"=c("GENE1"),"SNP"=c("SNP10")) ## 0 values in each person for this SNP
testexp_snpgene_noSNP<-data.table("GENE"=c("GENE1"),"SNP"=c("SNP11")) ## no SNP in the raw data
testexp_snpgene_pass<-data.table("GENE"=c("GENE1"),"SNP"=c("SNP4")) ## we'll return on row

##  gene_annot_counts(recodedgen, testexp_snpgene_noSNP)
test_that("function testing",{
              expect_error(gene_annot_counts(-4, 1))
	expect_is(gene_annot_counts(recodedgen, snpgene),'data.table')
expect_null(gene_annot_counts(recodedgen, testexp_snpgene))
expect_null(gene_annot_counts(recodedgen, testexp_snpgene_noSNP))

 expect_equal(nrow(gene_annot_counts(recodedgen, testexp_snpgene_pass)), 1)
})
