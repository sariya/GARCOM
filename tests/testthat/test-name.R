context("data integrity")

test_gene_coord<-data.table::data.table(GENE=c("ABC","XYG","alpha"),"START"=c(10,200,320),"END"=c(101,250,350))
test_snp_pos<-data.table::data.table(SNP=c("SNP1","SNP2","SNP3"),"BP"=c(101,250,350))
test_snp_gene<-data.table::data.table(SNP=c("SNP1","SNP2","SNP3"),"GENE"=c("ABC","BRCA1","gamma"))

test_that("data types correct works", {

expect_is(test_gene_coord$START, 'numeric')
expect_is(test_gene_coord$END, 'numeric')
expect_is(test_snp_pos$BP, 'numeric')
})

test_that("column names works", {

 expect_named(test_gene_coord, c("GENE","START","END"))
 expect_named(test_snp_pos, c("SNP","BP"))
 expect_named(test_snp_gene, c("SNP","GENE"))

})

## https://katherinemwood.github.io/post/testthat/

##https://b-rodrigues.github.io/fput/unit-testing.html