context("check gene_annot_counts function")

test_that("function testing",{
              expect_error(gene_annot_counts(-4, 1))
	expect_is(gene_annot_counts(recodedgen, snpgene),'data.table')
})
