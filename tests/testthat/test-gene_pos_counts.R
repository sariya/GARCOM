context("check gene_pos_counts function")

test_that("function testing",{
              expect_error(gene_pos_counts(-4, 1,1))
	expect_is(gene_pos_counts(recodedgen,  snppos, genecoord),'data.table')
})


