test_that("merging works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  double_sce <- rbind(sce, sce)
  deduped_sce <- sum_duplicate_genes(double_sce)

  expect_equal(dim(deduped_sce), dim(sce))

  expect_setequal(rownames(deduped_sce), rownames(sce))


  expect_equal(
    counts(deduped_sce)[rownames(sce), ], # need to make the row order identical
    2 * counts(sce) |> Matrix::drop0() # original matrices may have explicit 0s
  )
  expect_equal(
    assay(deduped_sce, "spliced")[rownames(sce), ],
    2 * assay(sce, "spliced") |> Matrix::drop0()
  )
})
