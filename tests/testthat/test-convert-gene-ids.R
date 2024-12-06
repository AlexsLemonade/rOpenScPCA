test_that("basic gene symbol conversion works", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")

  gene_symbols <- ensembl_to_symbol(ensembl_ids)
  expect_equal(gene_symbols, c("TP53", "MYCN"))
})

test_that("gene symbol conversion works with unexpected ids", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323", "foobar")

  expect_message(gene_symbols <- ensembl_to_symbol(ensembl_ids))
  expect_equal(gene_symbols, c("TP53", "MYCN", "foobar"))

  expect_no_warning(gene_symbols_na <- ensembl_to_symbol(ensembl_ids, leave_na = TRUE))
  expect_equal(gene_symbols_na, c("TP53", "MYCN", NA))
})

test_that("gene symbol conversion works for 10x references", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2020")
  expect_equal(gene_symbols, c("TP53", "MYCN"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2024")
  expect_equal(gene_symbols, c("TP53", "MYCN"))
})

test_that("gene symbol conversion works for 'unique' gene symbols", {
  ensembl_ids <- c("ENSG00000015479", "ENSG00000269226")

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "scpca", unique = FALSE)
  expect_equal(gene_symbols, c("MATR3", "TMSB15B"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "scpca", unique = TRUE)
  expect_equal(gene_symbols, c("MATR3.1", "TMSB15B.1"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2020", unique = FALSE)
  expect_equal(gene_symbols, c("MATR3", "TMSB15B"))

  gene_symbols <- ensembl_to_symbol(ensembl_ids, reference = "10x2020", unique = TRUE)
  expect_equal(gene_symbols, c("MATR3.1", "TMSB15B.1"))
})

test_that("gene symbol conversion works using an SCE reference", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  ensembl_ids <- c("ENSG00000015479", "ENSG00000269226")

  expect_message(gene_symbols <- ensembl_to_symbol(ensembl_ids, sce = sce, unique = FALSE))
  expect_equal(gene_symbols, c("MATR3", "TMSB15B"))

  expect_message(gene_symbols <- ensembl_to_symbol(ensembl_ids, sce = sce, unique = TRUE))
  expect_equal(gene_symbols, c("MATR3.1", "TMSB15B.1"))
})

test_that("gene symbol conversion in seurat compatibility mode works", {
  ensembl_ids <- c("ENSG00000141510", "ENSG00000134323")
  gene_symbols <- ensembl_to_symbol(ensembl_ids, unique = FALSE, seurat_compatible = TRUE)
  expect_equal(gene_symbols, c("TP53", "MYCN"))

  ensembl_ids <- c("ENSG00000285609", "ENSG00000252254", "ENSG00000283274")
  gene_symbols <- ensembl_to_symbol(ensembl_ids, unique = FALSE, seurat_compatible = TRUE)
  expect_equal(gene_symbols, c("5S-rRNA", "Y-RNA", "5-8S-rRNA"))
})



test_that("conversion of a full sce object works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  gene_symbols <- rowData(sce)$gene_symbol
  names(gene_symbols) <- rowData(sce)$gene_ids
  gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]

  expect_message(converted_sce <- sce_to_symbols(sce)) |> expect_message() # two messages expected here
  expect_equal(rownames(converted_sce), unname(gene_symbols))

  # check that hvg and PCA were converted too.
  expected_hvg <- gene_symbols[metadata(sce)$highly_variable_genes]
  expect_equal(metadata(converted_sce)$highly_variable_genes, expected_hvg)

  rotation_ids <- rownames(attr(reducedDim(converted_sce, "PCA"), "rotation"))
  expected_rotation_ids <- gene_symbols[rownames(attr(reducedDim(sce, "PCA"), "rotation"))]
  expect_equal(rotation_ids, expected_rotation_ids)
})


test_that("conversion of an sce object using a reference works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  gene_symbols <- rowData(sce)$gene_symbol
  names(gene_symbols) <- rowData(sce)$gene_ids
  gene_symbols[is.na(gene_symbols)] <- names(gene_symbols)[is.na(gene_symbols)]

  # testing with the ScPCA reference, which should be the same as internal table
  converted_sce <- sce_to_symbols(sce, reference = "scpca")
  expect_equal(rownames(converted_sce), unname(gene_symbols))

  # check that hvg and PCA were converted too.
  expected_hvg <- gene_symbols[metadata(sce)$highly_variable_genes]
  expect_equal(metadata(converted_sce)$highly_variable_genes, expected_hvg)

  rotation_ids <- rownames(attr(reducedDim(converted_sce, "PCA"), "rotation"))
  expected_rotation_ids <- gene_symbols[rownames(attr(reducedDim(sce, "PCA"), "rotation"))]
  expect_equal(rotation_ids, expected_rotation_ids)
})
