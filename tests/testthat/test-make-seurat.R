test_that("SCE to Seurat with no id conversion works", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  seurat_obj <- sce_to_seurat(sce, use_symbols = FALSE)

  expect_equal(dim(seurat_obj), dim(sce))
  expect_setequal(rownames(seurat_obj), rownames(sce))
  expect_setequal(colnames(seurat_obj), colnames(sce))

  expect_equal(names(seurat_obj@assays), c("RNA", "spliced"))
  expect_equal(Seurat::DefaultAssay(seurat_obj), "RNA")

  # array contents
  expect_contains(
    slotNames(seurat_obj[["RNA"]]),
    c("counts", "data", "scale.data", "var.features", "meta.features")
  )

  expect_equal(seurat_obj[["RNA"]]$counts, counts(sce))
  expect_equal(seurat_obj[["RNA"]]$data, logcounts(sce))

  expect_setequal(
    seurat_obj[["RNA"]]@var.features,
    metadata(sce)$highly_variable_genes
  )

  expect_equal(names(seurat_obj@reductions), c("pca", "umap"))
  expect_equal(
    dim(seurat_obj[["pca"]]),
    dim(reducedDim(sce, "PCA"))
  )
  expect_equal(
    dim(seurat_obj[["umap"]]),
    dim(reducedDim(sce, "UMAP"))
  )
})


test_that("SCE to Seurat with id conversion works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  seurat_obj <- sce_to_seurat(sce, use_symbols = TRUE, reference = "sce")

  new_genes <- suppressMessages(
    ensembl_to_symbol(rownames(sce), unique = TRUE, sce = sce)
  )
  new_genes <- gsub("_", "-", new_genes)

  hv_genes <- suppressMessages(
    ensembl_to_symbol(metadata(sce)$highly_variable_genes, unique = TRUE, sce = sce)
  )
  hv_genes <- gsub("_", "-", hv_genes)

  expect_equal(dim(seurat_obj), dim(sce))
  expect_setequal(colnames(seurat_obj), colnames(sce))
  expect_setequal(rownames(seurat_obj), new_genes)

  expect_equal(names(seurat_obj@assays), c("RNA", "spliced"))
  expect_equal(Seurat::DefaultAssay(seurat_obj), "RNA")

  # array contents
  expect_contains(
    slotNames(seurat_obj[["RNA"]]),
    c("counts", "data", "scale.data", "var.features", "meta.features")
  )

  expect_setequal(
    seurat_obj[["RNA"]]@var.features,
    hv_genes
  )
  expect_equal(names(seurat_obj@reductions), c("pca", "umap"))
  expect_equal(names(seurat_obj@reductions), c("pca", "umap"))
  expect_equal(
    dim(seurat_obj[["pca"]]),
    dim(reducedDim(sce, "PCA"))
  )
  expect_equal(
    dim(seurat_obj[["umap"]]),
    dim(reducedDim(sce, "UMAP"))
  )
})

test_that("SCE to Seurat with id conversion and 10x reference works as expected", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  seurat_obj <- sce_to_seurat(sce, use_symbols = TRUE, reference = "10x2020")

  new_genes <- suppressMessages(
    ensembl_to_symbol(rownames(sce), unique = TRUE, reference = "10x2020")
  )
  new_genes <- gsub("_", "-", new_genes)

  hv_genes <- suppressMessages(
    ensembl_to_symbol(metadata(sce)$highly_variable_genes, unique = TRUE, reference = "10x2020")
  )
  hv_genes <- gsub("_", "-", hv_genes)

  expect_equal(dim(seurat_obj), dim(sce))
  expect_setequal(colnames(seurat_obj), colnames(sce))
  expect_setequal(rownames(seurat_obj), new_genes)

  expect_equal(names(seurat_obj@assays), c("RNA", "spliced"))
  expect_equal(Seurat::DefaultAssay(seurat_obj), "RNA")

  # array contents
  expect_contains(
    slotNames(seurat_obj[["RNA"]]),
    c("counts", "data", "scale.data", "var.features", "meta.features")
  )

  expect_setequal(
    seurat_obj[["RNA"]]@var.features,
    hv_genes
  )
  expect_equal(names(seurat_obj@reductions), c("pca", "umap"))
  expect_equal(names(seurat_obj@reductions), c("pca", "umap"))
  expect_equal(
    dim(seurat_obj[["pca"]]),
    dim(reducedDim(sce, "PCA"))
  )
  expect_equal(
    dim(seurat_obj[["umap"]]),
    dim(reducedDim(sce, "UMAP"))
  )
})

test_that("conversion works with altExps", {
  sce <- readRDS(test_path("data", "scpca_sce.rds"))
  altsce <- SingleCellExperiment(assays = list(counts = counts(sce), logcounts = logcounts(sce)))
  altExps(sce) <- list(
    alt1 = altsce,
    alt2 = altsce
  )

  seurat_obj <- sce_to_seurat(sce, use_symbols = FALSE)

  expect_equal(dim(seurat_obj), dim(sce))
  expect_setequal(rownames(seurat_obj), rownames(sce))
  expect_setequal(colnames(seurat_obj), colnames(sce))

  expect_setequal(names(seurat_obj@assays), c("RNA", "spliced", "alt1", "alt2"))
  expect_equal(Seurat::DefaultAssay(seurat_obj), "RNA")
})
