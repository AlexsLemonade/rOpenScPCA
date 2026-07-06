set.seed(2024)
sce <- splatter::simpleSimulate(nGenes = 1000, verbose = FALSE) |>
  scater::logNormCounts() |>
  scater::runPCA(ncomponents = 10)
test_mat <- reducedDim(sce, "PCA")


cluster_df <- calculate_clusters(test_mat)

test_that("calculate_silhouette works as expected", {
  df <- calculate_silhouette(test_mat, cluster_df)

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "silhouette_width", "silhouette_other")
  )
  expect_equal(df$cell_id, rownames(test_mat))
  expect_equal(df$cluster, cluster_df$cluster)
  expect_vector(df$silhouette_width, ptype = numeric())
  expect_s3_class(df$silhouette_other, "factor")
})


test_that("calculate_silhouette works as expected with non-default cluster column name", {
  cluster_df <- cluster_df |>
    dplyr::rename(clusters = cluster)
  df <- calculate_silhouette(test_mat, cluster_df, cluster_col = "clusters")

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "silhouette_width", "silhouette_other")
  )
  expect_equal(df$clusters, cluster_df$clusters)
})


test_that("calculate_silhouette works as expected with non-default cell id column name", {
  cluster_df <- cluster_df |>
    dplyr::rename(barcodes = cell_id)
  df <- calculate_silhouette(test_mat, cluster_df, cell_id_col = "barcodes")

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "silhouette_width", "silhouette_other")
  )
  expect_equal(df$cluster, cluster_df$cluster)
})


test_that("calculate_silhouette throws a warning when there is only 1 cluster", {
  cluster_df <- cluster_df |>
    dplyr::mutate(cluster = 1)

  expect_warning({
    df <- calculate_silhouette(test_mat, cluster_df)
  })

  expected_output <- cluster_df |>
    dplyr::mutate(
      silhouette_width = NA,
      silhouette_other = NA
    )
  expect_equal(expected_output, df)
})


test_that("calculate_purity works as expected", {
  df <- calculate_purity(test_mat, cluster_df)

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "purity", "maximum_neighbor")
  )
  expect_equal(df$cell_id, rownames(test_mat))
  expect_equal(df$cluster, cluster_df$cluster)
  expect_vector(df$purity, ptype = numeric())
  expect_s3_class(df$maximum_neighbor, "factor")
})


test_that("calculate_purity works as expected with non-default cluster column name", {
  cluster_df <- cluster_df |>
    dplyr::rename(clusters = cluster)
  df <- calculate_purity(test_mat, cluster_df, cluster_col = "clusters")

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "purity", "maximum_neighbor")
  )
  expect_equal(df$clusters, cluster_df$clusters)
})


test_that("calculate_purity works as expected with non-default cell id column name", {
  cluster_df <- cluster_df |>
    dplyr::rename(barcodes = cell_id)

  df <- calculate_purity(test_mat, cluster_df, cell_id_col = "barcodes")

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "purity", "maximum_neighbor")
  )
  expect_equal(df$cluster, cluster_df$cluster)
})


test_that("calculate_purity throws a warning when there is only 1 cluster", {
  cluster_df <- cluster_df |>
    dplyr::mutate(cluster = 1)

  expect_warning({
    df <- calculate_purity(test_mat, cluster_df)
  })

  expected_output <- cluster_df |>
    dplyr::mutate(
      purity = NA,
      maximum_neighbor = NA
    )
  expect_equal(expected_output, df)
})


test_that("calculate_stability works as expected with defaults", {
  df <- calculate_stability(test_mat, cluster_df)

  expected_names <- colnames(cluster_df)[!(colnames(cluster_df) %in% c("cell_id", "cluster"))]
  expect_setequal(
    colnames(df),
    c("replicate", "ari", expected_names)
  )
  expect_equal(df$replicate, 1:20) # checks rows too
  expect_vector(df$ari, ptype = numeric())
})


test_that("calculate_stability works as expected with different replicates", {
  df <- calculate_stability(test_mat, cluster_df, replicates = 2)

  expect_equal(nrow(df), 2)
})


test_that("calculate_stability works as expected with object and pc_name", {
  reducedDimNames(sce) <- "my_pca"

  df <- calculate_stability(
    sce,
    cluster_df,
    replicates = 2,
    pc_name = "my_pca"
  )
  expect_equal(nrow(df), 2)
})


test_that("calculate_stability warnings argument works", {
  expect_silent({
    df <- calculate_stability(test_mat, cluster_df)
  })

  # only run the next test on lower versions of bluster
  # 1.16.0 is the Bioconductor 3.20 version, in which warnings were turned off by default.
  # previous versions of Bioconductor will print a warning here
  skip_if(packageVersion("bluster") >= "1.16.0")

  # set 1 replicate, otherwise there will be 20 warnings
  expect_warning({
    df <- calculate_stability(test_mat, cluster_df, replicates = 1, seed = 1, warnings = TRUE)
  })
})


test_that("calculate_stability works as expected with non-default cluster column name", {
  cluster_df <- cluster_df |>
    dplyr::rename(clusters = cluster)

  df <- calculate_stability(test_mat, cluster_df, cluster_col = "clusters")

  expected_names <- colnames(cluster_df)[!(colnames(cluster_df) %in% c("cell_id", "clusters"))]
  expect_setequal(
    colnames(df),
    c("replicate", "ari", expected_names)
  )
})


test_that("calculate_stability works as expected with non-default cell id name", {
  cluster_df <- cluster_df |>
    dplyr::rename(barcodes = cell_id)

  df <- calculate_stability(test_mat, cluster_df, cell_id_col = "barcodes")

  expected_names <- colnames(cluster_df)[!(colnames(cluster_df) %in% c("barcodes", "cluster"))]
  expect_setequal(
    colnames(df),
    c("replicate", "ari", expected_names)
  )
})


test_that("calculate_stability errors as expected", {
  # cluster_df too short
  expect_error({
    calculate_stability(test_mat, cluster_df[1:5, ])
  })

  # cluster_df too long
  cluster_df_extra <- cluster_df |>
    tibble::add_row(cell_id = "extra_barcode")
  expect_error({
    calculate_stability(test_mat, cluster_df_extra)
  })
})
