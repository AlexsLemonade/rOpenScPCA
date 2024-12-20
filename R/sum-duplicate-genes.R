#' Sum counts for genes with duplicate names in a SingleCellExperiment object.
#'
#' @description Genes with the same name are merged by summing their raw expression counts.
#' When multiple Ensembl gene IDs are associated with the same gene symbol,
#' identifier conversion can result in duplicate gene names. This function
#' resolves such duplicates by summing the expression values for each duplicate
#' gene name, which may be justified if the different Ensembl gene IDs share
#' substantial sequence identity, which could make separate quantification of
#' the two genes less reliable.
#'
#' The rowData for the summed SingleCellExperiment object is updated to reflect
#' the new set of gene names. In each case, the first row for any duplicated id
#' is retained. This may mean that for gene symbols that correspond to multiple
#' Ensembl ids, the first Ensembl id is retained and the others are dropped.
#'
#' If requested, the log-normalized expression values are recalculated,
#' otherwise that matrix is left blank.
#'
#' PCA and UMAP are also recalculated if requested (which requires recalculating
#' the normalized expression). If not recalculating, the original reduced
#' dimensions are retained as they are likely to remain similar in most cases.
#'
#'
#' @param sce a SingleCellExperiment object with duplicated row names
#' @param normalize a logical indicating whether to normalize the expression
#'   values. Default is TRUE.
#' @param recalculate_reduced_dims a logical indicating whether to recalculate
#'   PCA and UMAP. If FALSE, the input reduced dimensions are copied over. If
#'   TRUE, the highly variable genes are also recalculated with the new values
#'   stored in metadata. Default is FALSE.
#' @param cell_set_size Because very large matrices can cause trouble, we break
#'  the summing into submatrices with this many cells. Default is 5000, which
#'  should be fine for most datasets.
#'
#' @return a SingleCellExperiment object
#' @export
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
#' @examples
#' \dontrun{
#' # sum across duplicated gene names
#' summed_sce <- sum_duplicate_genes(sce)
#'
#' # sum across duplicated gene names and recalculate PCA and UMAP
#' summed_sce <- sum_duplicate_genes(sce, recalculate_reduced_dims = TRUE)
#' }
#'
sum_duplicate_genes <- function(sce,
                                normalize = TRUE,
                                recalculate_reduced_dims = FALSE,
                                cell_set_size = 5000) {
  stopifnot(
    "sce must be a SingleCellExperiment object" = is(sce, "SingleCellExperiment"),
    "normalize must be a logical" = is.logical(normalize),
    "recalculate_reduced_dims must be a logical" = is.logical(recalculate_reduced_dims),
    "normalize can not be FALSE if recalculate_reduced_dims is TRUE" = normalize || !recalculate_reduced_dims,
    "cell_set_size should be an integer (much) greater than 1" = cell_set_size > 1 && cell_set_size %% 1 == 0
  )
  if (normalize) {
    stopifnot(
      "Package `scran` must be installed if `normalize = TRUE` is set." =
        requireNamespace("scran", quietly = TRUE)
    )
  }
  stopifnot(
    "Package `scater` must be installed if `recalculate_reduced_dims = TRUE` is set." =
      requireNamespace("scater", quietly = TRUE)
  )

  if (!any(duplicated(rownames(sce)))) {
    message("No duplicated gene names found in the SingleCellExperiment object, returning original object")
    return(sce)
  }

  # calculate the reduced matrices----------------------------------------------
  unique_rows <- unique(rownames(sce)) # new row names

  # break up counts and assay matrices to avoid large non-sparse matrices during calculation
  cell_sets <- seq(0, ncol(sce) - 1) %/% cell_set_size
  # make sure the last set is not length one
  cell_sets[ncol(sce)] <- cell_sets[ncol(sce) - 1]

  counts <- unique(cell_sets) |>
    purrr::map(\(x) {
      rowsum(counts(sce)[, cell_sets == x], rownames(sce)) |>
        as("sparseMatrix") |>
        base::`[`(unique_rows, ) # reorder rows (faster with sparse)
    }) |>
    purrr::reduce(cbind)

  if ("spliced" %in% assayNames(sce)) {
    spliced <- unique(cell_sets) |>
      purrr::map(\(x) {
        rowsum(assay(sce, "spliced")[, cell_sets == x], rownames(sce)) |>
          as("sparseMatrix") |>
          base::`[`(unique_rows, )
      }) |>
      purrr::reduce(cbind)
    assays <- list(counts = counts, spliced = spliced)
  } else {
    assays <- list(counts = counts)
  }

  # regenerate rowData, using first row for each duplicate
  row_data <- rowData(sce)[unique_rows, ]

  if (recalculate_reduced_dims) {
    reduced_dims <- list()
  } else {
    reduced_dims <- reducedDims(sce)
  }



  # Build the new SingleCellExperiment object
  summed_sce <- SingleCellExperiment(
    assays = assays,
    rowData = row_data,
    colData = colData(sce),
    metadata = metadata(sce),
    # if we are not recalculating reduced dimensions, copy over previous (likely similar)
    reducedDims = reduced_dims,
    altExps = altExps(sce)
  )
  # remove and replace existing Feature stats
  rowData(summed_sce)$mean <- NULL
  rowData(summed_sce)$detected <- NULL
  summed_sce <- scuttle::addPerFeatureQCMetrics(summed_sce)

  # Add normalized values if requested
  if (normalize) {
    try(
      {
        # try to cluster similar cells
        # clustering may fail if < 100 cells in dataset
        suppressWarnings({
          qclust <- scran::quickCluster(summed_sce)
          summed_sce <- scran::computeSumFactors(summed_sce, clusters = qclust)
        })
      },
      silent = TRUE
    )
    summed_sce <- scuttle::logNormCounts(summed_sce)
  }

  # recalculate PCA and UMAP if requested, using the same dimensions as before (if available)
  if (recalculate_reduced_dims) {
    if ("PCA" %in% reducedDimNames(sce)) {
      pca_dim <- ncol(reducedDim(sce, "PCA"))
    } else {
      pca_dim <- min(50, ncol(sce) - 1) # always reduce dimensions at least 1!
    }

    if (is.null(metadata(sce)$highly_variable_genes)) {
      n_hvg <- 2000 # same default as ScPCA
    } else {
      n_hvg <- length(metadata(sce)$highly_variable_genes)
    }
    hv_genes <- scran::modelGeneVar(summed_sce) |>
      scran::getTopHVGs(n = n_hvg)
    metadata(summed_sce)$highly_variable_genes <- hv_genes # store the new HVGs

    summed_sce <- scater::runPCA(summed_sce, subset_row = hv_genes, ncomponents = pca_dim) |>
      scater::runUMAP()
  }

  return(summed_sce)
}
