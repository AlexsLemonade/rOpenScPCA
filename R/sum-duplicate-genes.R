#' Sum counts for genes with duplicate names in a SingleCellExperiment object.
#'
#' Genes with the same name are merged by summing their raw expression counts.
#' If requested, the log-normalized expression values are recalculated, otherwise
#' this is left blank.
#'
#'
#' @param sce a SingleCellExperiment object with duplicated row names
#' @param normalize a logical indicating whether to normalize the expression
#'   values. Default is TRUE
#' @param recalculate_reduced_dims a logical indicating whether to recalculate
#'   PCA and UMAP. If FALSE, the input reduced dimensions are copied over. If
#'   TRUE, the highly variable genes are also recalculated with the new values
#'   stored in metadata. Default is FALSE
#'
#' @return a SingleCellExperiment object
#' @export
#'
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#'
sum_duplicate_genes <- function(sce, normalize = TRUE, recalculate_reduced_dims = FALSE) {
  stopifnot(
    "sce must be a SingleCellExperiment object" = is(sce, "SingleCellExperiment"),
    "normalize must be a logical" = is.logical(normalize),
    "recalculate_reduced_dims must be a logical" = is.logical(recalculate_reduced_dims),
    "normalize can not be FALSE if recalculate_reduced_dims is TRUE" = normalize || !recalculate_reduced_dims
  )
  if (normalize) {
    stopifnot(
      "Package `scran` must be installed if `normalize = TRUE` is set." =
        requireNamespace("scran", quietly = TRUE),
      "Package `scuttle` must be installed if `normalize = TRUE` is set." =
        requireNamespace("scuttle", quietly = TRUE)
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

  # calculate the reduced matrices
  counts <- rowsum(counts(sce), rownames(sce)) |> as("sparseMatrix")
  if ("spliced" %in% assayNames(sce)) {
    spliced <- rowsum(assay(sce, "spliced"), rownames(sce)) |> as("sparseMatrix")
    assays <- list(counts = counts, spliced = spliced)
  } else {
    assays <- list(counts = counts)
  }


  if (recalculate_reduced_dims) {
    reduced_dims <- list()
  } else {
    reduced_dims <- reducedDims(sce)
  }


  # Build the new SingleCellExperiment object
  summed_sce <- SingleCellExperiment(
    assays = assays,
    colData = colData(sce),
    metadata = metadata(sce),
    # if we are not recalculating reduced dimensions, copy over previous (likely similar)
    reducedDims = reduced_dims,
    altExps = altExps(sce)
  )

  # Add normalized values if requested
  if (normalize) {
    try({
      # try to cluster similar cells
      # clustering may fail if < 100 cells in dataset
      suppressWarnings({
        qclust <- scran::quickCluster(summed_sce)
        summed_sce <- scran::computeSumFactors(summed_sce, clusters = qclust)
      })
    })
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
