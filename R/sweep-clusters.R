#' Calculate clusters across a set of parameters
#'
#' This function can be used to perform reproducible clustering while varying a set of parameters.
#' Multiple values can be provided for any of:
#'  - The algorithm (`algorithm`)
#'  - The weighting scheme (`weighting`)
#'  - Number of nearest neighbors (`nn`)
#'  - The resolution parameter (`resolution`)
#'  - The objective function parameter (`objective_function`)
#'
#' For each algorithm specified, all parameters possible to use with that
#' algorithm will be systematically varied. This function does not accept additional
#' parameters besides those listed above.
#' Note that defaults for some arguments may differ from the `bluster::NNGraphParam()` defaults.
#' Specifically, the clustering algorithm defaults to "louvain" and the weighting scheme to "jaccard"
#' to align with common practice in scRNA-seq analysis.
#'
#' @param x An object containing PCs that clustering can be performed in. This can be either
#'   a SingleCellExperiment object, a Seurat object, or a matrix where columns are PCs and
#'   rows are cells. If a matrix is provided, it must have row names of cell ids (e.g., barcodes).
#' @param algorithm Clustering algorithm to use. Must be one of "louvain" (default), "walktrap",
#'   or "leiden".
#' @param weighting Weighting scheme(s) to consider when sweeping parameters.
#' Provide a vector of unique values to vary this parameter. Options include "jaccard" (default),
#'   "rank", or "number"
#' @param nn Number of nearest neighbors to consider when sweeping parameters.
#'  Provide a vector of unique values to vary this parameter. The default is 10.
#' @param resolution Resolution parameter used by Louvain and Leiden clustering only.
#'   Provide a vector of unique values to vary this parameter. The default is 1.
#' @param objective_function Leiden-specific parameter for whether to use the
#'   Constant Potts Model ("CPM"; default) or "modularity". Provide a vector of unique values
#'   to vary this parameter.
#' @param seed Random seed to set for clustering.
#' @param threads Number of threads to use. The default is 1.
#' @param pc_name Name of principal components slot in provided object. This argument is only used
#'   if a SingleCellExperiment or Seurat object is provided. If not provided, the SingleCellExperiment
#'   object name will default to "PCA" and the Seurat object name will default to "pca".
#'
#' @return A list of data frames from performing clustering across all parameter combinations.
#'   Columns include `cluster_set` (identifier column for results from a single clustering run),
#'   `cell_id`, and `cluster`. Additional columns represent algorithm parameters and include at least:
#'   `algorithm`, `weighting`, and `nn`. Louvain and Leiden clustering will also include `resolution`,
#'   and Leiden clustering will further include `objective_function`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # perform Louvain clustering with Jaccard weighting (defaults),
#' # varying the nearest neighobor parameter, and set a seed for reproducibility
#' cluster_df <- sweep_clusters(
#'   sce_object,
#'   nn = c(10, 15, 20, 25),
#'   seed = 11
#' )
#'
#' # perform Louvain clustering, with Jaccard and rank weighting, and
#' # varying the nearest neighbor and resolution parameters.
#' cluster_df <- sweep_clusters(
#'   sce_object,
#'   algorithm = "louvain",
#'   weighting = c("jaccard", "rank"),
#'   nn = c(10, 15, 20, 25),
#'   resolution = c(0.5, 1),
#'   seed = 11
#' )
#'
#' # perform walktrap and Louvain clustering with Jaccard weighting, and
#' # varying the nearest neighbors for both algorithms, and resolution for Louvain.
#' cluster_df <- sweep_clusters(
#'   sce_object,
#'   algorithm = c("walktrap", "louvain"),
#'   weighting = "jaccard",
#'   nn = c(10, 15, 20, 25),
#'   resolution = c(0.5, 1),
#'   seed = 11
#' )
#' }
sweep_clusters <- function(
    x,
    algorithm = "louvain",
    weighting = "jaccard",
    nn = 10,
    resolution = 1, # Louvain or Leiden
    objective_function = "CPM", # Leiden only
    threads = 1,
    seed = NULL,
    pc_name = NULL) {
  # check and prepare matrix
  pca_matrix <- prepare_pc_matrix(x, pc_name = pc_name)

  # Collect all specific inputs into a single list
  sweep_params <- tidyr::expand_grid(
    algorithm = unique(algorithm),
    weighting = unique(weighting),
    nn = unique(nn),
    resolution = unique(resolution),
    objective_function = unique(objective_function)
  ) |>
    # set unused parameters for each algorithm to default; this will allow duplicates to be removed by distinct()
    dplyr::mutate(
      resolution = ifelse(algorithm %in% c("louvain", "leiden"), resolution, 1),
      objective_function = ifelse(algorithm == "leiden", objective_function, "CPM")
    ) |>
    dplyr::distinct()

  sweep_results <- sweep_params |>
    purrr::pmap(
      \(algorithm, weighting, nn, resolution, objective_function) {
        calculate_clusters(
          pca_matrix,
          algorithm = algorithm,
          weighting = weighting,
          nn = nn,
          resolution = resolution,
          objective_function = objective_function,
          threads = threads,
          seed = seed
        )
      }
    )

  return(sweep_results)
}
