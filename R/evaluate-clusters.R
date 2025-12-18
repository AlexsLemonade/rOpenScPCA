#' Calculate the silhouette width of clusters
#'
#' This function uses the `bluster::approxSilhouette()` function to calculate the
#' silhouette width for a clustering result. These results can be used downstream to
#' calculate the average silhouette width, a popular metric for cluster evaluation.
#'
#' @param x Either a matrix of principal components (PCs), or a SingleCellExperiment
#'  or Seurat object containing PCs. If a matrix is provided, rows should be cells
#'  and columns should be PCs, and row names should be cell ids (e.g., barcodes).
#' @param cluster_df A data frame that contains at least two columns: one representing
#'   unique cell ids, and one containing cluster assignments. By default, these columns
#'   should be named `cell_id` and `cluster` respectively, though this can be customized.
#'   The cell id column's values should match either the PC matrix row names, or the
#'   SingleCellExperiment/Seurat object cell ids. Typically this data frame will be
#'   output from the `rOpenScPCA::calculate_clusters()` function.
#' @param cluster_col The name of the column in `cluster_df` which contains cluster
#'   assignments. The default is "cluster".
#' @param cell_id_col The name of the column in `cluster_df` which contains unique cell
#'   ids. The default is "cell_id".
#' @param pc_name Optionally, the name of the PC matrix in the object. Not used if a
#'   matrix is provided. If the name is not provided, the name "PCA" is used for
#'   SingleCellExperiment objects, and "pca" for Seurat objects.
#'
#' @return Expanded `cluster_df` data frame with additional columns `silhouette_width`,
#'   the cell's silhouette width, and `silhouette_other`, the closest cluster other
#'   than the one to which the given cell was assigned. For more information,
#'   see documentation for `bluster::approxSilhouette()`.
#'   If there is only one cluster in the provided data, silhouette width can not be
#'   calculated. The returned data frame will have `NA` values in the
#'   `silhouette_width` and `silhouette_other` columns.
#'
#' @importFrom stats setNames
#'
#' @export
#' @examples
#' \dontrun{
#' # calculate silhouette width for clusters stored in a data frame
#' cluster_df <- calculate_silhouette(sce_object, cluster_df)
#' }
calculate_silhouette <- function(
    x,
    cluster_df,
    cluster_col = "cluster",
    cell_id_col = "cell_id",
    pc_name = NULL) {
  expected_df_names <- c(cell_id_col, cluster_col)
  stopifnot(
    "The cell id column name must be length of 1." = length(cell_id_col) == 1,
    "The cluster column name must be length of 1." = length(cluster_col) == 1,
    "Expected columns not present in cluster_df." =
      all(expected_df_names %in% colnames(cluster_df))
  )

  if (length(unique(cluster_df[[cluster_col]])) == 1) {
    warning("There is only 1 cluster in this data. Silhouette width will not be calculated.")
    cluster_df <- cluster_df |>
      dplyr::mutate(
        silhouette_width = NA,
        silhouette_other = NA
      )
    return(cluster_df)
  }


  silhouette_df <- x |>
    prepare_pc_matrix(pc_name) |>
    bluster::approxSilhouette(cluster_df[[cluster_col]]) |>
    as.data.frame() |>
    # note this gets renamed later as needed
    tibble::rownames_to_column("cell_id") |>
    dplyr::rename(
      "silhouette_width" = "width",
      "silhouette_other" = "other"
    )

  # join with cluster_df in this direction, so that columns in cluster_df come first
  # ensure provided cluster_df column names are used as well
  silhouette_df <- cluster_df |>
    dplyr::inner_join(silhouette_df, by = setNames(c("cell_id", "cluster"), c(cell_id_col, cluster_col)))

  return(silhouette_df)
}




#' Calculate the neighborhood purity of clusters
#'
#' This function uses the `bluster::neighborPurity()` function to calculate the
#' neighborhood purity values for a clustering result.
#'
#' @param x Either a matrix of principal components (PCs), or a SingleCellExperiment
#'  or Seurat object containing PCs. If a matrix is provided, rows should be cells
#'  and columns should be PCs, and row names should be cell ids (e.g., barcodes).
#' @param cluster_df A data frame that contains at least two columns: one representing
#'   unique cell ids, and one containing cluster assignments. By default, these columns
#'   should be named `cell_id` and `cluster` respectively, though this can be customized.
#'   The cell id column's values should match either the PC matrix row names, or the
#'   SingleCellExperiment/Seurat object cell ids. Typically this data frame will be
#'   output from the `rOpenScPCA::calculate_clusters()` function.
#' @param cluster_col The name of the column in `cluster_df` which contains cluster
#'   assignments. The default is "cluster".
#' @param cell_id_col The name of the column in `cluster_df` which contains unique cell
#'   ids. The default is "cell_id".
#' @param pc_name Optionally, the name of the PC matrix in the object. Not used if a
#'   matrix is provided. If the name is not provided, the name "PCA" is used for
#'   SingleCellExperiment objects, and "pca" for Seurat objects.
#' @param ... Additional arguments to pass to `bluster::neighborPurity()`
#'
#' @return Expanded `cluster_df` data frame with the additional columns `purity`,
#'   the cell's neighborhood purity, and `maximum_neighbor`, the cluster with the
#'   highest proportion of observations neighboring the given cell. For more
#'   information see documentation for `bluster::neighborPurity()`.
#'   If there is only one cluster in the provided data, neighborhood purity can not be
#'   calculated. The returned data frame will have `NA` values in the `purity` and
#'   `maximum_neighbor` columns.
#'
#' @export
#' @examples
#' \dontrun{
#' # calculate neighborhood purity for clusters stored in a data frame
#' cluster_df <- calculate_purity(sce_object, cluster_df)
#' }
calculate_purity <- function(
    x,
    cluster_df,
    cluster_col = "cluster",
    cell_id_col = "cell_id",
    pc_name = NULL,
    ...) {
  expected_df_names <- c(cell_id_col, cluster_col)
  stopifnot(
    "The cell id column name must be length of 1." = length(cell_id_col) == 1,
    "The cluster column name must be length of 1." = length(cluster_col) == 1,
    "Expected columns not present in cluster_df." =
      all(expected_df_names %in% colnames(cluster_df))
  )

  if (length(unique(cluster_df[[cluster_col]])) == 1) {
    warning("There is only 1 cluster in this data. Neighborhood purity will not be calculated.")
    cluster_df <- cluster_df |>
      dplyr::mutate(
        purity = NA,
        maximum_neighbor = NA
      )
    return(cluster_df)
  }

  purity_df <- x |>
    prepare_pc_matrix(pc_name) |>
    bluster::neighborPurity(cluster_df[[cluster_col]], ...) |>
    as.data.frame() |>
    tibble::rownames_to_column(cell_id_col) |>
    dplyr::rename(
      "maximum_neighbor" = "maximum"
    )

  # join with cluster_df in this direction, so that columns in
  # cluster_df come first
  purity_df <- cluster_df |>
    dplyr::inner_join(purity_df, by = cell_id_col)

  return(purity_df)
}



#' Calculate cluster stability using the Adjusted Rand Index (ARI)
#'
#' This function generates and clusters, using provided parameters, bootstrap
#' replicates calculates the Adjusted Rand Index (ARI) between each set of bootstrapped
#' clusters and the original provided clusters. ARI measures similarity between different
#' cluster results, where a value of 0 indicates an entirely random relationship between
#' results, and a value of 1 indicates perfect agreement.
#'
#' When assessing stability, you should specify the same clustering parameters here as
#' were used to calculate the original clusters.
#'
#' Note that this function will also make use of `bluster::clusterRows()` with the
#' `bluster::NNGraphParam()` function on a principal components matrix. Note that defaults
#' for some arguments may differ from the `bluster::NNGraphParam()` defaults.
#' Specifically, the clustering algorithm defaults to "louvain" and the weighting scheme
#' to "jaccard" to align with common practice in scRNA-seq analysis.
#'
#'
#' @param x An object containing PCs that clusters were calculated from. This can be
#'   either a SingleCellExperiment object, a Seurat object, or a matrix where columns
#'   are PCs and rows are cells. If a matrix is provided, it must have row names of cell
#'   ids (e.g., barcodes).
#' @param cluster_df A data frame that contains at least two columns: one representing
#'   unique cell ids, and one containing cluster assignments. By default, these columns
#'   should be named `cell_id` and `cluster` respectively, though this can be customized.
#'   The cell id column's values should match either the PC matrix row names, or the
#'   SingleCellExperiment/Seurat object cell ids. Typically this data frame will be
#'   output from the `rOpenScPCA::calculate_clusters()` function.
#' @param cluster_col The name of the column in `cluster_df` which contains cluster
#'   assignments. The default is "cluster".
#' @param cell_id_col The name of the column in `cluster_df` which contains unique cell
#'   ids. The default is "cell_id".
#' @param replicates Number of bootstrap replicates to perform. The default is 20.
#' @param seed Random seed
#' @param pc_name Optionally, the name of the PC matrix in the object. Not used if a
#'   matrix is provided. If the name is not provided, the name "PCA" is used for
#'   SingleCellExperiment objects, and "pca" for Seurat objects.
#' @param warnings Whether warnings related to distance ties when calculating bootstrap
#'   clusters should be printed. The default is FALSE.
#' @param ... Additional arguments to pass to `calculate_clusters()` which calculates
#'   bootstrapped clusters. Usually, these will be the same arguments used to generate
#'   the original clusters.
#'
#' @return Data frame with columns `replicate` and `ari`, representing the given bootstrap replicate
#'   and its ARI value, respectively, and columns representing clustering algorithm parameters which
#'   include at least `algorithm`, `weighting`, and `nn`. Louvain and leiden clustering will also
#'   include `resolution`, and leiden clustering will further include `objective_function`.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # First, cluster PCs from a SingleCellExperiment object using default parameters
#' # and setting a seed for reproducibility
#' cluster_df <- calculate_clusters(sce_object, seed = 11)
#' # Second, calculate cluster stability using default parameters
#' stability_df <- calculate_stability(sce_object, cluster_df, seed = 11)
#'
#'
#' # First, cluster PCs from a SingleCellExperiment object using default parameters
#' # and setting a seed for reproducibility
#' cluster_df <- calculate_clusters(sce_object, seed = 11)
#' # Second, calculate cluster stability using default parameters and 50 replicates
#' stability_df <- calculate_stability(
#'   sce_object,
#'   cluster_df,
#'   replicates = 50,
#'   seed = 11
#' )
#'
#'
#' # First, cluster PCs from a SingleCellExperiment object using the leiden
#' # algorithm and a resolution of 0.1
#' cluster_df <- calculate_clusters(
#'   sce_object,
#'   algorithm = "leiden",
#'   resolution = 0.1,
#'   seed = 11
#' )
#' # Second, calculate cluster stability using the same parameters as were used
#' # for the initial clustering
#' stability_df <- calculate_stability(
#'   sce_object,
#'   cluster_df,
#'   algorithm = "leiden",
#'   resolution = 0.1,
#'   seed = 11
#' )
#' }
calculate_stability <- function(
    x,
    cluster_df,
    cluster_col = "cluster",
    cell_id_col = "cell_id",
    replicates = 20,
    seed = NULL,
    pc_name = NULL,
    warnings = FALSE,
    ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ensure we have a matrix
  pca_matrix <- prepare_pc_matrix(x, pc_name = pc_name)

  # ensure pca matrix and cluster df compatibility
  stopifnot(
    "The cluster dataframe must have the same number of rows as the PCA matrix." =
      nrow(pca_matrix) == nrow(cluster_df),
    "Cell ids in the cluster dataframe must match the PCA matrix rownames." =
      length(setdiff(rownames(pca_matrix), cluster_df[[cell_id_col]])) == 0
  )

  # Check columns
  expected_df_names <- c(cell_id_col, cluster_col)
  stopifnot(
    "The cell id column name must be length of 1." = length(cell_id_col) == 1,
    "The cluster column name must be length of 1." = length(cluster_col) == 1,
    "Expected columns not present in cluster_df." =
      all(expected_df_names %in% colnames(cluster_df))
  )

  # Extract vector of clusters, ensuring same order as pca_matrix
  rownames(cluster_df) <- cluster_df[[cell_id_col]]
  clusters <- cluster_df[rownames(pca_matrix), cluster_col]

  # calculate ARI for each cluster result bootstrap replicate
  all_ari_df <- 1:replicates |>
    purrr::map(
      \(i) {
        sample_cells <- sample(nrow(pca_matrix), replace = TRUE)
        resampled_pca <- pca_matrix[sample_cells, ]
        original_clusters <- clusters[sample_cells]

        resampled_df <- withCallingHandlers(
          calculate_clusters(resampled_pca, ...),
          warning = \(w) {
            if (!warnings) tryInvokeRestart("muffleWarning")
          }
        )

        ari <- pdfCluster::adj.rand.index(resampled_df$cluster, original_clusters)

        # return df with ari and clustering parameters
        ari_df <- resampled_df |>
          dplyr::slice(1) |>
          # these column names come directly out of calculate_clusters; they are not customized
          dplyr::select(!dplyr::all_of(c("cell_id", "cluster"))) |>
          dplyr::mutate(
            # define this variable here to ensure it's numeric
            replicate = i,
            ari = ari,
            # ensure these columns come first
            .before = "algorithm"
          )

        return(ari_df)
      }
    ) |>
    dplyr::bind_rows()


  return(all_ari_df)
}

#' Evaluate cluster results
#'
#' This wrapper function can be used to evaluate clusters from a single clustering
#' calculation or a list calculated using `sweep_clusters()` function.
#' Input should be be a data frame from a single `calculate_clusters()` call
#' or a list of data frames with the resulting clusters from all parameter
#' combinations provided to the `sweep_clusters()` function.
#' Output is a list of results or a single data.frame depending on what was provided.
#' Evaluation statistics are added in the form of columns to the original cluster
#' object provided.
#'
#' @param x An object containing PCs that clusters were calculated from. This can be
#'   either a SingleCellExperiment object, a Seurat object, or a matrix where columns
#'   are PCs and rows are cells. If a matrix is provided, it must have row names of cell
#'   ids (e.g., barcodes).
#' @param cluster_results A single data frame or list of data frames obtained from
#'  `rOpenScPCA::calculate_clusters()` or `rOpenScPCA::sweep_clusters()` respectively.
#'   Each data frame in the list should contain
#'   at least two columns: one representing unique cell ids, and one containing
#'   cluster assignments. By default, these columns should be named `cell_id` and
#'   `cluster` respectively, though this can be customized. The cell id column's
#'   values should match either the PC matrix row names, or the
#'   SingleCellExperiment/Seurat object cell ids.
#' @param metrics Which metrics should be collected? Options are one or both of "purity" or "silhouette".
#' Default is to collect both purity and silhouette.
#'
#' @return An updated list of data frames with additional columns from running evaluation with
#' `rOpenScPCA::calculate_silhouette()` and/or `rOpenScPCA::calculate_purity()`, based on
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Calculate Principal Components
#' pca_matrix <- reducedDim(sce_object, "PCA")
#'
#' # We can calculate a single data frame of cluster results:
#' cluster_df <- calculate_clusters(
#'   pca_matrix,
#'   algorithm = "leiden",
#'   resolution = 0.1,
#'   seed = 11
#' )
#'
#' # Then we can evaluate these cluster stats with
#' # calculate_cell_cluster_metrics:
#' sweep_list_evaled <- calculate_cell_cluster_metrics(
#'   x = pca_matrix,
#'   cluster_results = cluster_df
#' )
#'
#' ######### Evaluate a list of multiple cluster results ############
#' # If we obtain a list of clusters like so...
#' sweep_list <- sweep_clusters(
#'   sce_object,
#'   weighting = "jaccard",
#'   nn = c(10, 15, 25),
#'   resolution = c(0.75, 1),
#'   seed = 9
#' )
#'
#' # Then we can evaluate these cluster stats with
#' # calculate_cell_cluster_metrics:
#' sweep_list_evaled <- calculate_cell_cluster_metrics(
#'   x = pca_matrix,
#'   cluster_results = sweep_list
#' )
#'
#'
#'
#' # If we only want one or the other metric calculated we can use the metrics
#' # argument which can be either 'silhouette' or 'purity':
#' sweep_list_evaled <- calculate_cell_cluster_metrics(
#'   x = pca_matrix,
#'   cluster_results = sweep_list,
#'   metrics = "purity"
#' )
#' }
#'
calculate_cell_cluster_metrics <- function(x,
                                           cluster_results,
                                           metrics = c("purity", "silhouette")) {
  supported_evals <- c("purity", "silhouette")

  if (is.data.frame(cluster_results)) {
    cluster_results <- list(cluster_results)
  }
  # Check input arguments
  stopifnot(
    "`cluster_results` must be a data.frame or list containing data.frame objects" =
      is.list(cluster_results) && all(sapply(cluster_results, is.data.frame)),
    " Cluster `evals` that are supported are only 'purity' and 'silhouette'" =
      all(metrics %in% supported_evals)
  )

  evaled_list <- cluster_results |>
    purrr::map(
      \(df) {
        if ("purity" %in% metrics) {
          df <- calculate_purity(
            x = x,
            cluster_df = df
          )
        }
        if ("silhouette" %in% metrics) {
          df <- calculate_silhouette(
            x = x,
            cluster_df = df
          )
        }
        return(df)
      }
    )

  # if there's only one data.frame input then let's just output a data frame
  if (length(evaled_list) == 1) {
    result <- evaled_list[[1]]
  } else {
    result <- evaled_list
  }

  return(result)
}
