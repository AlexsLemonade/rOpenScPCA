#' Convert an SCE object to Seurat
#'
#' @description Converts an ScPCA SingleCellExperiment (SCE) object to Seurat format. This is
#' primarily a wrapper around Seurat::as.Seurat() with some additional steps to
#' include ScPCA metadata and options for converting the feature index from
#' Ensembl gene ids to gene symbols.
#'
#' If present, reduced dimensions from SCE objects will be included, renamed to
#' match Seurat default naming.
#'
#'
#' @param sce The input SCE object
#' @param use_symbols Logical indicating whether the Seurat object uses gene
#'   symbols for indexing. Default is TRUE.
#' @param reference If `use_symbols` is TRUE, the reference to use for gene
#'   symbols. One of `sce`, `scpca`, `10x2020`, `10x2024`, where `sce` uses the
#'   symbols stored in the `gene_symbol` column of the SCE object's row data,
#'   and the others use standardized translation tables. See `ensembl_to_symbol()`
#'   for more information. Default is `sce`.
#' @param dedup_method Method to handle duplicated gene symbols. If `unique`,
#'   the gene symbols will be made unique following standard Seurat procedures.
#'   If `sum`, the expression values for each gene symbols will be summed.
#' @param seurat_assay_version Whether to create Seurat `v5` or `v3` assay
#'   objects. Default is `v5`.
#'
#' @return a Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' # convert an ScPCA SCE to Seurat renaming with gene symbols
#' seurat_obj <- sce_to_seurat(sce)
#'
#' # convert SCE to Seurat keeping the current names (usually Ensembl)
#' seurat_obj <- sce_to_seurat(sce, use_symbols = FALSE)
#'
#' # convert SCE to Seurat, merging duplicated gene names
#' seurat_obj <- sce_to_seurat(sce, dedup_method = "sum")
#'
#' # convert SCE to Seurat with v3 assay objects
#' seurat_obj <- sce_to_seurat(sce, seurat_assay_version = "v3")
#' }
#'
sce_to_seurat <- function(
    sce,
    use_symbols = TRUE,
    reference = c("sce", "scpca", "10x2020", "10x2024"),
    dedup_method = c("unique", "sum"),
    seurat_assay_version = c("v5", "v3")) {
  reference <- match.arg(reference)
  dedup_method <- match.arg(dedup_method)
  seurat_assay_version <- match.arg(seurat_assay_version)

  stopifnot(
    "Package `Seurat` must be installed for SCE to Seurat conversion." =
      requireNamespace("Seurat", quietly = TRUE),
    "Package `SeuratObject` must be installed for SCE to Seurat conversion." =
      requireNamespace("SeuratObject", quietly = TRUE),
    "`sce` must be a SingleCellExperiment object." = is(sce, "SingleCellExperiment"),
    "`sce` must contain both `gene_ids` and `gene_symbol` columns in the row data if it is being used as a reference." =
      reference != "sce" || all(c("gene_ids", "gene_symbol") %in% names(rowData(sce))),
    "`sce` should only include a single library. Merged SingleCellExperiment objects are not (yet) supported." =
      length(metadata(sce)$library_id) <= 1
  )

  if (use_symbols) {
    sce <- suppressMessages(sce_to_symbols(
      sce,
      reference = reference,
      unique = (dedup_method == "unique"),
      seurat_compatible = TRUE
    ))
  }

  if (dedup_method == "sum") {
    sce <- sum_duplicate_genes(sce)
  }

  mainExpName(sce) <- "RNA"

  # convert reducedDims to Seurat standard naming
  reducedDimNames(sce) <- tolower(reducedDimNames(sce))
  reducedDims(sce) <- reducedDims(sce) |>
    purrr::map(\(mat) {
      # insert underscore between letters and numbers in axis labels
      rd_names <- gsub("([A-Za-z]+)([0-9]+)", "\\1_\\2", colnames(mat)) |>
        tolower()
      colnames(mat) <- rd_names
      if (!is.null(attr(mat, "rotation"))) {
        colnames(attr(mat, "rotation")) <- rd_names
      }
      return(mat)
    })

  # pull out any altExps to handle manually
  alt_exps <- altExps(sce)
  altExps(sce) <- NULL

  # Let Seurat do initial conversion to capture most things
  if ("logcounts" %in% assayNames(sce)) {
    data_name <- "logcounts"
  } else {
    data_name <- NULL
  }
  sobj <- Seurat::as.Seurat(sce, data = data_name)

  # update identity values, using cluster if present
  sce_meta <- metadata(sce)
  sobj$orig.ident <- sce_meta$library_id
  if (is.null(sce$cluster)) {
    Seurat::Idents(sobj) <- sce_meta$library_id
  } else {
    Seurat::Idents(sobj) <- sce$cluster
  }

  # add variable features if present
  if (!is.null(sce_meta$highly_variable_genes)) {
    sobj[["RNA"]]@var.features <- sce_meta$highly_variable_genes
    sce_meta$highly_variable_genes <- NULL # remove from metadata to avoid redundancy
  }

  # convert and set functions depending on assay version requested
  if (seurat_assay_version == "v3") {
    create_seurat_assay <- SeuratObject::CreateAssayObject
    suppressWarnings(sobj[["RNA"]] <- as(sobj[["RNA"]], "Assay"))
  } else {
    create_seurat_assay <- SeuratObject::CreateAssay5Object
    suppressWarnings(sobj[["RNA"]] <- as(sobj[["RNA"]], "Assay5"))
  }

  # add spliced data if present
  if ("spliced" %in% assayNames(sce)) {
    sobj[["spliced"]] <- create_seurat_assay(counts = assay(sce, "spliced"))
  }

  # add altExps as needed.
  for (alt_exp_name in names(alt_exps)) {
    alt_exp <- alt_exps[[alt_exp_name]]

    if (!"counts" %in% assayNames(alt_exp)) {
      warning(
        "The altExp `", alt_exp_name,
        "` does not have a `counts` assay, so it will be skipped."
      )
      next
    }

    if (nrow(alt_exp) == 1 && seurat_assay_version == "v5") {
      warning(
        "The altExp `", alt_exp_name, "` has only one feature;",
        " for Seurat v5 compatibility, a dummy feature will be added with all zero counts."
      )
      alt_exp <- add_dummy_feature(alt_exp)
    }

    # check name compatibility for Seurat
    alt_exp_rownames <- rownames(alt_exp)
    if (any(grepl("_", alt_exp_rownames))) {
      warning(
        "Replacing underscores ('_') with dashes ('-') in feature names from ",
        alt_exp_name,
        " for Seurat compatibility."
      )
      rownames(alt_exp) <- gsub("_", "-", alt_exp_rownames)
    }

    sobj[[alt_exp_name]] <- create_seurat_assay(counts = counts(alt_exp))
    if ("logcounts" %in% assayNames(alt_exp)) {
      sobj[[alt_exp_name]]$data <- logcounts(alt_exp)
    }
  }

  # add sample-level metadata after removing non-vector objects
  sobj@misc <- sce_meta |>
    purrr::discard(is.object) |>
    purrr::discard(is.list)

  # scale data because many Seurat functions expect this
  sobj <- Seurat::ScaleData(sobj, assay = "RNA", verbose = FALSE)

  return(sobj)
}


#' Add a dummy feature to a SingleCellExperiment object with all zeros
#'
#' @param sce A SingleCellExperiment object. Should contain only one row/feature.
#'   If larger, it will be returned unmodified.
#' @param feature_name The name to give to the new dummy feature
#'
#' @returns A SingleCellExperiment object with a dummy feature added (all zeros)
add_dummy_feature <- function(sce, feature_name = "null-feature") {
  # add a dummy feature an SCE object
  if (nrow(sce) > 1) {
    message(
      "No need for a dummy feature in a SingleCellExperiment object",
      " with more than one row. Returning original object."
    )
    return(sce)
  }
  if (feature_name %in% rownames(sce)) {
    stop("A feature named `", feature_name, "` already exists in the SingleCellExperiment object.")
  }

  # duplicate the existing row, rename
  sce <- sce[c(1, 1), ]
  rownames(sce) <- c(rownames(sce)[1], feature_name)
  # set all rowData to NA, all assays to 0
  rowData(sce)[feature_name, ] <- NA
  assays(sce) <- assays(sce) |> purrr::map(\(x) {
    x[feature_name, ] <- 0
    return(x)
  })
  return(sce)
}
