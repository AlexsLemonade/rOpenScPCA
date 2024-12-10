#' Convert an SCE object to Seurat
#'
#' Converts an ScPCA SingleCellExperiment (SCE) object to Seurat format. This is
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
#' @param seurat_assay_version Whether to create Seurat `v3` or `v5` assay
#'   objects. Default is `v3`.
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
#' # convert SCE to Seurat with v5 assay objects
#' seurat_obj <- sce_to_seurat(sce, seurat_assay_version = "v5")
#' }
#'
sce_to_seurat <- function(
    sce,
    use_symbols = TRUE,
    reference = c("sce", "scpca", "10x2020", "10x2024"),
    dedup_method = c("unique", "sum"),
    seurat_assay_version = c("v3", "v5")) {
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
      reference != "sce" || all(c("gene_ids", "gene_symbol") %in% names(rowData(sce)))
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
  if (seurat_assay_version == "v5") {
    create_seurat_assay <- SeuratObject::CreateAssay5Object
    suppressWarnings(sobj[["RNA"]] <- as(sobj[["RNA"]], "Assay5"))
  } else {
    create_seurat_assay <- SeuratObject::CreateAssayObject
    suppressWarnings(sobj[["RNA"]] <- as(sobj[["RNA"]], "Assay"))
  }

  # add spliced data if present
  if ("spliced" %in% assayNames(sce)) {
    sobj[["spliced"]] <- create_seurat_assay(counts = assay(sce, "spliced"))
  }

  # add altExps as needed.
  for (alt_exp_name in names(alt_exps)) {
    alt_exp <- alt_exps[[alt_exp_name]]
    stopifnot(
      "All altExps must contain a `counts` assay." = "counts" %in% assayNames(alt_exp)
    )

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
