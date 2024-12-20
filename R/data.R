# nolint start

#' Conversion table for Ensembl gene ids and gene symbols
#'
#'
#' @description This table includes the mapping for gene ids to gene symbols from different
#' reference genome gene annotation lists.
#' Included are the original gene symbols and the modified gene symbols that
#' are created when running the `make.unique()` function, as is done when
#' importing data using Seurat.
#'
#' @format
#' A data frame with 7 columns:
#' \describe{
#'   \item{gene_ids}{Ensembl gene ids}
#'   \item{gene_symbol_scpca}{The gene symbol used in the ScPCA reference}
#'   \item{gene_symbol_scpca_unique}{The gene symbol from the ScPCA reference, after `make.unique()`}
#'   \item{gene_symbol_10x2020}{The gene symbol used in the 2020 10x human genome reference}
#'   \item{gene_symbol_10x2020_unique}{The gene symbol from the 2020 10x human genome reference, after `make.unique()`}
#'   \item{gene_symbol_10x2024}{The gene symbol used in the 2024 10x human genome reference}
#'   \item{gene_symbol_10x2024_unique}{The gene symbol from the 2024 10x human genome reference, after `make.unique()`}
#' }
"scpca_gene_reference"
