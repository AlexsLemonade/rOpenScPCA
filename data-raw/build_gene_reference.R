#!/usr/env Rscript

# This script downloads and stores reference gene lists from 10x Genomics
# datasets, creating a final table of Ensembl ids and corresponding symbols
# including the symbols that would be created on read by Seurat by application
# of the make.unique() function.
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(dplyr)
})

## Download the reference data -------------------------------------------------

setwd(here::here())

# Read in the test data ScPCA SCE object and extract the row data:
genes_scpca <- readRDS(testthat::test_path("data", "scpca_sce.rds")) |>
  rowData() |>
  as.data.frame() |>
  # Use Ensembl ID if gene symbol is missing, then make unique
  mutate(
    gene_symbol_scpca = ifelse(is.na(gene_symbol), gene_ids, gene_symbol),
    gene_symbol_scpca_unique = make.unique(gene_symbol_scpca)
  ) |>
  select(gene_ids, gene_symbol_scpca, gene_symbol_scpca_unique)


# Download and read in a 2020 10x reference dataset and extract the gene symbols.
# Note that the 2020 Cell Ranger reference does not use Ensembl gene IDs for
# missing symbols, but the 2024 reference does.
url_10x2020 <- "https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5" # nolint
temp_10x2020 <- tempfile(fileext = ".h5")
download.file(url_10x2020, temp_10x2020, mode = "wb")
on.exit(unlink(temp_10x2020), add = TRUE) # delete when done

genes_10x2020 <- DropletUtils::read10xCounts(temp_10x2020) |>
  rowData() |>
  as.data.frame() |>
  filter(Type == "Gene Expression") |>
  rename(
    gene_ids = ID,
    gene_symbol_10x2020 = Symbol
  ) |>
  # add unique column
  mutate(gene_symbol_10x2020_unique = make.unique(gene_symbol_10x2020)) |>
  select(gene_ids, gene_symbol_10x2020, gene_symbol_10x2020_unique)

# Download and read in a 2024 10x reference dataset and extract the gene symbols.
url_10x2024 <- "https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x/5k_Human_Donor2_PBMC_3p_gem-x_5k_Human_Donor2_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5" # nolint
temp_10x2024 <- tempfile(fileext = ".h5")
download.file(url_10x2024, temp_10x2024, mode = "wb")
on.exit(unlink(temp_10x2024), add = TRUE) # delete when done


genes_10x2024 <- DropletUtils::read10xCounts(temp_10x2024) |>
  rowData() |>
  as.data.frame() |>
  filter(Type == "Gene Expression") |>
  rename(
    gene_ids = ID,
    gene_symbol_10x2024 = Symbol
  ) |>
  mutate(gene_symbol_10x2024_unique = make.unique(gene_symbol_10x2024)) |>
  select(gene_ids, gene_symbol_10x2024, gene_symbol_10x2024_unique)

# Join the gene lists ----------------------------------------------------------
scpca_gene_reference <- genes_scpca |>
  full_join(genes_10x2020, by = "gene_ids") |>
  full_join(genes_10x2024, by = "gene_ids")

## Add the table to package data -----------------------------------------------
usethis::use_data(
  scpca_gene_reference,
  version = 3,
  overwrite = TRUE,
  compress = "xz"
)
