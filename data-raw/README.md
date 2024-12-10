This directory contains scripts to download and preprocess data used within the rOpenScPCA package.
There are also notebooks that explore some of the prepared datasets for exploration.

## Test data

- The `download_test_data.R` script downloads the test data used in this package from the OpenScPCA test data bucket.


## Building gene references

- The `build_gene_references.R` script creates a table of Ensembl id to gene symbol references.
The initial table of gene ids and gene symbols is extracted from an example ScPCA-formatted SCE object (`rOpenScPCA/tests/testthat/data/scpca_sce.rds`).
This is combined with the reference information extracted from example 10x Genomics datasets.
The full table is saved in `data/scpca_gene_reference.rda` (overwriting any previous file).

- The `explore_gene_references.Rmd` notebook explores the resulting gene references table a bit to see where some of the conversions differ among different references.
