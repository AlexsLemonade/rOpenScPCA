#!/usr/env Rscript

# Downloads a test data file from OpenScPCA and places it in the test directory

setwd(here::here())
url_base <- "https://openscpca-test-data-release-public-access.s3.us-east-2.amazonaws.com/test/SCPCP000001/SCPCS000001/" # nolint

test_processed_path <- testthat::test_path("data", "scpca_sce.rds")
processed_url <- paste0(url_base, "SCPCL000001_processed.rds")
download.file(processed_url, test_processed_path, mode = "wb")

test_filtered_path <- testthat::test_path("data", "filtered_sce.rds")
filtered_url <- paste0(url_base, "SCPCL000001_filtered.rds")
download.file(filtered_url, test_filtered_path, mode = "wb")
