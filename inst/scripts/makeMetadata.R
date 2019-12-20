## validate with `ExperimentHubData::makeExperimentHubMetadata()`
## (above pkg directory)

main.data <- data.frame(
  Title = c(
    "PBMC, 68k scRNA-seq data",
    "PBMC, 68k scRNA-seq sample (column) annotation",
    "PBMC, 68k scRNA-seq gene (row) annotation",
    "PBMC, Donor A scRNA-seq data",
    "PBMC, Donor A scRNA-seq sample (column) annotation",
    "PBMC, Donor A scRNA-seq gene (row) annotation",
    "PBMC, Donor B scRNA-seq data",
    "PBMC, Donor B scRNA-seq sample (column) annotation",
    "PBMC, Donor B scRNA-seq gene (row) annotation",
    "PBMC, Donor C scRNA-seq data",
    "PBMC, Donor C scRNA-seq sample (column) annotation",
    "PBMC, Donor C scRNA-seq gene (row) annotation",
    "PBMC, 33k scRNA-seq data",
    "PBMC, 33k scRNA-seq sample (column) annotation",
    "PBMC, 33k scRNA-seq gene (row) annotation",
    "PBMC, 3k scRNA-seq data",
    "PBMC, 3k scRNA-seq sample (column) annotation",
    "PBMC, 3k scRNA-seq gene (row) annotation",
    "PBMC, 6k scRNA-seq data",
    "PBMC, 6k scRNA-seq sample (column) annotation",
    "PBMC, 6k scRNA-seq gene (row) annotation",
    "PBMC, 4k scRNA-seq data",
    "PBMC, 4k scRNA-seq sample (column) annotation",
    "PBMC, 4k scRNA-seq gene (row) annotation",
    "PBMC, 8k scRNA-seq data",
    "PBMC, 8k scRNA-seq sample (column) annotation",
    "PBMC, 8k scRNA-seq gene (row) annotation",
    "PBMC, 5k CITE-Seq data",
    "PBMC, 5k CITE-Seq sample (column) annotation",
    "PBMC, 5k CITE-Seq gene (row) annotation"
  ),
  Description = c(as.vector(t(outer(
    paste("Single-cell RNA-seq data for human PBMC data.",
          c("68k", "Donor A", "Donor B", "Donor C",
            "33k", "3k", "6k", "4k", "8k"), "dataset."),
    c(
      "Full rectangular, block-compressed format, 1GB block size.",
      "Inferred sample descriptions.",
      "Gene descriptions."
    ), FUN = paste
  ))),
  as.vector(t(outer(
    paste("CITE-Seq data for human PBMC data.",
          c("5k"), "dataset."),
    c(
      "Full rectangular, block-compressed format, 1GB block size.",
      "Inferred sample descriptions.",
      "Gene descriptions."
    ), FUN = paste
  )))
  ),
  RDataPath = as.vector(t(outer(
    paste0("TENxPBMCData/v1.0.0/",
    c("pbmc68k", "frozen_pbmc_donor_a", "frozen_pbmc_donor_b",
      "frozen_pbmc_donor_c", "pbmc33k", "pbmc3k", "pbmc6k", "pbmc4k", "pbmc8k", "pbmc5k-CITEseq")),
    c("_rectangular.h5", "_colData.rds", "_rowData.rds"), FUN = paste0))),
  BiocVersion="3.8",
  Genome=c(rep("hg19", 21), rep("GRCh38", 9)),
  SourceType="tar (of mtx and tsv)",
  SourceUrl=rep(paste0("http://cf.10xgenomics.com/samples/cell-exp/",
                   c("1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz",
                   "1.1.0/frozen_pbmc_donor_a/frozen_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz",
                   "1.1.0/frozen_pbmc_donor_b/frozen_pbmc_donor_b_filtered_gene_bc_matrices.tar.gz",
                   "1.1.0/frozen_pbmc_donor_c/frozen_pbmc_donor_c_filtered_gene_bc_matrices.tar.gz",
                   "1.1.0/pbmc33k/pbmc33k_filtered_gene_bc_matrices.tar.gz",
                   "1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
                   "1.1.0/pbmc6k/pbmc6k_filtered_gene_bc_matrices.tar.gz",
                   "2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz",
                   "2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz",
                   "3.0.2/5k_pbmc_protein_v3/5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz"
                   
  )), each = 3),
  SourceVersion=rep(c("0c8c546cbadc8ea8eb45be6cfe040852",
                      "9d0a60d7d51d45e16e7c04cdb509879f",
                      "96bd95909d0409cce922800890805b6b",
                      "14654894a9f20704cef5db6301a0d793",
                      "fef79ae6064f704888df538aba3f1d5c",
                      "f76d73bf7a582aaf32b798bcae9574d3",
                      "4c9a44cd9386234b3168016945b5b072",
                      "f61f4deca423ef0fa82d63fdfa0497f7",
                      "cbb3a7553ece7eaddeb8a56df781ccb0",
                      "bad4139cc6ef9154806d2e879a2a2435"), each = 3),
  Species="Homo sapiens",
  TaxonomyId="9606",
  Coordinate_1_based=TRUE,
  DataProvider="10X Genomics",
  Maintainer="Kasper D. Hansen <kasperdanielhansen@gmail.com>",
  RDataClass="character",
  DispatchClass=rep(c("H5File", "Rds", "Rds"), 10),
  Chemistry=c(rep("Gemcode", 12), rep("Chromium_v1", 9), rep("Chromium_v2", 6), rep("Chromium_v3", 3)),
  stringsAsFactors = FALSE
)

write.csv(file="metadata.csv", main.data, row.names=FALSE)

