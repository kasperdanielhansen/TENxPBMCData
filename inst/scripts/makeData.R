library(DropletUtils)
library(org.Hs.eg.db)
library(BiocFileCache)
library(HDF5Array)
library(beachmat)
library(Matrix)

bfc <- BiocFileCache(ask=FALSE)
set.seed(1000)

#########################

create_sce <- function(url, version = c(1, 2),
                       genome_build, sample_name, coldata, large=FALSE) {
    path <- bfcrpath(bfc, url)
    tmp <- tempfile()
    untar(path, exdir=tmp)
    if(version == 1) {
        subdir <- "filtered_matrices_mex"
    } else if(version == 2) {
        subdir <- "filtered_gene_bc_matrices"
    } else if(version == 3) {
      subdir <- "filtered_feature_bc_matrix"
      genome_build <- ""
    }
    sce <- read10xCounts(file.path(tmp, subdir, genome_build))
    symb <- mapIds(org.Hs.eg.db, keys=rownames(sce),
                   keytype="ENSEMBL", column="SYMBOL")
    colnames(rowData(sce)) <- c("ENSEMBL_ID", "Symbol_TENx")
    rowData(sce)$Symbol <- symb

    colData(sce)$Sequence <- sub("-.*", "", colData(sce)$Barcode)
    colData(sce)$Library <- as.integer(sub(".*-", "", colData(sce)$Barcode))
    colData(sce)$Sample <- sample_name
    colData(sce)$Cell_ranger_version <- coldata["version"]
    colData(sce)$Tissue_status <- coldata["tissue_type"]
    colData(sce)$Barcode_type <- coldata["barcode_type"]
    colData(sce)$Chemistry <- coldata["chemistry"]
    colData(sce)$Sequence_platform <- coldata["sequencer"]
    colData(sce)$Individual <- coldata["donor"]
    colData(sce)$Date_published <- coldata["published_date"]

    if(large) {
      mat <- matrix(0L, nrow = nrow(sce), ncol = ncol(sce))
      mat[which(counts(sce)>0)] <- as.integer(counts(sce)[which(counts(sce)>0)])
      counts(sce) <- mat
    } else {
      counts(sce) <- as.matrix(counts(sce))
    }

    mode(counts(sce)) <- "integer"

    return(sce)
}

save_sce <- function(sce, sample_name) {
  saveRDS(
    as.data.frame(rowData(sce)),
    paste0(sample_name, "_rowData.rds")
  )

  saveRDS(
    as.data.frame(colData(sce)),
    paste0(sample_name, "_colData.rds")
  )

  # Converting into a `HDF5Matrix` object

  options(DelayedArray.block.size=1e9) # 1GB block size.
  mat.out <- writeHDF5Array(
    counts(sce),
    filepath=paste0(sample_name, "_rectangular.h5"),
    name="counts",
    chunkdim=beachmat::getBestChunkDims(dim(sce))
  )
}

# pbmc4k
url <- "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz"
genome_build <- "GRCh38"
version <- 2
sample_name <- "pbmc4k"
coldata <- c(tissue_type = NA,
             version = "v2.1",
             chemistry = "Chromium_v2",
             barcode_type = "Chromium",
             sequencer = "Hiseq4000",
             donor = "HealthyDonor1",
             published_date = "2017-11-08")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "pbmc4k")

# pbmc8k
url <- "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz"
genome_build <- "GRCh38"
version <- 2
sample_name <- "pbmc8k"
coldata <- c(tissue_type = NA,
             version = "v2.1.0",
             chemistry = "Chromium_v2",
             barcode_type = "Chromium",
             sequencer = "Hiseq4000",
             donor = "HealthyDonor1",
             published_date = "2017-11-08")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "pbmc8k")

# pbmc33k
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc33k/pbmc33k_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 2
sample_name <- "pbmc33k"
coldata <- c(tissue_type = NA,
             chemistry = "Chromium_v1",
             version = "v2.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "HealthyDonor2",
             published_date = "2016-09-29")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "pbmc33k")

# pbmc3k
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 2
sample_name <- "pbmc3k"
coldata <- c(tissue_type = NA,
             chemistry = "Chromium_v1",
             version = "v1.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "HealthyDonor2",
             published_date = "2016-05-26")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "pbmc3k")

# pbmc6k
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc6k/pbmc6k_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 1
sample_name <- "pbmc6k"
coldata <- c(tissue_type = NA,
             chemistry = "Chromium_v1",
             version = "v1.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "HealthyDonor2",
             published_date = "2016-07-31")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "pbmc6k")

# pbmc68k
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 1
sample_name <- "pbmc68k"
coldata <- c(tissue_type = "Fresh",
             chemistry = "GemCode",
             version = "v1.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "DonorA",
             published_date = "2016-07-24")
sce <- create_sce(url, version, genome_build, sample_name, coldata, large=TRUE)
save_sce(sce, "pbmc68k")

# frozen_pbmc_donor_a
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/frozen_pbmc_donor_a/frozen_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 1
sample_name <- "frozen_pbmc_donor_a"
coldata <- c(tissue_type = "Frozen",
             chemistry = "GemCode",
             version = "v1.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "DonorA",
             published_date = "2016-07-24")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "frozen_pbmc_donor_a")

# frozen_pbmc_donor_b
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/frozen_pbmc_donor_b/frozen_pbmc_donor_b_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 1
sample_name <- "frozen_pbmc_donor_b"
coldata <- c(tissue_type = "Frozen",
             chemistry = "GemCode",
             version = "v1.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "DonorB",
             published_date = "2016-07-24")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "frozen_pbmc_donor_b")

# frozen_pbmc_donor_c
url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/frozen_pbmc_donor_c/frozen_pbmc_donor_c_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"
version <- 1
sample_name <- "frozen_pbmc_donor_c"
coldata <- c(tissue_type = "Frozen",
             chemistry = "GemCode",
             version = "v1.1.0",
             barcode_type = "GemCode",
             sequencer = "NextSeq500",
             donor = "DonorC",
             published_date = "2016-07-24")
sce <- create_sce(url, version, genome_build, sample_name, coldata)
save_sce(sce, "frozen_pbmc_donor_c")
