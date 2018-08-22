## This is skeleton

library(DropletUtils)
library(org.Hs.eg.db)
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
set.seed(1000)

#########################

create_sce <- function(url, version = c("v1.1", "v2.1"),
                       genome_build, coldata) {
    path <- bfcrpath(bfc, url)
    tmp <- tempfile()
    untar(path, exdir=tmp)
    if(version == "v1.1") {
        subdir <- "filtered_matrices_mex"
    } else if(version == "v2.1") {
        subdir <- "filtered_gene_bc_matrices"
    }
    sce <- read10xCounts(file.path(tmp, subdir, genome_build))
    symb <- mapIds(org.Hs.eg.db, keys=rownames(sce),
                   keytype="ENSEMBL", column="SYMBOL")
    colnames(rowData(sce)) <- c("ENSEMBL", "Symbol_TENx")
    rowData(sce)$Symbol <- symb

    colData(sce)$Sequence <- sub("-.*", "", colData(sce)$Barcode)
    colData(sce)$Library <- sub(".*-", "", colData(sce)$Barcode)
    colData(sce) <- cbind(colData(sce), coldata)
    return(sce)
}

url <- "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz"
genome_build <- "GRCh38"

create_sce(url, genome_build, coldata = NULL)
