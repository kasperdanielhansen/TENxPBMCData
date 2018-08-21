## This is skeleton

library(DropletUtils)
library(org.Hs.eg.db)
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
set.seed(1000)

#########################

create_sce <- function(url, genome_build) {
    path <- bfcrpath(bfc, url)
    tmp <- tempfile()
    untar(path, exdir=tmp)
    sce <- read10xCounts(file.path(tmp, "filtered_matrices_mex",
                                   genome_build))
    symb <- mapIds(org.Hs.eg.db, keys=rownames(sce),
                   keytype="ENSEMBL", column="SYMBOL")
    colnames(rowData(sce)) <- c("ENSEMBL", "Symbol_TENx")
    rowData(sce)$Symbol <- symb
    return(sce)
}

url <- "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"
genome_build <- "hg19"

create_sce(url, genome_build)
