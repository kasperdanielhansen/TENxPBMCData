## This is skeleton

library(scran)
library(scater)
library(DropletUtils)

library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)    
set.seed(1000)

#########################

path.68 <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz")
tmp.68 <- tempfile()
untar(path.68, exdir=tmp.68)
sce.68 <- read10xCounts(file.path(tmp.68, "filtered_matrices_mex/hg19/")) 

library(org.Hs.eg.db)
symb <- mapIds(org.Hs.eg.db, keys=rownames(sce.68), keytype="ENSEMBL", column="SYMBOL")
rowData(sce.68)$Symbol <- symb


# Pre-processing the 4K T-cell dataset.
bfc <- BiocFileCache(ask=FALSE)    
path.4 <- bfcrpath(bfc, "http://cf.10xgenomics.com/samples/cell-exp/2.1.0/t_4k/t_4k_filtered_gene_bc_matrices.tar.gz")
tmp.4 <- tempfile()
untar(path.4, exdir=tmp.4)
sce.4 <- read10xCounts(file.path(tmp.4, "filtered_gene_bc_matrices/GRCh38/"))

symb <- mapIds(org.Hs.eg.db, keys=rownames(sce.4), keytype="ENSEMBL", column="SYMBOL")
rowData(sce.4)$Symbol <- symb

