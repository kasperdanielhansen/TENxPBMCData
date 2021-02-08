TENxPBMCData <- function(dataset = c("pbmc4k", "pbmc68k", "frozen_pbmc_donor_a",
                                     "frozen_pbmc_donor_b", "frozen_pbmc_donor_c",
                                     "pbmc33k", "pbmc3k", "pbmc6k",
                                     "pbmc8k", "pbmc5k-CITEseq"),
                         as.sparse = TRUE)
{
    ## Download HDF5 (dense assay) and RDS (row and column
    ## annotations) files from ExperimentHub, compose into a
    ## SingleCellExperiment object.
    dataset <- match.arg(dataset)
    hub <- ExperimentHub()
    version <- "v1.0.0"
    base <- paste0("TENxPBMCData/", version, "/")
    
    ## row and column data
    rdatapath <- paste0(base, dataset, "_rowData.rds")
    rowData <- query(hub, rdatapath)[[1]]
    
    suppressMessages({
      rdatapath <- paste0(base, dataset, "_colData.rds")
      colData <- query(hub, rdatapath)[[1]]
      
      ## HDF5, from ExperimentHub:
      rdatapath <- paste0(base, dataset, "_rectangular.h5")
      h5file <- query(hub, rdatapath)[[1]]
    })
    h5array <- HDF5Array(h5file, "counts", as.sparse = as.sparse)
    
    sce <- SingleCellExperiment(
      list(counts = h5array), rowData = rowData, colData = colData
    )
    
    # if multiple Types present in rowData, make altExp for SCE
    if (length(unique(rowData(sce)$Type)) > 1) {
      sce <- splitAltExps(sce, rowData(sce)$Type)
      
      # modify rowData of altExp for CITE-seq data
      rowData(altExp(sce))$ENSEMBL_ID <- NULL
      rowData(altExp(sce))$Symbol <- rownames(rowData(altExp(sce)))
    }
    
    return(sce)
}
