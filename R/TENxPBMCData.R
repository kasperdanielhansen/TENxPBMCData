TENxPBMCData <- function(dataset = c("pbmc68k", "frozen_pbmc_donor_a",
                                     "frozen_pbmc_donor_b", "frozen_pbmc_donor_c",
                                     "pbmc33k", "pbmc3k", "pbmc6k", "pbmc4k",
                                     "pbmc8k"))
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
    h5array <- HDF5Array(h5file, "counts")
    
    SingleCellExperiment(
      list(counts = h5array), rowData = rowData, colData = colData
    )
}
