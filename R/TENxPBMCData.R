TENxPBMCData <- function(sample_name = c("pbmc4k"))
  ## Download HDF5 (dense assay) and RDS (row and column
  ## annotations) files from ExperimentHub, compose into a
  ## SingleCellExperiment object.
{
    hub <- ExperimentHub()
    base <- sample_name
    
    ## row and column data
    rdatapath <- paste0("TENxPBMCData/", base, "_rowData.rds")
    rowData <- query(hub, rdatapath)[[1]]
    
    suppressMessages({
      rdatapath <- paste0("TENxPBMCData/", base, "_colData.rds")
      colData <- query(hub, rdatapath)[[1]]
      
      ## HDF5, from ExperimentHub:
      rdatapath <- paste0("TENxPBMCData/", base, "_rectangular.h5")
      h5file <- query(hub, rdatapath)[[1]]
    })
    h5array <- HDF5Array(h5file, "counts")
    
    SingleCellExperiment(
      list(counts = h5array), rowData = rowData, colData = colData
    )
}
