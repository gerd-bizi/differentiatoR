#' Assign Sample Metadata
#'
#' This function assigns metadata (e.g., timepoint, trial) to samples.
#'
#' @param counts A counts matrix with samples as columns.
#' @param metadata A data frame containing sample metadata.
#' @return A SummarizedExperiment object with counts and metadata.
#' @examples
#' se <- assignSampleMetadata(counts, metadata)
#' @export
assignSampleMetadata <- function(counts, metadata) {
  library(SummarizedExperiment)
  if (!all(colnames(counts) %in% metadata$sample_id)) {
    stop("Sample IDs in counts and metadata do not match.")
  }
  metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    colData = metadata
  )
  return(se)
}