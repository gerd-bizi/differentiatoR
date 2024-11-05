#' Assign Sample Metadata
#'
#' This function assigns metadata (e.g., timepoint, trial) to samples.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @param counts A counts matrix with samples as columns.
#' @param metadata A data frame containing sample metadata.
#' @return A SummarizedExperiment object with counts and metadata.
#' @examples
#' counts <- matrix(rpois(100, 10), nrow = 10)
#' colnames(counts) <- paste0("sample", 1:10)
#' metadata <- data.frame(
#'  sample_id = colnames(counts),
#' timepoint = rep(c("Day1", "Day2"), each = 5)
#' )
#' se <- assignSampleMetadata(counts, metadata)
#' @export
assignSampleMetadata <- function(counts, metadata) {
  if (!all(colnames(counts) %in% metadata$sample_id)) {
    stop("Sample IDs in counts and metadata do not match.")
  }
  metadata <- metadata[match(colnames(counts), metadata$sample_id), ]
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(counts)),
    colData = metadata
  )
  return(se)
}