#' Cluster Samples Within a Timepoint
#'
#' This function clusters samples at a given timepoint and identifies top genes.
#'
#' @param se A SummarizedExperiment object containing counts and metadata.
#' @param timepoint The timepoint to analyze.
#' @param n_genes Number of top genes to identify (default 10).
#' @return A list containing cluster assignments and top genes.
#' @importFrom SummarizedExperiment assay colData
#' @importFrom stats prcomp kmeans
#' @examples
#' raw_counts <- matrix(rpois(100, 10), nrow = 10)
#' normalized_counts <- normalizeCounts(raw_counts, method = "TMM")
#' metadata <- data.frame(
#'   sample_id = colnames(normalized_counts),
#'   timepoint = rep(c("Day1", "Day2"), each = 5)
#' )
#' se <- assignSampleMetadata(normalized_counts, metadata)
#' cluster_res <- clusterWithinTimepoint(se, timepoint = "Day1")
#' @export
clusterWithinTimepoint <- function(se, timepoint, n_genes = 10) {
  counts <- SummarizedExperiment::assay(se)
  metadata <- SummarizedExperiment::colData(se)
  idx <- metadata$timepoint == timepoint
  counts_sub <- counts[, idx]
  
  # Check if counts_sub is a matrix
  if (!is.matrix(counts_sub)) {
    counts_sub <- as.matrix(counts_sub)
  }
  
  # Perform PCA
  pca_res <- stats::prcomp(t(counts_sub), scale. = TRUE)
  
  # Cluster samples
  cluster_res <- stats::kmeans(pca_res$x[, 1:5], centers = 2)
  
  # Identify top genes contributing to PC1
  loading_scores <- abs(pca_res$rotation[, 1])
  top_genes <- names(sort(loading_scores, decreasing = TRUE))[1:n_genes]
  
  result <- list(
    clusters = cluster_res$cluster,
    top_genes = top_genes
  )
  return(result)
}