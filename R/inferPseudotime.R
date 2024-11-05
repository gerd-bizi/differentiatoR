# #' Infer Pseudotime Trajectory
# #'
# #' This function infers pseudotime ordering of samples.
# #'
# #' @param se A SummarizedExperiment object with counts and metadata.
# #' @param reduced_dims A data frame with reduced dimensions (from performClustering).
# #' @param cluster_assignments Optional cluster labels for samples.
# #' @return A vector of pseudotime values for each sample.
# #' @importFrom SummarizedExperiment assay
# #' @importFrom SingleCellExperiment SingleCellExperiment reducedDims<-
# #' @importFrom slingshot slingshot slingPseudotime
# #' @importFrom S4Vectors SimpleList
# #' @examples
# #' raw_counts <- matrix(rpois(100, 10), nrow = 15)
# #' normalized_counts <- normalizeCounts(raw_counts, method = "TMM")
# #' metadata <- data.frame(
# #'  sample_id = colnames(normalized_counts),
# #' timepoint = rep(c("Day1", "Day2", "Day3"), each = 5)
# #' )
# #' se <- assignSampleMetadata(normalized_counts, metadata)
# #' reduced_dims <- performClustering(normalized_counts, method = "PCA")
# #' pseudotime <- inferPseudotime(se, reduced_dims)
# #' @export
# inferPseudotime <- function(se, reduced_dims, cluster_assignments = NULL) {
#     counts <- SummarizedExperiment::assay(se)
#     if (is.null(cluster_assignments)) {
#         num_samples <- nrow(reduced_dims)
#         num_clusters <- min(3, num_samples - 1)
#         cluster_assignments <- stats::kmeans(reduced_dims, centers = num_clusters)$cluster
#     }
#     sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
#     SingleCellExperiment::reducedDims(sce) <- S4Vectors::SimpleList(PCA = as.matrix(reduced_dims))
#     sce <- slingshot::slingshot(sce, clusterLabels = cluster_assignments, reducedDim = 'PCA')
#     pseudotime <- slingshot::slingPseudotime(sce)
#     return(pseudotime)
# }