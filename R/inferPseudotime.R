#' Infer Pseudotime Trajectory
#'
#' This function infers pseudotime ordering of samples.
#'
#' @param se A SummarizedExperiment object with counts and metadata.
#' @param reduced_dims A data frame with reduced dimensions (from performClustering).
#' @param cluster_assignments Optional cluster labels for samples.
#' @return A vector of pseudotime values for each sample.
#' @examples
#' pseudotime_values <- inferPseudotime(se, dr_results)
#' @export
inferPseudotime <- function(se, reduced_dims, cluster_assignments = NULL) {
    library(slingshot)
    counts <- assay(se)
    if (is.null(cluster_assignments)) {
        cluster_assignments <- kmeans(reduced_dims, centers = 3)$cluster
    }
    sce <- SingleCellExperiment(assays = list(counts = counts))
    reducedDims(sce) <- SimpleList(PCA = as.matrix(reduced_dims))
    sce <- slingshot(sce, clusterLabels = cluster_assignments, reducedDim = 'PCA')
    pseudotime <- slingPseudotime(sce)
    return(pseudotime)
}