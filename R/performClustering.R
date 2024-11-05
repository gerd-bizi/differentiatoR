#' Perform Dimensionality Reduction and Clustering
#'
#' This function performs PCA, t-SNE, or UMAP on normalized counts.
#'
#' @importFrom stats prcomp
#' @importFrom Rtsne Rtsne
#' @importFrom umap umap
#' @param counts A normalized counts matrix (genes x samples).
#' @param method The dimensionality reduction method: "PCA", "tSNE", or "UMAP".
#' @param n_components Number of components to reduce to (default is 2).
#' @return A data frame with reduced dimensions for each sample.
#' @examples
#' raw_counts <- matrix(rpois(100, 10), nrow = 10)
#' normalized_counts <- normalizeCounts(raw_counts, method = "TMM")
#' dr_results <- performClustering(normalized_counts, method = "UMAP")
#' @export
performClustering <- function(counts, method = "PCA", n_components = 2) {
    if (method == "PCA") {
        pca_res <- stats::prcomp(t(counts))
        dr_df <- as.data.frame(pca_res$x[, 1:n_components])
        colnames(dr_df) <- paste0("PC", 1:n_components)
    } else if (method == "tSNE") {
        tsne_res <- Rtsne::Rtsne(t(counts), dims = n_components)
        dr_df <- as.data.frame(tsne_res$Y)
        colnames(dr_df) <- paste0("tSNE", 1:n_components)
    } else if (method == "UMAP") {
        umap_res <- umap::umap(t(counts), n_components = n_components)
        dr_df <- as.data.frame(umap_res$layout)
        colnames(dr_df) <- paste0("UMAP", 1:n_components)
    } else {
        stop("Method must be 'PCA', 'tSNE', or 'UMAP'.")
    }
    return(dr_df)
}