#' Normalize Counts Using CPM or TMM
#'
#' This function normalizes raw count data using CPM or TMM methods.
#'
#' @importFrom edgeR cpm DGEList calcNormFactors
#' @param counts A matrix or data frame of raw counts (genes x samples).
#' @param method The normalization method to use: "CPM" or "TMM".
#' @return A normalized counts matrix.
#' @examples
#' raw_counts <- matrix(rpois(100, 10), nrow = 10)
#' normalized_counts <- normalizeCounts(raw_counts, method = "TMM")
#' @export
normalizeCounts <- function(counts, method = "TMM") {
    if (!method %in% c("CPM", "TMM")) {
        stop("Method must be 'CPM' or 'TMM'.")
    }
    if (method == "CPM") {
        norm_counts <- edgeR::cpm(counts)
    } else if (method == "TMM") {
        dge <- edgeR::DGEList(counts = counts)
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        norm_counts <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
    }
    log_norm_counts <- log2(norm_counts + 0.1)
    return(log_norm_counts)
}