#' Normalize Counts Using CPM or TMM
#'
#' This function normalizes raw count data using CPM or TMM methods.
#'
#' @param counts A matrix or data frame of raw counts (genes x samples).
#' @param method The normalization method to use: "CPM" or "TMM".
#' @return A normalized counts matrix.
#' @examples
#' normalized_counts <- normalizeCounts(raw_counts, method = "TMM")
#' @export
normalizeCounts <- function(counts, method = "TMM") {
    if (!method %in% c("CPM", "TMM")) {
        stop("Method must be 'CPM' or 'TMM'.")
    }
    # Add 0.1 to all counts to avoid log(0) errors
    counts <- counts + 0.1
    if (method == "CPM") {
        library(edgeR)
        counts_cpm <- cpm(counts)
        return(counts_cpm)
    } else if (method == "TMM") {
        library(edgeR)
        dge <- DGEList(counts = counts)
        dge <- calcNormFactors(dge, method = "TMM")
        counts_tmm <- cpm(dge, normalized.lib.sizes = TRUE)
        return(counts_tmm)
    }
}