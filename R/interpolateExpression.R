#' Interpolate Gene Expression at a Selected Timepoint
#'
#' This function interpolates gene expression levels at a specified timepoint.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @param se A SummarizedExperiment object with counts and metadata.
#' @param timepoint The timepoint at which to interpolate expression.
#' @return A numeric value of the interpolated expression level.
#' @examples
#' counts <- matrix(rpois(100, 15), nrow = 15)
#' colnames(counts) <- paste0("sample", 1:15)
#' metadata <- data.frame(
#' sample_id = colnames(counts),
#' timepoint = rep(c(0, 5, 10), each = 5)
#' )
#' se <- assignSampleMetadata(counts, metadata)
#' interp_expr <- interpolateExpression(se, timepoint = 2.5)
#' @export
interpolateExpression <- function(se, timepoint) {
    counts <- SummarizedExperiment::assay(se)
    metadata <- as.data.frame(colData(se))

    # Figure out the closest timepoints
    timepoints <- unique(metadata$timepoint)
    closest_timepoints <- c(min(timepoints), max(timepoints))
    if (timepoint < min(timepoints) | timepoint > max(timepoints)) {
        stop("Timepoint is outside the range of available timepoints.")
    }
    if (timepoint %in% timepoints) {
        stop("Timepoint already exists in the data.")
    }

    # First, fit the expression levels to a third-degree polynomial, accounting for
    # the timepoint as a factor
    fit <- lm(counts ~ poly(metadata$timepoint, 3), data = metadata)

    # Predict the expression level at the specified timepoint
    pred <- predict(fit, newdata = data.frame(timepoint = timepoint))

    return(pred)
}