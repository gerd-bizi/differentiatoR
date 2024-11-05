#' Calculate Temporal Smoothness Metric and Return Individual Scores
#'
#' This function computes the Temporal Smoothness Metric (TSM) for each pair of consecutive time points.
#'
#' @param cluster_assignments A vector of cluster labels.
#' @param time_points A vector of time points corresponding to each sample.
#' @return A data frame with the temporal smoothness scores for each pair of time points.
#' @examples
#' cluster_assignments <- c(1, 1, 2, 2, 3, 3, 1, 1, 2, 2)
#' time_points <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
#' tsm_df <- calculateTemporalSmoothness(cluster_assignments, time_points)
#' @export
calculateTemporalSmoothness <- function(cluster_assignments, time_points) {
  unique_times <- sort(unique(time_points))
  n_times <- length(unique_times)
  ts_scores <- numeric(n_times - 1)
  comparisons <- character(n_times - 1)
  
  for (i in 1:(n_times - 1)) {
    t1 <- unique_times[i]
    t2 <- unique_times[i + 1]
    comparisons[i] <- paste(t1, "to", t2, sep = " ")
    
    idx1 <- which(time_points == t1)
    idx2 <- which(time_points == t2)
    
    clusters1 <- cluster_assignments[idx1]
    clusters2 <- cluster_assignments[idx2]
    
    # Create contingency table
    ct <- table(clusters1, clusters2)
    n_same <- sum(diag(ct))
    total <- sum(ct)
    ts_scores[i] <- n_same / total
  }
  
  # Create a data frame with the results
  tsm_df <- data.frame(
    Comparison = comparisons,
    TemporalSmoothness = ts_scores
  )
  
  return(tsm_df)
}