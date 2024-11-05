#' Perform Differential Expression Analysis with edgeR
#'
#' This function performs DE analysis between all pairwise timepoints using the edgeR package.
#'
#' @importFrom edgeR DGEList filterByExpr calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom stats model.matrix
#' @param se A SummarizedExperiment object with counts and metadata.
#' @return A list of DE results data frames for each pairwise comparison.
#' @examples
#' counts <- matrix(rpois(100, 10), nrow = 10)
#' colnames(counts) <- paste0("sample", 1:10)
#' metadata <- data.frame(
#' sample_id = colnames(counts),
#' timepoint = rep(c("Day1", "Day2", "Day3"), each = 5)
#' )
#' se <- assignSampleMetadata(counts, metadata)
#' de_results <- performDEAnalysis(se)
#' @export
performDEAnalysis <- function(se) {
  # Extract counts and metadata
  counts <- SummarizedExperiment::assay(se)
  metadata <- as.data.frame(colData(se))
  
  # Ensure 'timepoint' is a factor
  if (!"timepoint" %in% colnames(metadata)) {
    stop("The colData of the SummarizedExperiment must contain a 'timepoint' column.")
  }
  metadata$timepoint <- factor(metadata$timepoint)
  
  # Create DGEList object
  dge <- edgeR::DGEList(counts = counts, samples = metadata)
  
  # Filter lowly expressed genes
  keep <- edgeR::filterByExpr(dge, group = metadata$timepoint)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Recalculate library sizes after filtering
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  
  # Design matrix for all timepoints
  design <- stats::model.matrix(~ 0 + timepoint, data = metadata)
  colnames(design) <- levels(metadata$timepoint)
  
  # Estimate dispersions
  dge <- edgeR::estimateDisp(dge, design)
  
  # Fit the generalized linear model (GLM)
  fit <- edgeR::glmQLFit(dge, design)
  
  # Perform pairwise comparisons
  timepoints <- levels(metadata$timepoint)
  results_list <- list()
  
  for (i in 1:(length(timepoints) - 1)) {
    for (j in (i + 1):length(timepoints)) {
      tp1 <- timepoints[i]
      tp2 <- timepoints[j]
      
      # Create contrast vector
      contrast <- rep(0, length(timepoints))
      names(contrast) <- timepoints
      contrast[tp1] <- 1
      contrast[tp2] <- -1
      
      # Perform quasi-likelihood F-test
      qlf <- edgeR::glmQLFTest(fit, contrast = contrast)
      
      # Get top tags (DE results)
      de_result <- edgeR::topTags(qlf, n = Inf)
      de_table <- as.data.frame(de_result)
      de_table$gene <- rownames(de_table)
      
      # Add to results list
      comparison_name <- paste(tp1, "_vs_", tp2, sep = "")
      results_list[[comparison_name]] <- de_table
    }
  }
  
  return(results_list)
}