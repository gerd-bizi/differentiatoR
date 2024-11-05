#' Plot Volcano Plot for DE Results
#'
#' This function generates a volcano plot for DE results.
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_manual theme_minimal labs
#' @param de_result A data frame containing DE results for a comparison.
#' @param p_cutoff Significance cutoff for adjusted p-value (default 0.05).
#' @param logfc_cutoff Cutoff for log2 fold change (default 1).
#' @return A ggplot object of the volcano plot.
#' @examples
#' de_results <- performDEAnalysis(se)
#' plotVolcano(de_results[["Day1_vs_Day2"]])
#' plotVolcano(de_results[["Day1_vs_Day2"]])
#' @export
plotVolcano <- function(de_result, p_cutoff = 0.05, logfc_cutoff = 1) {
  de_result$Significant <- with(de_result, ifelse(
    padj < p_cutoff & abs(log2FoldChange) > logfc_cutoff, "Yes", "No"
  ))
  plot <- ggplot2::ggplot(de_result, aes(x = log2FoldChange, y = -log10(padj))) +
    ggplot2::geom_point(aes(color = Significant)) +
    ggplot2::scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
  return(plot)
}