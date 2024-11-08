% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotVolcano.R
\name{plotVolcano}
\alias{plotVolcano}
\title{Plot Volcano Plot for DE Results}
\usage{
plotVolcano(de_result, p_cutoff = 0.05, logfc_cutoff = 1)
}
\arguments{
\item{de_result}{A data frame containing DE results for a comparison.}

\item{p_cutoff}{Significance cutoff for adjusted p-value (default 0.05).}

\item{logfc_cutoff}{Cutoff for log2 fold change (default 1).}
}
\value{
A ggplot object of the volcano plot.
}
\description{
This function generates a volcano plot for DE results.
}
\examples{
de_results <- performDEAnalysis(se)
plotVolcano(de_results[["Day1_vs_Day2"]])
plotVolcano(de_results[["Day1_vs_Day2"]])
}
