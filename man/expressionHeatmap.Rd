\name{expressionHeatmap}
\alias{expressionHeatmap}
\title{
Heatmap of Gene Expression Matrix
}
\description{
Turn a matrxi of gene expression into a heatmap plot
}
\usage{
expressionHeatmap(m, Ngenes = NULL, minIntensity = 0, delta = 1, 
		heatcolors = redgreen(75), ...)
}
\arguments{
  \item{m}{
Matrix of gene expression data, as from \code{\link{expressionFileSetToMatrix}}
}
  \item{Ngenes}{
Limit the plot to a subset of genes, keeping just the \code{Ngenes} with the greatest variance
}
  \item{minIntensity}{
Minimum gene expression value, passed to \code{\link{expressionMatrixToMvalue}}
}
  \item{delta}{
Gene expression offset value, passed to \code{\link{expressionMatrixToMvalue}}
}
  \item{heatcolors}{
Vector of colors for the heatmap
}
  \item{\dots}{
Other arguments, passed to the heatmap plotting function
}
}
\value{
Invisibly, the return value from \code{\link{heatmap}}
}
