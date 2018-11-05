\name{plotTranscripts}
\alias{plotTranscripts}
\alias{plotRatios}
\title{
Scatterplot or MA plot
}
\description{
Scatter or MA plots of trancriptomes or differential expression (DE) files
}
\usage{
plotTranscripts(file1, file2, fid1 = basename(file1), fid1 = basename(file1), 
		geneColumn = "GENE_ID", value1Column = "RPKM_M", 
		value2Column = value1Column, cex = 1, units1 = "RPKM", units2 = units1,
		offset = 1, keepIntergenics = FALSE, 
		label = "Plot", plotType = c("Scatter", "MA"), 
		marker.genes = NULL, marker.col = 1, marker.cex = 1, 
		marker.labels = TRUE, marker.pch = 1,
		minYmax = NULL, sep = "\t", sym.asp = TRUE, hideZero = FALSE,
		showAllGenes = FALSE)

plotRatios(file, geneColumn = "GENE_ID", value1Column = "RPKM_1_M", value2Column = "RPKM_2_M", 
		cex = 1, units = "RPKM", offset = 1, keepIntergenics = FALSE, 
		label = "Plot", plotType = c("Scatter", "MA"), 
		marker.genes = NULL, marker.col = 1, marker.cex = 1, 
		marker.labels = TRUE, marker.pch = 1,
		minYmax = NULL, sep = "\t", sym.asp = TRUE, lab1 = "", lab2 = "",
		hideZero = FALSE)

plotRatiosTwoFiles(file1, file2, fid1 = basename(file1), fid1 = basename(file1), 
		geneColumn = "GENE_ID", valueColumn = "LOG_2_FOLD_M", 
		cex = 1, units = "Log2 Fold", offset = 0, keepIntergenics = FALSE, 
		label = "Plot", marker.genes = NULL, marker.col = 1, marker.cex = 1, 
		marker.labels = TRUE, marker.pch = 1,
		minYmax = NULL, sep = "\t", sym.asp = TRUE, hideZero = FALSE)

}
\arguments{
  \item{file, file1, file2}{
character string, full pathname to existing transcript or DE file
}
  \item{fid1, lab1}{
sampleID text label for file 1, added to X axis
}
  \item{fid2, lab2}{
sampleID text label for file 2, added to Y axis
}
  \item{geneColumn}{
name of the column that contains GeneIDs
}
  \item{valueColumn}{
name of the column(s) that contains expression (or DE) data
}
  \item{cex}{
expansion factor for plotted points
}
  \item{units}{
character string of expression units, added to plot labels
}
  \item{offset}{
linear offset added to expression data, to prevent divide by zero during log transform, etc.
}
  \item{keepIntergenics}{
logical, keep the explicit 'non-genes' or drop them from the plot
}
  \item{label}{
text for the main plot label
}
  \item{plotType}{
type of plot. Either scatterplot or the MA plot typical of microarrays
}
  \item{marker.genes}{
character vector of GeneIDs to be highlighted, or \code{"identify"} to select by mouse
clicks.
}
  \item{marker.col}{
color for highlighted gene points
}
  \item{marker.cex}{
expansion factor for highlighted genes
}
  \item{marker.labels}{
logical, add the GeneIDs for the marker genes to the plot
}
  \item{marker.pch}{
'pch' parameter for the symbol for the marker genes
}
  \item{minYmax}{
explicit control for the Y axis upper, useful for making making identical
scaling over multiple plots
}
  \item{hideZero}{
logical, should genes with no expression be omitted.
}
}
\details{
Function to standardize plotting of transcriptomes and DE result files.  Suitable for
comparing RNA-seq and microarray results.  The 'Transcript' form can combine any 2
files of expression data.  The 'Ratios' form takes its 2 transcripts from different 
columns within one differential expression file. The 'RatiosTwoFiles' form take
differential expression fold change values from two separate files.
}
\value{
In addition to a plot, these function return a list of data (subject to change) that
currently includes:
  \item{x,y }{ the coordinates of the plotted points}
  \item{id }{ the GeneIDs of those points}
  \item{Pearson_R }{ the correlation coefficient}
  \item{Spearman_Rho }{ the correlation coefficient}
}
