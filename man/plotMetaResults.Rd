\name{plotMetaResults}
\alias{plotMetaResults}
\title{
Deviation Plots for Meta Results
}
\description{
Generate a few visual assessments of Meta Results variance among the DE tools
}
\usage{
plotMetaResults(file, columnNames = c("RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
		columnColors = rainbow(length(columnNames), end = 0.72), label = "")
}
\arguments{
  \item{file}{
Full path name to an existing file containing a meta results or meta rank table
}
  \item{columnNames}{
Character vector of column names to use in meta results assessment
}
  \item{columnColors}{
colors for each column
}
  \item{label}{
Optional text to append to plot top legend
}
}
\details{
This function assesses the variance between each DE tools call for gene rankings against the consensus
result, and then presents the data visually in a few plot styles.
}
\value{
Plot files are written to the current directory.
}
\note{
In the line plot image, we expect a 'fish' shaped result, where all tools agree somewhat about the most
up-regulated and down-regulated genes but have higher variance for the non differentially expressed genes 
in the middle of the file.
}
