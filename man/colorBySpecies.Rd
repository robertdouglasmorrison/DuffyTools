\name{colorBySpecies}
\alias{colorBySpecies}
\title{
Assign colors by SpeciesID
}
\description{
Assigns unique colors by SpeciesID, for multi-species plotting.
}
\usage{
colorBySpecies(speciesSet, palette = c("cyan", "springgreen", "gold", "hotpink"), 
	intergenic = NULL, intergenic.color = "brown")
}
\arguments{
  \item{speciesSet}{
a vector of SpeciesIDs, or similar identifiers to factor on
}
  \item{palette}{
a vector of colors to choose from, at least as long as the levels of \code{speciesSet}
}
  \item{intergenic}{
an optional boolean condition, of the same length as \code{speciesSet}, that flags which
IDs will be colored as 'intergenic'
}
  \item{intergenic.color}{
the color for 'intergenic'
}
}
\details{
a utility for coloring features of multi-species images, like chromosome or gene pie and bar charts.
Lets the speciesID of each feature be color coded without explicitly naming them.
}
\value{
a LIST, with
  \item{colors }{ a vector of colors of the same length as \code{speciesSet}}
  \item{legend.colors }{ a named vector of the unique colors and their SpeciesIDs, for use in figure legends}
}
\author{
Bob Morrison
}
