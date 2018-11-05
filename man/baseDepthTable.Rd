\name{baseDepthTable}
\alias{baseDepthTable}
\alias{baseDepthVector}
\alias{baseDepthTableSubset}
\alias{baseDepthVectorToTable}
\alias{baseDepthTableToVector}
\alias{sub.baseDepthTable}
\alias{sub.baseDepthTableSubset}
\alias{merge.baseDepthTables}
\alias{factorBaseLocations}
\alias{emptyBaseDepthTable}
\alias{emptyBaseDepthVector}
\title{
Base Depth Tables and Vectors
}
\description{
Routines and objects to manage Base Depth data, such as Wiggle tracks
}
\usage{
baseDepthTable(start, stop, weight = 1)
baseDepthVector(start, stop, weight = 1)

baseDepthTableToVector(baseTbl)
baseDepthVectorToTable(baseVec)

baseDepthTableSubset(baseTbl, start, stop)

sum.baseDepthTable(baseTbl)
sum.baseDepthTableSubset(baseTbl, start, stop)

merge.baseDepthTables(baseTbl1, baseTbl2)

emptyBaseDepthTable()
emptyBaseDepthVector()

factorBaseLocations(x)
}
\arguments{
  \item{start}{
an integer vector of starting locations
}
  \item{stop}{
an integer vector of stopping locations
}
  \item{weight}{
a numeric vector of weights (depths)
}
  \item{baseTbl}{
a data frame of base depths with elements START, STOP, DEPTH
}
  \item{baseVec}{
a numeric vector of base depths, with integer names representing 
the base locations, sorted by name (location)
}
  \item{x}{
an integer vector of base locations, potentially with duplicates and unordered
}
}
\details{
The Base Depth objects and functions implement the idea of wiggle tracks, where
short nucleotide sequences pile up on chromosomal locations.  The table form is
compact, while the vector form allows for easy merging, adding, etc.  Functions for
subsetting and counting bases are provided.
}
\value{
The default constructors \code{emptyBaseDepthTable} \code{emptyBaseDepthVector} 
return empty objects of the appropriate type.

The standard constructors turn a set of start/stop pairs with explicit depths
(weights) into either type of base depth object.

The converters turn a table into a vector or vice versa.

Subsetting returns a (usually) smaller table containing just the locations within
the given start,stop boundaries.

The 'sum' functions return the total number of bases in the table or subset.

the 'merge' function adds the depths from the 2 tables into a new merged table.
}
\note{
These routines track the depth (count) of bases at locations.  They do not
store any base identity information whatsoever.
}
\seealso{
  \code{\link{mergeIntegerTables}}, for merging (adding) base depth vectors.
}
