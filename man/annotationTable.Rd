\name{AnnotationTable}
\alias{AnnotationTable}
\alias{readAnnotationTable}
\alias{getAnnotationValue}
\alias{getAnnotationTrue}
\alias{whichAnnotationTrue}
\title{
Extract facts from an Annotation Table
}
\description{
Functions to query an Annotation Table for keyed values about one or more named sample datasets.
}
\usage{
readAnnotationTable(file, neededColumns = NULL, sep = "\t", verbose = FALSE)
getAnnotationValue(tbl, key, columnArg, notfound = "", verbose = TRUE)
getAnnotationTrue(tbl, key, columnArg, notfound = FALSE, verbose = TRUE)
whichAnnotationTrue(tbl, columnArg, negate = FALSE, verbose = TRUE)
}
\arguments{
  \item{file}{
Name of a sample annotation file.
}
  \item{neededColumns}{
Optional character vector of required column names that the annotation table must contain.
}
  \item{tbl}{
Either an annotation table previously read via 'readAnnotationTable' or a character string 
file name to open and read from.
}
  \item{key}{
The keyword value that determines which one row the the annotation to query.  This is required to be an
entry in the first (leftmost) column of the annotation table.
}
  \item{columnArg}{
The name of one column in the annotation table being queried.
}
  \item{notfound}{
The value to be used if that column is not found in the annotation table.
}
  \item{negate}{
Should the column of logical values be negated before evaluating 'which'.
}
  \item{verbose}{
report the details of the search
}
}
\details{
This is the principal means of querying sample specific values that affect the 
runtime behavior of many Duffy Lab package tools.  For settings that are not 
sample specific, see \code{\link{OptionsTable}}.  

The annotation table is a text file with tab separated fields.  Each row represents
one sample or dataset or experiment, where the first (leftmost) column is typically
labeled 'SampleID' and that column is the sample identifier.  The 
file may include a column called 'Include' or 'Exclude' that
controls which rows of the table are actually read in from the file.  This makes
it easy to control what subset of samples are seen by the analysis tools.

Any text in the annotation table file can have embedded environment variables,
which will be run-time substituted at the time they are 'gotten'.

The truth of a value uses the general heuristic that {F, False, N, No, 0, and the 
empty string ""} resolve to \code{FALSE} and everything else is \code{TRUE}. 
See \code{\link{as.TRUEorFALSE}}.
}
\value{
For \code{\link{readAnnotationTable}} a table of sample annotation values, with 'excluded' 
rows removed.

For \code{\link{getAnnotationValue}} the one entry from the table at:  row='key', col='arg'
as a characer string, or the 'notfound' value.

For \code{\link{getAnnotationTrue}} the one entry from the table at:  row='key', col='arg' 
as a logical value, or the 'notfound' value.

For \code{\link{whichAnnotationTrue}} a vector of row numbers where column 'arg' resolves 
to a logical TRUE value, after possible negation.
}
\author{
Bob Morrison
}
\seealso{
  \code{\link{readOptionsTable}},  \code{\link{getOptionValue}},  \code{\link{getOptionTrue}},
  for corresponding functions for program options that are not sample specific.
}
