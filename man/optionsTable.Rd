\name{optionsTable}
\alias{OptionsTable}
\alias{readOptionsTable}
\alias{getOptionValue}
\alias{getOptionTrue}
\title{
Extract Facts from an Options Table
}
\description{
Functions to query an Options Table for values that affect program behaviour
}
\usage{
readOptionsTable(file)

getOptionValue(tbl, arg, notfound = NA, verbose = TRUE, speciesID = NULL)
getOptionTrue(tbl, arg, notfound = NA, verbose = TRUE, speciesID = NULL)
}
\arguments{
  \item{file}{
Name of an options file.
}
  \item{tbl}{
Either an options table previously read via 'readOptionsTable' or a character string 
file name to open and read from.
}
  \item{arg}{
The named option to be found in the table.
}
  \item{notfound}{
The value to be used if that option is not found in the options table,
where \code{NA} means the option is required to be in the table.
}
  \item{verbose}{
Report the details of the search to standard out.
}
  \item{speciesID}{
A SpeciesID for specifying an option that has different values for different organisms.
Treated as a suffix to \code{arg}, so if \code{speciesID='MTb'} and \code{arg='DE.minimumRPKM'},
the function will look for an option named \code{'DE.minimumRPKM.MTb'}.
}
}
\details{
This is the principal means of querying values that affect the general
runtime behavior of many Duffy Lab package tools.  For settings that are 
sample specific, see \code{\link{AnnotationTable}}.  

The options table is a text file with tab separated fields.  Each row represents
one named option, where the first (leftmost) column is the option name and
the second column is the value for that option.

Any text in the value column of the options table file can have embedded 
environment variables, which will be run-time substituted at the time they are 'gotten'.

The truth of a value uses the general heuristic that {F, False, N, No, 0, and the 
empty string ""} resolve to \code{FALSE} and everything else is \code{TRUE}. 
See \code{\link{as.TRUEorFALSE}}.
}
\value{
For \code{\link{readOptionsTable}} a table of option values, with 2 columns 
\code{OptionName, Value}.

For \code{\link{getOptionValue}} the one entry from the table at:  OptionName='arg'
as a characer string, or the 'notfound' value.

For \code{\link{getOptionTrue}} the one entry from the table at:  OptionName='arg' 
as a logical value, or the 'notfound' value.
}
\author{
Bob Morrison
}

\seealso{
  \code{\link{readAnnotationTable}},  \code{\link{getAnnotationValue}},  
\code{\link{getAnnotationTrue}}, for the corresponding functions for querying sample 
specific values from the Annotation table.
}
