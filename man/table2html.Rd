\name{table2html}
\alias{table2html}
\title{
Convert Data Frame to HTML file
}
\description{
Convert any data frame into a HTML file, augmented with hyperlinks to
arbitrary content.
}
\usage{
table2html(x, fileout, title = "", maxRows = NULL, extraHTMLtextFile = NULL, 
		linkColumnNames = "GENE_ID", linkPaths = ".", linkExtensions = ".png",
		linkTargetColumnNames=linkColumnNames)
}
\arguments{
  \item{x}{
a data frame to be turned to HTML
}
  \item{fileout}{
the full pathname to the newly created HTML file
}
  \item{title}{
the title for the HTML document.  May contain hypertext markup.
}
  \item{maxRows}{
the maximum number of rows of \code{x} to include in the resulting file.
}
  \item{extraHTMLtextFile}{
an optional text file on content, prepended to the HTML file prior to the
table
}
  \item{linkColumnNames}{
character vector of table column names, that will be converted to hyperlinks.  By default
the content of the cell is both the base filename of the created hyperlink and the text
displayed in the cell.
}
  \item{linkPaths}{
character vector of pathname prefixes that become part of the hyperlink anchor tag
}
  \item{linkExtensions}{
character vector of pathname suffixes that becomre part of the hyperlink anchor tag
}
  \item{linkTargetColumnNames}{
optional character vector of table column names, that allow the hyperlink for a cell to be built
from text other than what is displayed in the cell.  Most often used to make GeneID plot hyperlinks
for columns that are not gene identifiers.
}
}
\details{
This function is used to turn trancript and DE files into clickable HTML with links
to gene plots.  Care should be taken to assure that any target plots have the 
corrent filenames and path locations, to properly match the above \code{link...} 
arguments.  It is sometimes advisable to edit the column names and formatting
of \code{x} before turning it to HTML to improve the final appearance in the 
browser window.
}
\value{
A file of HTML is written to disk.
}
\references{
  uses much of the original CGIwithR routines by David Firth.
}
\author{
Bob Morrison
}
