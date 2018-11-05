\name{tempFolder}
\alias{tempFolder}
\alias{removeTempFolder}
\title{
Create a Temporary Folder on an available file system tree
}
\description{
Searches for a named file system, and creates a temporary folder
}
\usage{
tempFolder(foldername, tmpRoot = "/tmp")
removeTempFolder(path)
}
\arguments{
  \item{foldername}{
The wanted leaf node name for the created folder path
}
  \item{tmpRoot}{
A character string of length 1, containing 1 or more possible file system temporary storage paths, separated by commas
}
  \item{path}{
The full pathname to the temp folder to be removed
}
}
\details{
This tries to allow creation of temporary storage in a OS independant manner, by allowing a multiplicity of 
named temporary storage paths.  The first 'found' path is selected.  No attempt to assess the quantity of 
free storage space is attempted.
}
\value{
For \code{tempFolder()}, a full pathname to a valid file system location, after creating that folder and verifying its creation.

For \code{removeTempFolder()}, no return value, but the named folder subtree and all its contents is deleted.
}
