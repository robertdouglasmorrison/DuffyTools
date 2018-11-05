\name{setCurrentTarget}
\alias{setCurrentTarget}
\title{
Change the Current Target
}
\description{
Change the target species(s) for subsequent analysis, thereby
changing the set of species to be investigated.
}
\usage{
setCurrentTarget(targetID = NULL, optionsFile = NULL)
}
\arguments{
  \item{targetID}{
character string of one known TargetID
}
  \item{optionsFile}{
-or- the name of an \code{\link{OptionsTable}} to extract the
TargetID from
}
}
\details{
This function is the way to configure any of the Duffy Lab packages
to operate on a certain set of organisms.  See \code{\link{Targets}} for
full details.  If \code{optionFile} is specified, the option called
\code{targetID} is used as the value.
}
\value{
Returns the now-current TargetID, as a character string.  The main
effect of this routine is to alter the family of species for subsequent
analysis, by changing the information returned by the 
\code{getCurrentTarget...()} family of functions.
}
\examples{
  setCurrentTarget("MT_H37")
  getCurrentTargetSpecies()
  setCurrentTarget("HsPf")
  getCurrentTargetSpecies()
}
