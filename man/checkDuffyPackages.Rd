\name{checkDuffyPackages}
\alias{checkDuffyPackages}
\title{
Check for Packages used by DuffyTools 
}
\description{
Interogates the current R installation for packages used by DuffyTools, 
trying to install those not yet detected.
}
\usage{
checkDuffyPackages()
}
\note{
Makes use of several Bioconductor packages, which it tries to install via the biocLite() function.
}
