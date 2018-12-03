\name{BScountToPctParasitemia}
\alias{BScountToPctParasitemia}
\alias{pctParasitemiaToBScount}
\alias{BScountPerML}
\title{
Conversion Functions for Parasite Infection Measurements.
}
\description{
Conversion routines for quantifying parasite density in any of several standard methods,
as measured by assays like thin or thick smear by microscopy.
}
\usage{
BScountToPctParasitemia(BScount, WBCdenominator = 300, WBCperML = 8000, RBCperML = 5e+06)
pctParasitemiaToBScount(pctParasitemia, WBCdenominator = 300, WBCperML = 8000, RBCperML = 5e+06)
BScountPerML(BScount, WBCdenominator = 300, WBCperML = 8000)
}
\arguments{
  \item{BScount}{
Vector of parasite densities, given as the number of parasites per fixed number of white blood cells (WBC).
}
  \item{pctParasitemia}{
Vector of parasite densities, given as the percentage of red blood cells that are infected.
}
  \item{WBCdenominator}{
The fixed number of WBC used as the denominator when giving parasite density as a blood smear (BS) count.
}
  \item{WBCperML}{
The number of WBC per millileter (ML) of blood.  A value taken from an explicit assay such as 
a CBC machine is prefered.  Dufault is a typical average count.
}
  \item{RBCperML}{
The number of red blood cells (RBC) per millileter of blood.  A value taken from an explicit assay such as 
a CBC machine is prefered.  Dufault is a typical average count.
}
}
\details{
The functions are written such that CBC machine values can be used directly, where any invalid or missing
CBC values will be replace by default values.
}
\value{
A vector of converted parasite density measurements in the desired units.
}
\author{
Bob Morrison
}
