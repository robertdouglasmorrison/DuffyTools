\name{averageLineWithErrorBars}
\alias{averageLineWithErrorBars}
\title{
Add error bars to an existing plot
}
\description{
Add an average line and standard error bars to an existing plot
}
\usage{
averageLineWithErrorBars(x, y, average.FUN = mean, col = 1, error.col = 1, error.lty = 1, whisker = 0.2, ...)
}
\arguments{
  \item{x}{
vector of X axis locations
}
  \item{y}{
vecotr of Y axis values, to be averaged at each value of X
}
  \item{average.FUN}{
The function for calculating the average
}
  \item{col}{
Color for the average line
}
  \item{error.col}{
Color for the error bars
}
  \item{error.lty}{
Line type for the error bars
}
  \item{whisker}{
Size of the whiskers, in units of X
}
  \item{\dots}{
Additional arguments, passed to \code{\link{line}}
}
}
\details{
A rather hokey attempt to be able to superimpose error bars on a variety of plots.
}
