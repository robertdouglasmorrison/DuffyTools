\name{movingAverage}
\alias{movingAverage}
\alias{unitIntervalMA}
\title{
Moving Average
}
\description{
Apply a moving average window as a smoothing function to a numeric vector
}
\usage{
movingAverage(x, window = 9, at = names(x), FUN = mean, do.ends = TRUE)

unitIntervalMA(x, step=1)
}
\arguments{
  \item{x}{
numeric vector of data to be smoothed by a moving window
}
  \item{window}{
size of the window, in elements of \code{x}
}
  \item{at}{
numeric vector representing the locations that the values of \code{x} come from
}
  \item{FUN}{
averaging function to apply the the windowed subsets of \code{x}
}
  \item{do.ends}{
logical, should the ends of \code{x} that can't have a full window be droppped (FALSE) 
or approximated with smaller windows (TRUE)
}
  \item{step}{
increment stepsize, passed to \code{seq( ..., by=step)}
}
}
\value{
A numeric vector, where the values are the averages in the window, and the names are
the locations for the centers of the windows used.

To force the result of \code{movingAverage} to have fixed gaps between 
\code{at}, it can be passed to \code{unitIntervalMA}, which uses splines 
to interpolate to a new set of locations with a constant step size.
}
