# splineTools.R -- functions for best fit spline lines


spline.line <- function( x, y, n=length(x)*5, useLog=FALSE, ...) {

	if (useLog) {
		splineAns <- spline( x, log10(y+1), n=n)
		xSpline <- splineAns$x
		ySpline <- (10 ^ (splineAns$y)) - 1
		# spline in log space can look extreme...
		minY <- min(y, na.rm=T) * 0.1
		ySpline[ ySpline < minY] <- minY
	} else {
		splineAns <- spline( x, y, n=n)
		xSpline <- splineAns$x
		ySpline <- splineAns$y
	}

	lines( xSpline, ySpline, ...)
}
