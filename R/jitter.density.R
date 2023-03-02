# jitterDensity.R -- implement a variant of jitter that uses the density distribution 
#			of a second variable to scale the magnitude of the jittered range.

`jitter.density` <- function( x, y, factor=1, amount=NULL, density.exp=NULL, 
			n.density=128, cut.density=0, bw.density="nrd0") {
	 
	# assess the density distribution separately for each unique group of X values
	N <- length(x)
	if ( ! N) return(x)
	xFac <- factor( x)
	nGroups <- nlevels(xFac)
	xPtrs <- tapply( 1:length(x), xFac, NULL)
	notNA <- which( !is.na(x) & !is.na(y))

	## set up storage for the result, and do each group separately
	xOut <- x
	for ( i in 1:nGroups) {
		use <- intersect( which( xPtrs == i), notNA)
		myX <- x[use]
		myY <- y[use]
		nNow <- length(myX)
		if ( nNow < 5) {
			xOut[use] <- jitter( myX, factor=factor, amount=amount)
			next
		}
		# calculate the desity of the Y values
		tmp.dens <- density( x=myY, bw=bw.density, cut=cut.density, n=n.density)
		xvals <- tmp.dens$x
		yvals <- tmp.dens$y

		# perhaps soften the shape of the curve, and then scale to unit max
		if ( ! is.null( density.exp)) yvals <- yvals ^ density.exp
		yvals <- yvals/max(yvals)
		
		# find the closest density X point to each given Y value
		# and use that density Y point as the factor to jitter our original X value
		bestDensPt <- sapply( myY, function(yy) which.min(abs(xvals-yy)))

		# use that height of the density curve at that location as the 
		# width of the random interval
		upLim <- yvals[ bestDensPt] * factor
		lowLim <- -upLim
		myJit <- runif( nNow, lowLim, upLim)
		xOut[use] <- myX + myJit
	}
	return( xOut)
}

