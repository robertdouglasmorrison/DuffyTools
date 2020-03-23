# expressionMatrixToMvalule.R -- turn abundance into log 2 fold change


`expressionMatrixToMvalue` <- function( x, average.FUN=median, minIntensity=0, small.offset=1, verbose=T) {


	x <- as.matrix(x)

	if ( small.offset < 0) {
		cat( "'small.offset' cannot be negative.  Setting to zero..")
		small.offset <- 0
	}
	if ( minIntensity < 0) minIntensity <- 0
	if (any( x <= minIntensity)) {
		if (verbose) cat( "Clipping low abundance values at: ", minIntensity)
		x[ x < minIntensity] <- minIntensity
	}

	rowAvgs <- apply( x, MARGIN=1, FUN=average.FUN, na.rm=T)
	isZero <- which(rowAvgs == 0)
	mv <- x
	lapply( 1:nrow(x), function(i) {
			if (i %in% isZero && small.offset == 0) {
				mv[ i, ] <<- 0
			} else {
				mv[ i, ] <<- log2( (x[i, ] + small.offset) / (rowAvgs[i] + small.offset))
			}
			return()
		})

	return( mv)
}

