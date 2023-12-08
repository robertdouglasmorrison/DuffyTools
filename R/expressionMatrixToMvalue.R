# expressionMatrixToMvalule.R -- turn abundance into log 2 fold change


`expressionMatrixToMvalue` <- function( x, average.FUN=median, baselineColumns=NULL, 
					minIntensity=0, small.offset=1, verbose=T) {


	x <- as.matrix(x)

	# prep to make sure all values are within expected range
	if ( small.offset < 0) {
		cat( "'small.offset' cannot be negative.  Setting to zero..")
		small.offset <- 0
	}
	if ( minIntensity < 0) minIntensity <- 0
	if (any( x < minIntensity)) {
		if (verbose) cat( "Clipping low abundance values at: ", minIntensity)
		x[ x < minIntensity] <- minIntensity
	}

	# calculate the average expression for each gene row, allowing for use of a subset of columns
	if ( is.null(baselineColumns)) {
		rowAvgs <- apply( x, MARGIN=1, FUN=average.FUN, na.rm=T)
	} else {
		useColumns <- intersect( 1:ncol(x), as.integer(baselineColumns))
		if ( length(useColumns) < 1) {
			cat( "\nError: 'baselineColumns' must be integers in the range 1 to ncol(x)")
			stop( "invalid value for baseline column subsetting")
		}
		x2 <- x[ , useColumns, drop=F]
		rowAvgs <- apply( x2, MARGIN=1, FUN=average.FUN, na.rm=T)
	}

	# do the log2 ratio
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

