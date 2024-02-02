# expressionMatrixToMvalule.R -- turn abundance into log 2 fold change


`expressionMatrixToMvalue` <- function( x, average.FUN=median, baselineColumns=NULL, 
					groupNames=NULL, minIntensity=0, small.offset=1, verbose=T) {


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
		# new idea: if we are given group membership for all columns, use those to get the 
		# average of each group first, and then take the average of the groups.
		# this allows very uneven group counts to better reflect differences between groups
		# instead of differences between samples
		if ( is.null(groupNames)) {
			rowAvgs <- apply( x, MARGIN=1, FUN=average.FUN, na.rm=T)
		} else {
			if ( length(groupNames) != ncol(x)) stop( "'groupNames' length must match 'ncol(x)'")
			grpFac <- factor(groupNames)
			nGrps <- nlevels(grpFac)
			nGenes <- nrow(x)
			grpAvgs <- matrix( 0, nGenes, nGrps)
			for ( i in 1:nGenes) {
				grpAvgs[ i, ] <- tapply( x[i,], grpFac, FUN=average.FUN, na.rm=T)
			}
			rowAvgs <- apply( grpAvgs, MARGIN=1, FUN=average.FUN, na.rm=T)
		}
	} else {
		# expecting column numbers, but allow column names
		if ( is.character(baselineColumns)) {
			useColumns <- match( baselineColumns, colnames(x))
			useColumns <- useColumns[ !is.na(useColumns)]
		} else {
			useColumns <- intersect( 1:ncol(x), as.integer(baselineColumns))
		}
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

