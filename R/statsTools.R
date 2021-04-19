# statsTools.R -- various statistics functions


`cv` <- function( x, y=NULL, na.rm=FALSE) {


	# do either the 1-sample or 2-sample Coefficient fo Variance
	if ( is.null( y)) {

		# simple case
		return( sd( x, na.rm=na.rm) / mean( x, na.rm=na.rm))
	}

	# the 2-sample CV:  using the 'within subject standard deviation method (Jones & Payne, 1997)
	if ( length(y) != length(x)) stop( "X and Y must be of equal length")

	twoN <- ( 2 * length(x))
	myMean <- sum( x + y) / twoN
	mySD <- sqrt( sum( (x-y)^2) / twoN)

	return( mySD / myMean)
}


`sparse.t.test` <- function( x, y=NULL, min.obs=3, ..., min.random.value=NULL) {

	# make extra observed values using sensible constraints
	myX <- x[ ! is.na(x)]
	if ( length(myX) < min.obs) {
		meanX <- mean( myX)
		sdX <- min( sqrt( abs( meanX/2)), abs( meanX/2))
		extraX <- rnorm( min.obs, mean=meanX, sd=sdX)
		if ( ! is.null( min.random.value)) extraX[ extraX < min.random.value] <- min.random.value
		x <- c( x, extraX)
	}
	if ( is.null(y)) return( t.test( x, ...))

	myY <- y[ ! is.na(y)]
	if ( length(myY) < min.obs) {
		meanY <- mean( myY)
		sdY <- min( sqrt( abs( meanY/2)), abs( meanY/2))
		extraY <- rnorm( min.obs, mean=meanY, sd=sdY)
		if ( ! is.null( min.random.value)) extraY[ extraY < min.random.value] <- min.random.value
		y <- c( y, extraY)
	}
	return( t.test( x, y, ...))
}

