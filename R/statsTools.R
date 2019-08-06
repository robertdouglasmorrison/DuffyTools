# statsTools.R -- various statistics functions


cv <- function( x, y=NULL, na.rm=FALSE) {


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
		
