`erfc` <-
function( x) { 
	sqrt2 <- 2^0.5
	return( 2 * pnorm( x * sqrt2, lower=FALSE))
}

