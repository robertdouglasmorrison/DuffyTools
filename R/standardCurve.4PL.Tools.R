# standardCurve.4PL.Tools.R -- 
	
# very low level functions that implement the 4-parameter equation:
#   calculate OD (Y) given Concentration (X);  
#   or the inverse of calculate Concentration (X) given OD (Y);
#   called by nls() for doing curve fitting

# Using common definitions of the A/B/C/D variables.  These may need to be adjusted.
# A = lowest asymptote (low)
# B = slope at inflexion midpoint (Hill coefficient)
# C = inflexion point (midpoint)
# D = highest asymptote (high)


`forward.4pl` <- function(x, A, B, C, D) { return( D + ((A - D) / (1 + (x / C) ^ B)))}

`inverse.4pl` <- function( y, A, B, C, D) { return( (((1/(y-D) * (A-D)) - 1) ^ (1/B)) * C)}


`calculateStandardCurveABCD` <- function( odValues, concValues, plot=FALSE, label="", las=1, digits=4) {

	require( minpack.lm)

	odValues <- as.numeric( odValues)
	concValues <- as.numeric( concValues)
	useV <- which( ! is.na( odValues))
	odValues <- odValues[useV]
	concValues <- concValues[useV]
	NCV <- length(concValues)
	if ( length(odValues)!= NCV || NCV < 4) stop( "Error: not enough OD & Concentration values to fit")
	
	# we can now do the 4PL fit to deduce the best A/B/C/D variables.  
	# A = lowest asymptote (lowest OD)
	# B = slope at inflection midpoint (Hill's coefficient)
	# C = inflection point ('midpoint' along X-axis)
	# D = highest asymptote (highest OD)
	
	# start from some guesses about starting points and lower & upper bounds. These may need to be adjusted.
	A = 0.05
	B = 2
	C = median( concValues)
	D = 4
	starts <- c( "A"=A, "B"=B, "C"=C, "D"=D)
	lowers <- c( "A"=0.00001, "B"=0.00001, "C"=0.00001, "D"=1)
	uppers <- c( "A"=2, "B"=10, "C"=10, "D"=6)
	
	nlsAns <- nlsLM( odValues ~ forward.4pl( concValues, A, B, C, D), algorithm="LM", start=starts,
			lower=lowers, upper=uppers)
	
	nlsAns2 <- coef( summary( nlsAns))
	out1 <- round( nlsAns2[ , "Estimate"], digits=6)
	
	# also use this model to predict the OD values we should have seen for our inputs
	odPredicted <- predict( nlsAns)
	r2 <- cor( odValues, odPredicted) ^ 2
	
	out2 <- c( out1, "R^2"=r2)	
	out <- list( "estimate"=out2, model=nlsAns)
	
	if ( plot) {
		# make the plot have the expected direction, with high OD at the right
		ord <- order( odValues, decreasing=FALSE)
		plot.default( 1:NCV, odValues[ord], type="p", pch=1, col=1, xlab="Concentration", ylab="OD Values", 
				xaxt="n", main=paste("Standard Curve Fit: ", label))
		axis( side=1, at=1:NCV, round( concValues[ord], digits=digits), las=las)
		lines( 1:NCV, odValues[ord], lty=1, col=1)
		points( 1:NCV, odPredicted[ord], pch=2, col=4)
		lines( 1:NCV, odPredicted[ord], lty=2, col=4)
		legend( 'bottomright', c("Observed","Fitted"), pch=c(1,2), col=c(1,4), lwd=2)
		outNames <- sub( "(^[ABCD])", "    \\1", names(out2))
		legend( 'topleft', paste( outNames, "=", round(out2, digits=digits)))
	}

	return( out)
}
