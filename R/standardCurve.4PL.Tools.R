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


`calculateStandardCurveABCD` <- function( odValues, concValues, plot=FALSE, label="", las=1, digits=4, 
					pt.col=4, lin.col=1, fit.col=4, add.plot=FALSE, show.legend=!add.plot) {

	require( minpack.lm)

	# prep the OD and Concentration values
	odValues <- as.numeric( odValues)
	concValues <- as.numeric( concValues)
	useV <- which( ! is.na( odValues))
	odValues <- odValues[useV]
	concValues <- concValues[useV]
	NCV <- length(concValues)
	if ( length(odValues)!= NCV || NCV < 5) {
		cat( "\nError: not enough OD & Concentration values to fit: ", NCV)
		return(NULL)
	}
	
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
		ord <- order( concValues, odValues, decreasing=FALSE)
		concValues <- concValues[ord]
		odValues <- odValues[ord]
		odPredicted <- odPredicted[ord]
		concFac <- factor( concValues)
		NCONC <- nlevels(concFac)
		xPtr <- tapply( 1:NCV, concFac, FUN=NULL)
		avgODobs <- tapply( odValues, concFac, mean)
		avgODpred <- tapply( odPredicted, concFac, mean)

		# we may have replicates of the concentrations.  So deduce how many unique values

		# either create a new plot, or just add data to existing
		if ( ! add.plot) {
			plot.default( xPtr, odValues, type="p", pch=1, col=1, xlab="Concentration", ylab="OD Values", 
					xaxt="n", main=paste("Standard Curve Fit: ", label))
			axis( side=1, at=1:NCONC, round( as.numeric(levels(concFac)), digits=digits), las=las)
		}
		lines( 1:NCONC, avgODobs, lty=1, col=lin.col)
		points( xPtr, odPredicted, pch=2, col=pt.col)
		lines( 1:NCONC, avgODpred, lty=2, col=fit.col)
		if ( show.legend) {
			legend( 'bottomright', c("Observed","Fitted"), pch=c(1,2), col=c(1,4), lwd=2)
			outNames <- sub( "(^[ABCD])", "    \\1", names(out2))
			legend( 'topleft', paste( outNames, "=", round(out2, digits=digits)))
		}
	}

	return( out)
}
