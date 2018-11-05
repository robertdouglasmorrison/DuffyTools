# elapsedProcTime.R

`elapsedProcTime` <- function( t1, t2, N=1, what="Read") {

	# t1, t2 are proc.time() objects, we want certain elements from them
	computeTerms <- c( 1,2,4)
	wallTerms <- 3

	computeTime <- sum( t2[ computeTerms] - t1[ computeTerms])
	wallTime <- sum( t2[ wallTerms] - t1[ wallTerms])
	perReadTime <- computeTime / N


	`applyTimeUnits` <- function( x) {
		ans <- vector()
		for (i in 1:length(x)) {
			if ( is.na(x[i])) 
				ans[i] <- "NA"
			else if ( x[i] < 0.000001)
				ans[i] <- paste( formatC( x[i]/0.000000001, digits=2, format="f"), " nanoSec")
			else if ( x[i] < 0.001)
				ans[i] <- paste( formatC( x[i]/0.000001, digits=2, format="f"), "microSec")
			else if ( x[i] < 1)
				ans[i] <- paste( formatC( x[i]/0.001, digits=2, format="f"), "milliSec")
			else if ( x[i] < 60)
				ans[i] <- paste( formatC( x[i]/1, digits=2, format="f"), " Seconds")
			else if ( x[i] < 3600)
				ans[i] <- paste( formatC( x[i]/60, digits=2, format="f"), " Minutes")
			else 
				ans[i] <- paste( formatC( x[i]/3600, digits=2, format="f"), "   Hours")
		}
		return( ans)
	}

	out <- applyTimeUnits( c( wallTime, computeTime, perReadTime))
	names( out) <- c( "Total.Elapsed.Time", "Computation.Time", paste( "Time.per", what, sep="."))
	out <- as.list( out)
	out$Raw.Seconds <- wallTime
	return( out)
}

