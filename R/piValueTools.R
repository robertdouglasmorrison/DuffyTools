# piValueTools.R -- various statistics functions
		
# turn the blend of Fold Change and P-value into a single measure, as defined by:
# Xiao Y, et.al.  Novel Significance Score for Gene Selection and Ranking 
#  		(2014) Bioinformatics 30(6):801-807

# with some optional modifications to limit effect of extreme small P-values

`piValue` <- function( log2fold, pvalue, max.correction=4) {

	# trap zeros that go undefined
	pvalue[ pvalue < 1e-100] <- 1e-100
	pvalue[ pvalue > 0.99] <- 0.99

	# do the minus log 10 step
	log10p <- -log10( pvalue)

	# before applying, catch unwanted effect of absurdly low P values
	if ( ! is.null( max.correction) && as.numeric(max.correction) >= 1) {
	    if ( any( log10p > max.correction)) {
		# apply the correction to any that are net improvements
		whoFix <- which( log10p > 1.0)
		temp10p  <- log10p[ whoFix]
		# apply a second log 10 transform to these correction terms, to better smooth the affect
		temp10.10p <- log10( temp10p)
		max.log10 <- log10( max.correction)
		scaleFac <- max.log10 / max(temp10.10p)
		# apply that scaling and under the second log10 step
		temp10.10p <- temp10.10p * scaleFac
		temp10p <- 10 ^ temp10.10p
		# ok, put these back
		log10p[whoFix] <- temp10p
	    }
	    # also correct the other way, to prevent really high P values from completely zeroing the FC
	    if ( any( log10p < 1/max.correction)) {
		# apply this correction to any that are net punishments
		whoFix <- which( log10p < 1.0)
		temp10p  <- log10p[ whoFix]
		# apply a second log 10 transform to these correction terms, to better smooth the affect
		temp10.10p <- log10( temp10p)
		min.log10 <- log10( 1/max.correction)
		scaleFac <- min.log10 / min(temp10.10p)
		# apply that scaling and under the second log10 step
		temp10.10p <- temp10.10p * scaleFac
		temp10p <- 10 ^ temp10.10p
		# ok, put these back
		log10p[whoFix] <- temp10p
	    }
	}

	# apply the pi value formula of Xiao
	out <- log2fold * log10p
	return( out)
}


`piValueRankOrder` <- function( log2fold, pvalue, max.correction=4) {

	piVal <- piValue( log2fold, pvalue, max.correction=max.correction)
	return( order( piVal, decreasing=TRUE))
}

