# enrichment.R

# calcualate the probability of seeing a certain number of genes in common
`enrichment` <- function( nMatch, nYourSet, nTotal, nTargetSubset) {

	x <- nMatch
	m <- nTargetSubset
	n <- nTotal - m
	k <- nYourSet

	expected <- (nTargetSubset / nTotal) * nYourSet
	enrich <- x / expected

	# get the entire prob. dist. and sum up both halves
	allPs <- dhyper( 0:k, m, n, k)
	# the probability of 0 to X genes
	lowerTailPvalue <- sum( allPs[1:(x+1)])
	# the probability of X up to K genes
	upperTailPvalue <- sum( allPs[(x+1):length(allPs)])

	out <- list( "nWhiteBalls"=m, "nBlackBalls"=n, "nDrawn"=k, "nWhiteDrawn"=x, 
			"nExpected"=expected, "Enrichment"=enrich, 
			"P_atLeast_N"=upperTailPvalue, "P_atMost_N"=lowerTailPvalue)

	return( out)
}


enrichment.Nway <- function( nSets=2, nMatch, nDrawn, nTotal, nSimulations=1000000) {

	nTotal <- as.integer( nTotal)
	nDrawn <- as.integer( nDrawn)

	tallyEnrichment <- function( nOverlap, nGood, nSimulation) {

		# trim storage and tabulate
		if ( length( nOverlap) > nGood) length( nOverlap) <- nGood
		dist <- table( nOverlap)
		expected <- sum( nOverlap) / nSimulation

		yes <-sum( nOverlap >= nMatch)
		no <- nSimulations - yes

		return( list( "nExpected"=expected, 
			"P_atLeast_N"=( yes / nSimulations), 
			"P_atMost_N"=( no / nSimulations), 
			"distribution"=dist))
	}

	# storage to hold how many in common from all successful trials
	nOverlap <- vector( mode="numeric", length=nSimulations/10)
	nGood <- 0

	# do those trials
	base::lapply( 1:nSimulations, function(x) {
		picks <- base::unlist( base::lapply( 1:nSets, FUN=function(x) sample.int( n=nTotal, size=nDrawn)))
		hits <- sum( tabulate(picks) == nSets)
		if ( hits > 0) {
			nGood <<- nGood + 1
			nOverlap[ nGood] <<- hits
		}

		if ( x %% 10000 == 0) {
			ans <- tallyEnrichment( nOverlap, nGood, x)
			cat( "\rIter:", x, "  Expect:", ans$nExpected, "  P.atleastN:", ans$P_atLeast_N,
					"  P.atmostN:", ans$P_atMost_N,
					"  Dist: ", paste( names( ans$distribution), ans$distribution, sep="-",
					collapse=", "))
		}
		return()
	})

	# trim storage and tabulate
	return( tallyEnrichment( nOverlap, nGood, nSimulations))
}


simulate.enrichment.Nway <- function( nSets=3, nMatch=40, nDrawn=100, nTotal=5475, nSimulations=100000) {

	cat( "\nChoose", nDrawn, "from a pool of", nTotal, "into", nSets, "sets.")
	cat( "\nLikelihood of finding ", nMatch, "in common:")
	cat( "\nSimulating ", nSimulations, " random trials...")
	ans <- enrichment.Nway( nSets=nSets, nMatch=nMatch, nDrawn=nDrawn, nTotal=nTotal, nSimulations=nSimulations)
	cat( "\n")

	expect <- ans$nExpected
	hits <- ans$distribution
	ats <- as.numeric( names( hits))
	freqs <- as.numeric( hits)
	# add the number of 'zero' events into the model
	nHits <- sum( hits)
	if ( nHits < nSimulations) {
		nZero <- nSimulations - nHits
		ats <- c( 0, ats)
		freqs <- c( nZero, freqs)
	}
	if ( length( ats) < 4) cat( "\nMay not have enough points for NLS curve fitting...")
	
	# we need the entire interval from zero to 'nDrawn'
	xEvents <- 0:nDrawn
	normAns <- fit.normal( x=ats, y=freqs)
	yEvents <- normal( xEvents, mean=normAns$mean, sd=normAns$sd, height=normAns$height)

	# use a finer mesh for the drawn fit curve
	xShowMax <- max( nMatch+1, ats)
	xShow <- seq( 0, xShowMax, by=0.2)
	yShow <- normal( xShow, mean=normAns$mean, sd=normAns$sd, height=normAns$height)
	plot( xShow, yShow, main=paste( "Simulation Curve of 'Choose ", nDrawn, " from ",nTotal,"' overlap in ",nSets," datasets"), 
		xlab=paste("Number of genes found in 'Top ",nDrawn,"' of all ",nSets," datasets"), 
		ylab= "Number of Occurrances", type="l", lwd=2, col=2, log='')
	points( ats, freqs, cex=1.2, lwd=2, col=1)
	points( nMatch, yEvents[ xEvents == nMatch], pch=21, cex=1.8, col=1, bg=5)

	# we report area under the full probability curve
	pcts <- yEvents / sum( yEvents)
	who <- which( xEvents < nMatch)
	if ( nMatch < expect) {
		myP <- sum( pcts[who])
	} else {
		myP <- 1.0 - sum( pcts[who])
	}

	out <- c( ans, "fitted.p.value"=myP)
	return( out)
}

