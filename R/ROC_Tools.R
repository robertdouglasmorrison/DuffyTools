# ROC_Tools.R - generic ROC function


duffy.ROC <- function( goodScores, badScores, label="", visualize=TRUE, legend.cex=1) {

	# to test ROC, we need the known positives and negatives
	require( ROC)

	# trap/remove any NAs
	goodScores <- goodScores[ ! is.na(goodScores)]
	badScores <- badScores[ ! is.na(badScores)]

	# build the truth table and our measurements
	truth <- c( rep( 1, times=length(goodScores)), rep( 0, times=length( badScores)))
	values <- c( goodScores, badScores)
	rule <- function( x, threshold) { ifelse( x > threshold, 1, 0)}

	# do the ROC call
	rocAns <<- rocdemo.sca( truth, values, rule, markerLabel=label, caseLabel="Positive Response")

	#if ( visualize) plot( rocAns, main=paste( "ROC curve:    ", label), line=F, show.thresh=T)
	
	# get the specificity and sensitivity
	spec <- rocAns@spec
	cuts <- rocAns@cuts
	sens <- rocAns@sens
	N <- length( spec)

	mspec <- (1 - spec)
	msens <- sens

	# best spot is the one whose distance to diagonal is maximal
	dist <- abs( mspec - msens)
	best <- which.max( dist)
	cutAns <- cuts[best]
	aucAns <- AUC( rocAns)

	if (visualize) {
		par( mai=c(1,1,0.8,0.4))
		plot( mspec, msens, main=paste( "ROC curve:    ", label), xlim=c(0,1.1),
			xlab="1 - Specificity", ylab="Sensitivity", type="p", font.axis=2, font.lab=2)

		# try to show some values in a sensible way
		show <- seq( 1, N, length.out=20)
		ndigits <- floor( log10( max( cuts)))
		ndigShow <- if ( ndigits > 3) 0 else (4-ndigits)
		isMore <- which( cuts[show] >= cutAns)
		isLess <- which( cuts[show] < cutAns)
		rotateAngle <- c( 0, -90)
		posOffset <- c( 4, 1)
		text( mspec[show[isMore]], msens[show[isMore]], formatC( cuts[show[isMore]], format="f", 
				digits=ndigShow), pos=posOffset[1], srt=rotateAngle[1], offset=1.0)
		text( mspec[show[isLess]], msens[show[isLess]], formatC( cuts[show[isLess]], format="f", 
				digits=ndigShow), pos=posOffset[2], srt=rotateAngle[2], offset=1.5)

		# show the curve and the diagonal
		lines( c(0,1), c(0,1), col=1, lty=3)
		lines( mspec, msens, col=2, lwd=2)
	
		points( mspec[best], msens[best], pch=16, cex=2.5, col=4)
		diagPt <- (msens[best]+mspec[best])/2
		lines( c(diagPt,mspec[best]), c(diagPt,msens[best]), col=2, lty=3)

		if ( cutAns < 10) ndigShow <- 3
		if ( cutAns < 1) ndigShow <- 5
		formatType <- "f"
		if ( cutAns < 0.01) {
			ndigShow <- 2
			formatType <- "e"
		}
		cutAnsText <- formatC( cutAns, format=formatType, digits=ndigShow)
		aucAnsText <- round( aucAns, digits=4)

		legend( "bottomright", paste( c( "AUC (Area Under Curve) =", "Optimal Yes/No Cut Point ="),
				c( aucAnsText, cutAnsText)), bg='white', cex=legend.cex)
	}
	
	return( list( "cutpoint"=cutAns, "AUC"=aucAns))
}

