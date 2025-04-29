# ROC_Tools.R - generic ROC function


duffy.ROC <- function( goodScores, badScores, label="ROC curve: ", visualize=TRUE, legend.cex=1, label.cex=1,
			show.cut.value=FALSE) {

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
	sensAns <- sens[best]
	specAns <- spec[best]

	if (visualize) {
		#par( mai=c(1,1,0.8,0.4))
		plot( mspec, msens, main=label, xlim=c(0,1.1), ylim=c(0,1.04),
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
				digits=ndigShow), pos=posOffset[1], srt=rotateAngle[1], offset=1.0, cex=label.cex)
		text( mspec[show[isLess]], msens[show[isLess]], formatC( cuts[show[isLess]], format="f", 
				digits=ndigShow), pos=posOffset[2], srt=rotateAngle[2], offset=1.5, cex=label.cex)

		# Patrick had an ask, to know how many observations land at a single location on the ROC curve.
		# for now focus on the upper right corner where the smallest "bad" value is represented.
		# But how we ask depends on which direction "bad" went..
		if ( mean(goodScores) > mean(badScores)) {
			nMostBad <- sum( badScores == min(badScores)) 
		} else {
			nMostBad <- sum( goodScores == min(goodScores)) 
		}
		pctChange <- ((msens[1]-msens[2]) / max(mspec[1]-mspec[2], 0.1)) 
		mySRT <- round( pctChange * 45)
		if ( nMostBad > 4) text( mean(mspec[1:2]), mean(msens[1:2]), paste( "Nobs = ",nMostBad,sep=""), pos=3, srt=mySRT, offset=0.75, cex=label.cex*0.75)

		# show the curve and the diagonal
		lines( c(0,1), c(0,1), col=1, lty=3)
		lines( mspec, msens, col=2, lwd=2)
	
		points( mspec[best], msens[best], pch=13, cex=2.5, col=2, lwd=2)
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
		if (show.cut.value) {
			myPos <- if (mspec[best] > 0.2) 2 else 3
			text( mspec[best], msens[best], paste("Cut =",cutAnsText,sep=" "), cex=label.cex*1.2, col=1, 
				pos=myPos, offset=1.05)
		}
		aucAnsText <- round( aucAns, digits=3)
		sensAnsText <- round( sensAns, digits=3)
		specAnsText <- round( specAns, digits=3)

		legendText <- paste( c( "AUC (Area Under Curve) =", "Optimal Yes/No Cutpoint =", "Sensitivity =", "Specificity ="),
				c( aucAnsText, cutAnsText, sensAnsText, specAnsText), "  ")
		legend( "bottomright", legendText, bg='white', cex=legend.cex)
		dev.flush()
	}
	
	return( list( "cutpoint"=cutAns, "AUC"=aucAns, "sensitivity"=sensAns, "specificity"=specAns))
}

