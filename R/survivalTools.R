# survivalTools.R -- implement Kaplan-Meier curves

kaplanMeier <- function( groups, times, outcomes, eventOutcome="Y", col=2:(length(unique(groups))+1),
			makePlot=TRUE, xlab="Time  (weeks)", ylab="Survival", lwd=3, legend.cex=1.1,
			main="Kaplan-Meier: ", nFDR=0) {

	require( survival)

	# turn the outcome calls and followup time into censored survival data
	surv <- Surv( time=times, event=(outcomes == eventOutcome))
	survAns <- survfit( surv ~ groups)
	survDiff <- survdiff( surv ~ groups)

	# calc the P-value
	pval <- pchisq( survDiff$chisq, df=1, lower=F)

	# gather details for plotting
	grpNames <- sort( unique( groups))

	if (makePlot) {
		plot( survAns, col=col, lwd=lwd, xlab=xlab, ylab=ylab, main=main,
			ylim=c(0,1.03))
		legend( 'topright', grpNames, col=col, lwd=lwd, bg='white', cex=legend.cex)
		legend( 'bottomleft', paste( "P-value =", round( pval, digits=4)), bg='white', cex=legend.cex)
	}

	# perhaps doa FDR permutation test
	if (nFDR > 0) {
		pFDR <- vector( length=nFDR)
		N <- length(groups)
		for ( i in 1:nFDR) {
			fdr.groups <- groups[ sample(N)]
			survDiff2 <- survdiff( surv ~ fdr.groups)
			pFDR[i] <- pchisq( survDiff2$chisq, df=1, lower=F)
		}
		myFDR <- sum( pFDR <= pval) / nFDR
		legend( 'bottomright', paste( "FDR (Trials=",nFDR,") = ", round( myFDR, digits=4), sep=""), 
			bg='white', cex=legend.cex)
	}

	out <- list( "surv.diff"=survDiff, "p.value"=pval)
	return( invisible( out))
}

