# smfaTools.R -- low level SMFA helper functions

# functions to calculate the main 4 SMFA metrics:  TRA and TBA, etc.

# Notes from Jen Hume
# For us the equations that Olga uses are:
# TRA = (mean control oocyst number – mean test oocyst number) * 100 / mean control oocyst number
# TBA = (control % infection prevalence - test % infection prevalence) * 100 / % infection prevalence
SIG_DIGITS <- 5
SMFA.Metric.Names <- c( "Avg.Oocysts", "N.Dissected", "N.Zero.Oocysts", "Infect.Prevalence", "TRA",  "TBA",
			"Ctrl.Avg.Oocysts", "Ctrl.Prevalence")

nDissected <- function( oocystCounts) { return( sum( ! is.na( oocystCounts)))}

avgOocysts <- function( oocystCounts) { return( mean( as.numeric( oocystCounts), na.rm=T))}

nZeroOocysts <- function( oocystCounts) { return( sum( oocystCounts == 0, na.rm=T))}

prevalence <- function( oocystCounts) {
		ndiss <- nDissected( oocystCounts)
		nzero <- nZeroOocysts( oocystCounts)
		return( round( (ndiss - nzero) * 100 / ndiss, digits=SIG_DIGITS))
}

# Calculate TRA using a model instead of just means, to give confidence intervals
tra <- function( testOocystCounts, controlOocystCounts, statistics=FALSE) {

	# remove any invalid data
	testV <- as.numeric(testOocystCounts)
	testV <- testV[ ! is.na(testV)]
	ctrlV <- as.numeric(controlOocystCounts)
	ctrlV <- ctrlV[ ! is.na(ctrlV)]
	# set up as a model formula
	ctrl <- mean(ctrlV)
	test <- mean(testV)
	tra <- round( 100 * ( 1 - (test/ctrl)), digits=SIG_DIGITS)
	if ( ! statistics) return( tra)
	
	# implement the CI method 2 (WMW) from Swihart J Am Stats 2018
	testV[ testV == 0] <- 0.5
	ctrlV[ ctrlV == 0] <- 0.5
	wAns <- suppressWarnings( wilcox.test( log(testV), log(ctrlV), conf.int=T))
	ciBetas <- wAns$conf.int[2:1]
	traCI <- round( 100 * ( 1 - exp(ciBetas)), digits=SIG_DIGITS)
	out <- list( "tra"=tra, "confint"=traCI, "p.value"=wAns$p.value)
	return( out)
}

tba <- function( testOocystCounts, controlOocystCounts) {
	prevCntrl <- prevalence( controlOocystCounts)
	prevTest <- prevalence( testOocystCounts)
	ans <- (prevCntrl - prevTest) * 100 / prevCntrl
	return( round( ans, digits=SIG_DIGITS))
}

# do all ~8 metrics we want from one feed
smfaValues <- function( testOocystCounts, controlOocystCounts) {
	myAvgOocysts <- round( avgOocysts( testOocystCounts), digits=SIG_DIGITS_FINAL)
	myTRA <- round( tra( testOocystCounts, controlOocystCounts), digits=SIG_DIGITS_FINAL)
	myNdissect <- nDissected( testOocystCounts)
	myZeroOocyst <- nZeroOocysts( testOocystCounts)
	myPrevalence <- round( prevalence( testOocystCounts), digits=SIG_DIGITS_FINAL)
	myTBA <- round( tba( testOocystCounts, controlOocystCounts), digits=SIG_DIGITS_FINAL)
	ctrlAvgOocysts <- round( avgOocysts( controlOocystCounts), digits=SIG_DIGITS_FINAL)
	ctrlPrevalence <- round( prevalence( controlOocystCounts), digits=SIG_DIGITS_FINAL)
	
	out <- c( myAvgOocysts, myNdissect, myZeroOocyst, myPrevalence, myTRA, myTBA, ctrlAvgOocysts, ctrlPrevalence)
	names(out) <- SMFA.Metric.Names
	return(out)
}

# do all 8 metrics we want for all assays as a set
smfaValuesMatrix <- function( countsM, control) {
	N <- ncol(countsM)
	out <- matrix( NA, nrow=N, ncol=length(SMFA.Metric.Names))
	colnames(out) <- SMFA.Metric.Names
	rownames(out) <- colnames(countsM)
	for ( i in 1:N) out[ i, ] <- smfaValues( countsM[,i], countsM[,control])
	return(out)
}


oocystCountString <- function( oocystCounts) {
	countText <- as.character( oocystCounts)
	countText[ is.na( oocystCounts)] <- ""
	return( paste( countText, collapse=","))
}


# functions for doing IC80 regression fitting of SMFA and ELISA assay data

logMeanRatio <- function( tra) {

	# implement the LMR math of Kazutoyo Miura (Vaccine 2017)
	
	# step one: catch any 100% TRA values, as they break the log math
	bad <- which( tooLarge <- (tra >= 100))
	if ( length(bad)) {
		bigValid <- max( tra[ ! tooLarge], na.rm=T)
		smallDif <- 100 - bigValid
		imputedValue <- 100 - (smallDif/2)
		tra[ bad] <- imputedValue
	}
	
	# now do the log ratio step
	LMR <- log10( 100 / (100 - tra))
	return( LMR)
}


# do the correlation of TRA to ELISE EU, and perhaps plot the data as well
correlate.TRA.to.ELISA <- function( tra, eu, group=NULL, PLOT=TRUE, 
				col='black', bg='white', pch=1, icThreshold=80, 
				lwd=2, lty=2, legend.loc="topleft", label="Correlation of TRA & ELISA", 
				showIC=TRUE, ...) {

	# implement the math and figure style of Kazutoyo Miura (Vaccine 2017)
	N <- length(tra)
	if ( length(eu) != N) stop( "'tra' and 'eu' vectors must be same length")
	
	# allow doing correlations by subsets
	Ngrp <- 1
	if ( ! is.null(group)) {
		if ( length(group) != N) stop( "'group' must be same length as value vectors")
		if ( ! is.factor(group)) group <- factor(group)
		Ngrp <- nlevels(group)
		grpNames <- levels(group)
	} else {
		group <- rep.int( 1, N)
	}
	
	if (PLOT) {
		if ( length(col) < Ngrp) col <- rep( col, length.out=Ngrp)
		if ( length(bg) < Ngrp) bg <- rep( bg, length.out=Ngrp)
		if ( length(pch) < Ngrp) pch <- rep( pch, length.out=Ngrp)
	}

	# step 1: transform the input values to the desired units, after discarding any NA data
	naTRA <- which( is.na(tra))
	naEU <- which( is.na(eu))
	drops <- sort( union( naTRA, naEU))
	if ( length(drops)) {
		tra <- tra[ -drops]
		eu <- eu[ -drops]
		group <- group[ -drops]
	}
	lmr <- logMeanRatio( tra)
	sqrtEU <- sqrt( eu)
	
	# step 2: do the linear regression for each group
	icOut <- pvOut <- vector( length=Ngrp)
	for ( ig in 1:Ngrp) {
		use <- which( as.numeric(group) == ig)
		
		# do the regression for this subset of data
		smlLM <<- lm( lmr[use] ~ sqrtEU[use])
		smlCoef <- coef( summary( smlLM))
		b0 <- smlCoef[ 1, 1]
		b1 <- smlCoef[ 2, 1]
		pvOut[ig] <- smlCoef[ 2, 4]
		icValue <- ( ( log10( 100/(100-icThreshold)) - b0) / b1) ^ 2
		icOut[ig] <- icValue
		
		# draw the data
		if (PLOT) {
			oldMAI <- par("mai")
			if ( oldMAI[4] < 1) par( "mai"=c( oldMAI[1:3], 1))
			if ( ig == 1) {
				plot( 1, 1, type="n", xlab="ELISA EU  (sqrt scale)", ylab="% TRA", xlim=range(sqrtEU)*c(0.8,1.1),
						main=label, ylim=range( lmr, 0, 2.5), xaxt="n", yaxt="n", ...)
				elisaAts <- c( 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
				axis( side=1, at=sqrt(elisaAts), labels=elisaAts, ...)
				axis( side=4, at=pretty(lmr), ...)
				mtext( "Log Mean Ratio (LMR)", side=4, line=2.5, cex=1.25)
				traAts <- c( 0, 50, 80, 90, 95, 99)
				lmrAts <- logMeanRatio( traAts)
				axis( side=2, at=lmrAts, labels=traAts, ...)
			}
			points( jitter(sqrtEU[use],amount=0.025), jitter(lmr[use],amount=0.025), pch=pch[ig], col=col[ig], bg=bg[ig], ...)
			abline( reg=smlLM, col=col[ig], lwd=lwd, lty=lty)
		}
	}
	
	# done with the math and plotting
	if (PLOT && Ngrp >= 1) {
		legendText <- grpNames
		if ( showIC) {
			moreText <- paste( "IC", round(icThreshold), " = ", round(icOut), "EU  (P=", formatC(pvOut,format="e",digits=2), ")", sep="")
			legendText <- paste( legendText, moreText)
		}
		legend( legend.loc, legendText, pch=pch, col=col, pt.bg=bg)
	}
	
	out <- list( "IC"=icOut, "p.value"=pvOut)
	return(out)
}
