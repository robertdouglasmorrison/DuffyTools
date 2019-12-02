# chromatogramModeling.R -- pieces to do modeling of chromatograms


`standardizeChromatogram` <- function( chromoObj, peak.dist=11, call.Ns=TRUE, constant.height=FALSE) {

	# turn the chromatogram into a uniform set of fixed width peaks, so we can add/model between multiple chromatograms
	peak.dist <- as.integer( peak.dist)
	half.width <- floor( peak.dist/2)

	traceIn <- chromoObj$TraceM
	peaksIn <- chromoObj$PeakPosition
	NP <- length( peaksIn)
	callsIn <- names(peaksIn)

	# make new storage to hold these peaks
	traceSizeOut <- NP * peak.dist
	traceOut <- matrix( 0, nrow=traceSizeOut, ncol=4)
	colnames(traceOut) <- colnames(traceIn)
	# hard code their new true center locations
	peaksOut <- seq( half.width+1, traceSizeOut, by=peak.dist)
	names(peaksOut) <- callsOut <- callsIn

	obsHeight <- vector()

	# we we will step along, and interpolate as needed to fill in the new data
	for ( i in 1:NP) {
		# make the set of indices in both systems for the region centered on one peak, 
		# and out +/- 1 peak in both directions.  Units are center peak at zero.
		centerIn <- peaksIn[i]
		leftIn <- if (i > 1) peaksIn[i-1] else 1
		rightIn <- if (i < NP) peaksIn[i+1] else nrow(traceIn)
		pointsIn <- (leftIn : rightIn) - centerIn
		centerOut <- peaksOut[i]
		leftOut <- if (i > 1) peaksOut[i-1] else 1
		rightOut <- if (i < NP) peaksOut[i+1] else nrow(traceOut)
		pointsOut <- (leftOut : rightOut) - centerOut

		# when these are perfect match, no interpolation is needed
		if ( length(pointsIn) == length(pointsOut)) {
			for ( j in 1:4) {
				vIn <- traceIn[ leftIn:rightIn, j]
				traceOut[ leftOut:rightOut, j] <- pmax( traceOut[leftOut:rightOut,j] , vIn)
			}
		} else {
			# interpolate
			scaleLeft <- min(pointsOut) / min(pointsIn)
			scaleRight <- max(pointsOut) / max(pointsIn)
			xxScaledLeft <- pointsIn[ pointsIn < 0] * scaleLeft
			xxScaledRight <- pointsIn[ pointsIn > 0] * scaleRight
			xx <- c( xxScaledLeft, 0, xxScaledRight)
			# may be small chance of raondoff error at the edges...  Force perfect
			xx[1] <- round( xx[1])
			xx[length(xx)] <- round( xx[length(xx)])
			# now we have what we need to do spline interpolation
			for ( j in 1:4) {
				vIn <- traceIn[ leftIn:rightIn, j]
				splineAns <- spline( xx, vIn, xout=pointsOut)
				traceOut[ leftOut:rightOut, j] <- pmax( traceOut[leftOut:rightOut,j] , splineAns$y)
			}
		}

		# tally the observed volume under each peak, in case we will scale to constant volume later
		obsHeight[i] <- sum( traceOut[ (centerOut-1):(centerOut+1), ])
	}

	# the Trace matrix is updated to uniform spacing for peaks.
	# We may be asked to update the un-called bases
	dna <- chromoObj$DNA_Calls
	aa <- chromoObj$AA_Calls
	isN <- which( callsIn == "N")
	if ( call.Ns && length(isN)) {
		anyChanges <- FALSE
		# visit each 'N' site, and assess the best base to call, keep just the first
		for ( i in isN) {
			thisPeakSite <- peaksOut[i]
			baseCallAns <- baseCallOnePeak( peakSite=thisPeakSite, traceM=traceOut)[1]
			if ( names(baseCallAns)[1] != "N") {
				anyChanges <- TRUE
				callsOut[i] <- names(baseCallAns)[1]
			}
		}
		# redo all the sequence info if any got updated
		if (anyChanges) {
			names(peaksOut) <- callsOut
			newSeq <- paste( callsOut, collapse="")
			seqData <- chromatogramSequences( newSeq)
			dna <- seqData$DNA
			aa <- seqData$AA
		}
	}

	# if asked to provide constant height, do that as a second pass
	if ( constant.height) {
		idealHeight <- median( obsHeight)
		minHeight <- idealHeight * 0.10
		trace2 <- traceOut
		for ( i in 1:NP) {
			centerOut <- peaksOut[i]
			left <- centerOut - half.width
			right <- centerOut + half.width
			# adjust the observed volume under each peak, by scaling all 4 bases at once
			scaleFac <- idealHeight / max( obsHeight[i], minHeight)
			tmpM <- traceOut[ left:right, ]
			tmpM <- tmpM * scaleFac
			trace2[ left:right, ] <- tmpM
		}
		traceOut <- trace2
	}

	# since we did interpolation on the trace matrix, clean it's significant digits some
	traceOut <- round( traceOut, digits=2)
	
	# lastly, recall the peak confidence
	peakConfidence <- calcChromatogramPeakConfidence( traceOut, peaksOut)

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=dna, "AA_Calls"=aa, "Filename"=chromoObj$Filename)
	out
}


`syntheticChromatogram` <- function( seq, peak.dist=11, height=1000, width=2.5, center=0, traceOnly=F) {

	# create a synthetic chromatogram, given a DNA sequence and some constants
	peak.dist <- as.integer( peak.dist)
	half.width <- floor( peak.dist/2)
	NP <- nchar( seq)
	calls <- strsplit( toupper(seq), split="")[[1]]

	# create storage for this synthetic chromatogram
	traceSize <- NP * peak.dist
	traceM <- matrix( 0, nrow=traceSize, ncol=4)
	colnames(traceM) <- c( "A", "C", "G", "T")
	peaks <- seq( half.width+1, traceSize, by=peak.dist)
	names(peaks) <- calls

	# make one gaussian peak, that extends about 2 peaks to either side
	peak.tail <- peak.dist * 2
	myX <- -peak.tail : peak.tail

	nHt <- length( height)
	nWd <- length( width)
	nCtr <- length( center)
	useVariableGauss <- ( nHt > 1 && nWd > 1 && nHt == nWd)
	fixedGaussV <- gaussian( myX, center=center[1], width=width[1], height=height[1], floor=0)
	
	# determine which column of the trace gets intensity for each base call
	columnHit <- match( calls, colnames(traceM), nomatch=0)
	
	# add this to every peak
	for ( i in 1:NP) {
		# get where this next peak is centered
		thisCenter <- peaks[i]
		# get the extent of its tails, and trim the gaussian fit to match
		thisLeft <- thisCenter - peak.tail
		thisRight <- thisCenter + peak.tail
		if (useVariableGauss) {
			thisV <- gaussian( myX, center=center[i], width=width[i], height=height[i], floor=0)
		} else {
			thisV <- fixedGaussV
		}
		NV <- length( thisV)
		if (thisLeft < 1) {
			nDrop <- 1 - thisLeft
			thisLeft <- 1
			thisV <- thisV[ (nDrop+1) : NV]
			NV <- length( thisV)
		}
		if (thisRight > traceSize) {
			nDrop <- thisRight - traceSize
			thisRight <- traceSize
			thisV <- thisV[ 1 : (NV-nDrop)]
			NV <- length( thisV)
		}
		# ready to add this in...   if the base is "N", give 1/4 to all
		thisColumn <- columnHit[i]
		if ( thisColumn > 0) {
			traceM[thisLeft:thisRight,thisColumn] <- traceM[thisLeft:thisRight,thisColumn] + thisV
		} else {
			thisV <- thisV * 0.25
			for (k in 1:4) traceM[thisLeft:thisRight,k] <- traceM[thisLeft:thisRight,k] + thisV
		}
	}

	# make the other fields the chromatogram has
	if (traceOnly) {
		dna <- aa <- ""
	} else {
		seqData <- chromatogramSequences( seq)
		dna <- seqData$DNA
		aa <- seqData$AA
	}

	# since we did interpolation on the trace matrix, clean it's significant digits some
	traceM <- round( traceM, digits=2)

	# lastly, call the peak confidence
	peakConfidence <- calcChromatogramPeakConfidence( traceM, peaks)
	
	out <- list( "TraceM"=traceM, "PeakPosition"=peaks, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=dna, "AA_Calls"=aa, "Filename"="")
	out
}


#`synthChromatogram` <- function( seq, height=1000, width=2.5, center=0) {
`syntheticTraceMatrix` <- function( seq, height=1000, width=2.5, center=0) {

	# bares bones wrapper around synthetic chromatogram, given a DNA sequence and some constants
	# returns just the trace matrix
	ans <- syntheticChromatogram( seq, height=height, width=width, center=center, traceOnly=TRUE)
	return( ans$TraceM)
}


# special function for using Generalized Simulated Annealing (GenSA) on chromatograms
# GenSA needs a function that calculates the error residual of the current model
`genSA.Chromatogram.residual` <- function( terms, obs, seq) {

	# either using a fixed Ht/Wd model, or separate for each peak
	if ( (N <- length(terms)) == 3) {
		ht <- terms[1]
		wid <- terms[2]
		ctr <- terms[3]
	} else {
		n3 <- round( N / 3)
		ht <- terms[ 1:n3]
		wid <- terms[ (n3+1):(n3*2)]
		ctr <- terms[ (n3*2+1):N]
	}

	ansModel <- syntheticTraceMatrix( seq, height=ht, width=wid, center=ctr)
	diff <- ( obs - ansModel)
	resid <- sqrt( sum( diff * diff))
	return( resid)
}


`subtractChromatogram` <- function( chromoObj1, chromoObj2) {

	# implenment as Object 1 minus Object 2

	# get the data we need 
	traceM1 <- chromoObj1$TraceM
	traceM2 <- chromoObj2$TraceM
	if ( ! all( dim(traceM1) == dim(traceM2))) {
		cat( "\nError:  Chromatograms must be identical sizes..")
		return( NULL)
	}

	# subtract the second trace matrix (model) from the first (observed)
	traceOut <- traceM1 - traceM2
	traceOut[ traceOut < 0] <- 0

	# for now, just keep all other data from Object1
	# at some point, we should fix this to recall bases, confidence, etc....
	peaksOut <- chromoObj1$PeakPosition
	peakConf <- chromoObj1$PeakConfidence
	dna <- chromoObj1$DNA_Calls
	aa <- chromoObj1$AA_Calls

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "PeakConfidence"=peakConf, 
				"DNA_Calls"=dna, "AA_Calls"=aa, "Filename"=chromoObj1$Filename)
	out
}


`modelChromatogram` <- function( obsChromo, seq=obsChromo$DNA_Calls[1], fixedPeaks=TRUE, 
				effort=1) {

	# given an observed chromatogram, fit a sequence to it

	# get a rough sense of the amplitudes
	obsTrace <- obsChromo$TraceM
	NP <- length( obsChromo$PeakPosition)
	guess.height <- quantile( as.vector(obsTrace), 0.95)
	max.height <- guess.height * 5
	min.height <- guess.height * 0.01
	guess.width <- 3
	max.width <- 4
	min.width <- 1
	guess.center <- 0
	max.center <- 2
	min.center <- -2

	if ( fixedPeaks) {
		start.height <- guess.height
		start.width <- guess.width
		start.center <- guess.center
		lowerBounds <- c( min.height, min.width, min.center)
		upperBounds <- c( max.height, max.width, max.center)
		max.iter <- 100
		max.time <- 100
		max.calls <- 1000000
	} else {
		start.height <- rep.int( guess.height, NP)
		start.width <- rep.int( guess.width, NP)
		start.center <- rep.int( guess.center, NP)
		lowerBounds <- c( rep.int(min.height,NP), rep.int(min.width,NP), rep.int(min.center,NP)) 
		upperBounds <- c( rep.int(max.height,NP), rep.int(max.width,NP), rep.int(max.center,NP))
		max.iter <- round( 25 * NP * effort)
		max.time <- round( 10 * NP * effort)
		max.calls <- round( 20000 * NP * effort)
	}

	# set up for NLS
	#controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	#starts <- list( "gaussian.height"=guess.height, "gaussian.width"=guess.width)
	#x <- seq
	#y <- as.vector( obsTrace)
	#nlsAns <- nls( y ~ synthChromatogram( x, gaussian.height, gaussian.width), 
	#			start=starts, control=controlList, algorithm="port", lower=lowerBounds,
	#			upper=upperBounds)
	#return( nlsAns)

	# set up for Generalize Simulated Annealing
	require( GenSA)

	# we can stop if we explain 95% of the observed data
	stopValue <- sqrt( sum( obsTrace ^ 2)) * 0.05

	# GenSA wants one vector of parameter weights
	wts <- c( start.height, start.width, start.center)
	control.list <- list( "maxit"=max.iter, "threshold.stop"=stopValue, "smooth"=FALSE, "max.call"=max.calls,
				"max.time"=max.time, "trace.mat"=FALSE)

	# make the call
	ans <- GenSA( par=wts, lower=lowerBounds, upper=upperBounds, fn=genSA.Chromatogram.residual,
			control=control.list, obs=obsTrace, seq=seq)

	# extract and apply the optimal parameters
	resid <- ans$value
	pars <- ans$par
	iters <- ans$counts
	if ( fixedPeaks) {
		ht <- pars[1]
		wid <- pars[2]
		ctr <- pars[3]
	} else {
		npar <- round(length(pars)/3)
		ht <- pars[1:npar]
		wid <- pars[(npar+1):(npar*2)]
		ctr <- pars[(npar*2+1):(npar*3)]
	}

	# make the optimal model
	chromoOut <- syntheticChromatogram( seq, height=ht, width=wid, center=ctr, traceOnly=FALSE)

	# append the modeling details
	chromoOut$FitPeakHeight <- ht
	chromoOut$FitPeakWidth <- wid
	chromoOut$FitPeakCenter <- ctr
	chromoOut$Residual <- resid
	chromoOut$Iterations <- iters

	return( chromoOut)
}

