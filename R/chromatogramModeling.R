# chromatogramModeling.R -- pieces to do modeling of chromatograms


`standardizeChromatogram` <- function( chromoObj, peak.dist=11, call.Ns=TRUE, constant.height=FALSE) {

	# turn the chromatogram into a uniform set of fixed width peaks, so we can add/model between multiple chromatograms
	
	# allow being given a filename
	if ( is.character(chromoObj) && file.exists( chromoObj[1])) {
		chromoObj <- loadChromatogram( chromoObj)
	}
	
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
		useHalf <- max( 1, half.width-1)
		obsHeight[i] <- sum( traceOut[ (centerOut-useHalf):(centerOut+useHalf), ])
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
		# some data has a long tail of super low peaks, that should not have a strong impact
		idealHeight <- max( median( obsHeight), mean(obsHeight), quantile(obsHeight, 0.75))
		minHeight <- idealHeight * 0.025
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
		
		# small chance that the scaling caused extreme large values
		# crop at something sensible
		max.value.in <- quantile( traceIn, 0.99)
		traceOut[ traceOut > max.value.in] <- max.value.in
	}

	# since we did interpolation on the trace matrix, clean it's significant digits some
	traceOut <- round( traceOut, digits=2)
	traceOut[ traceOut < 0] <- 0
	
	# lastly, recall the peak confidence
	peakConfidence <- calcChromatogramPeakConfidence( traceOut, peaksOut)

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=dna, "AA_Calls"=aa, "Filename"=chromoObj$Filename)
	out
}


`syntheticChromatogram` <- function( seq, peak.dist=11, height=1000, width=2.5, center=0, traceOnly=F, 
									addJitter=FALSE) {

	# create a synthetic chromatogram, given a DNA sequence and some constants
	peak.dist <- as.integer( peak.dist)
	half.width <- floor( peak.dist/2)
	NP <- nchar( seq)
	calls <- strsplit( toupper(seq), split="")[[1]]
	if ( ! all( calls %in% c("A","C","G","T","N","-"))) {
		cat( "\nError:  Given some non-ACGT base calls:\n")
		print( table( calls))
		return( NULL)
	}

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
		} else if (calls[i] != "-" ) {
			# allow gaps to be nothing added, otherwise assume N
			thisV <- thisV * 0.25
			for (k in 1:4) {
				traceM[thisLeft:thisRight,k] <- traceM[thisLeft:thisRight,k] + thisV
			}
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
	traceM[ traceM < 0] <- 0
	
	# some synthetic data may be too perfect for the model fitting tools, so allow the
	# choice to add a bit of noise
	if (addJitter) {
		traceM <- jitter( traceM, factor=1, amount=(height/1000))
		traceM[ traceM < 0] <- 0
	}

	# lastly, call the peak confidence
	peakConfidence <- calcChromatogramPeakConfidence( traceM, peaks)
	
	out <- list( "TraceM"=traceM, "PeakPosition"=peaks, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=dna, "AA_Calls"=aa, "Filename"="")
	out
}


`syntheticTraceMatrix` <- function( seq, peak.dist=11, height=1000, width=2.5, center=0) {

	# bares bones wrapper around synthetic chromatogram, given a DNA sequence and some constants
	# returns just the trace matrix
	ans <- syntheticChromatogram( seq, peak.dist=peak.dist, height=height, width=width, center=center, traceOnly=TRUE)
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


`modelFitChromatogram` <- function( obsChromo, seq=obsChromo$DNA_Calls[1], fixedPeaks=TRUE, 
				effort=1, doStandardize=TRUE, doSubset=TRUE, algorithm=c("nls","GenSA")) {

	# given an observed chromatogram, fit a sequence to it and return a best-fit model
	# chromatogram, suitable for subtraction to generate a residual chromatogram

	# allow being given a filename
	if ( is.character(obsChromo) && file.exists( obsChromo[1])) {
		obsChromo <- loadChromatogram( obsChromo)
	}
	
	# the modeling needs a standardized starting object that matches the sequence roughly
	if (doStandardize) {
		obsChromo <- standardizeChromatogram( obsChromo)
	}
	if (doSubset) {
		obsChromo <- subsetChromatogram( obsChromo, seq=seq)
	}
	
	# get a rough sense of the amplitudes
	obsTrace <- obsChromo$TraceM
	NP <- length( obsChromo$PeakPosition)
	guess.height <- as.numeric( quantile( as.vector(obsTrace), 0.95))
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

	algorithm <- match.arg( algorithm)
	# set up for NLS
	if (algorithm == "nls" && fixedPeaks) {
		controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
		starts <- list( "height"=guess.height, "width"=guess.width, "center"=guess.center)
		cat( "\nDebug: \n")
		print( starts)
		print( lowerBounds)
		print( upperBounds)
		x <- seq
		y <- obsTrace
		nlsAns <<- nls( y ~ syntheticTraceMatrix( x, peak.dist=11, height, width, center), 
					start=starts, control=controlList, algorithm="port", lower=lowerBounds,
					upper=upperBounds)
		nlsAns2 <<- summary( nlsAns)

	} else {
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


`modelBlendChromatogram` <- function( obsChromo, seqs, plot.chromatograms=T, 
									synthetic.width=2.5, min.pct.plot=5, max.pval.plot=0.05,
									label="") {

	# given an observed chromatogram, fit 2 or more sequences to it and return the contribution of
	# how much of each sequence was needed to best fit the observed chromatogram.
	
	# This is akin to deconvolution, but component sequences must be given, not deduced automatically

	# allow being given a filename of a chromatogram
	if ( is.character(obsChromo) && file.exists( obsChromo[1])) {
		obsChromo <- loadChromatogram( obsChromo)
	}
	
	# default label is the file's name
	if ( label == "" && "Filename" %in% names(obsChromo)) {
		label <- sub( ".ab1$", "", basename( obsChromo$Filename))
	}
	
	# there is a bunch of setup steps, to assure the sequences and chromatogram can be directly compared
	setupOK <- TRUE
	# 1.  All the sequences must be exact same size.  Gap characters are allowed.
	if ( length(seqs) < 2) {
		cat( "\nMust give 2+ DNA sequences for 'modelBlend'")
		setupOK <- FALSE
	}
	if ( length( unique( nchar(seqs))) > 1) {
		# catch/allow gaps to be present without braking this rule
		seqs2 <- gsub( "-", "", seqs, fixed=T)
		if ( length( unique( nchar(seqs2))) > 1) {
			cat( "\nAll DNA sequences for 'modelBlend' must be same length.  Explicit gap characters are allowed.")
			setupOK <- FALSE
		} else {
			# removing the gaps put all at the same length, so use those version of the seqs
			seqs <- seqs2
		}
	}
	if ( is.null( names(seqs))) {
		cat( "\nDNA sequences for 'modelBlend' must have 'names' attributes")
		setupOK <- FALSE
	}
	if ( ! setupOK) return( NULL)
	NS <- length( seqs)
	seqs <- toupper( seqs)
	noGapSeqs <- gsub( "-", "", seqs, fixed=T)
	
	# 2.  We need to find the best matching sequence to know which to use to crop the chromatogram
	scores <- numeric(NS)
	stdChromo <- standardizeChromatogram( obsChromo, constant.height=TRUE)
	for ( i in 1:NS) {
		tmpChromo <- subsetChromatogram( stdChromo, seq=noGapSeqs[i])
		scores[i] <- tmpChromo$UnitScore
	}
	best <- which.max( scores)
	bestSeq <- seqs[ best]
	# debug:  look at the scores
	names(scores) <- names(seqs)
	tmpChromo <- subsetChromatogram( stdChromo, seq=bestSeq)	
	
	# 3. We need to know if we are doing Fwd or RevComp to get the sequences facing the right way
	editDist <- adist( bestSeq, tmpChromo$DNA_Calls)
	if ( which.min( editDist) == 2) {
		for (i in 1:NS) seqs[i] <- myReverseComplement( seqs[i])
		bestSeq <- seqs[ best]
	}
	
	# with strand known, we now need to move any gap bases to the end, to correctly represent what
	# that sequence's chromatogram would actually look like
	doGaps <- grep( "-", seqs, fixed=T)
	if ( length( doGaps)) {
		for (i in doGaps) {
			bp <- strsplit( seqs[i], split="")[[1]]
			isgap <- which( bp == "-")
			chunk1 <- bp[ -isgap]
			chunk2 <- bp[ isgap]
			seqs[i] <- paste( c(chunk1, chunk2), collapse="")
		}
		bestSeq <- seqs[ best]
	}
	
	# with the best starting sequence chosen, and Fwd/Rev known, and gaps moved to the ends, 
	# now we can extract that portion from the observed chromatogram,
	observedChromo <- subsetChromatogram( stdChromo, seq=bestSeq)	
	maxObsHeight <- quantile( observedChromo$TraceM, 0.99)
	observedTraceM <- observedChromo$TraceM
	
	# let's allow using the model fit algorithm to optimize the peak width we use
	if ( is.character( synthetic.width) && synthetic.width == "fit") {
		cat( "\nFitting best model element sequence to observed data..")
		fitAns <- modelFitChromatogram( observedChromo, seq=bestSeq, doStandardize=F, doSubset=F, algorithm="GenSA")
		synthetic.width <- as.numeric( fitAns$FitPeakWidth)
		cat( "    Optimal Peak Width =", synthetic.width)
	}
	
	# finally ready to create synthetic chromatograms for each sequence to be modeled
	synthChromoTraceMs <- synthChromos <- vector( mode="list")
	synthSizeError <- FALSE
	for ( i in 1:NS) {
		synthChromos[[i]] <- tmpChromo <- syntheticChromatogram( seqs[i], height=maxObsHeight, width=synthetic.width)
		synthChromoTraceMs[[i]] <- tmpChromo$TraceM
		# small chance that we have a problem making all the chromatograms be the exact same size
		if ( nrow(tmpChromo$TraceM) != nrow(observedTraceM)) synthSizeError <- TRUE
	}
	# if any flagged for bad size, die now
	if (synthSizeError) {
		cat( "\nError creating model chromatogram elements.  Unable to do model fit..")
		return( NULL)
	}
	
	# add in some noise sequences to be more realistic
	noiseSeqs <- vector()
	noiseBases <- c("A","C","G","T","N")
	NN <- length( noiseBases)
	noiseTraceMs <- noiseChromos <- vector( mode="list")
	for (i in 1:NN) {
		ch <- noiseBases[i]
		noiseStr <- paste( rep.int(ch,nchar(bestSeq)), collapse="")
		noiseSeqs[i] <-  noiseStr
		noiseChromos[[i]] <- tmpChromo <- syntheticChromatogram( noiseStr, height=maxObsHeight, 
														width=synthetic.width, addJitter=T)
		noiseTraceMs[[i]] <- tmpChromo$TraceM
	}
	names(noiseSeqs) <- paste( "Poly", noiseBases, sep="_")
	
	# join our real choices and the noise choices
	seqs <- c( seqs, noiseSeqs)
	NS <- length( seqs)
	synthChromos <- c( synthChromos, noiseChromos)
	synthChromoTraceMs <- c( synthChromoTraceMs, noiseTraceMs)
	
	# ready to do the NLS fitting...   set up the controls
	controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	starts <- list( "weights"=jitter( rep.int( 1/NS, NS)))
	lowerBounds <- rep.int( 0, NS)
	upperBounds <- rep.int( 10, NS)
	x <- synthChromoTraceMs
	y <- as.vector( observedChromo$TraceM)
	
	# define the called NLS blending function
	blendChromatograms <- function( x, weights) {
				N <- length( x)
				v <- as.vector( x[[1]]) * weights[1]
				if ( N > 1) for ( j in 2:N) {
					v2 <- as.vector( x[[j]]) * weights[j]
					v <- v + v2
				}
				return( v)
	}
	
	# call the NLS to do the fit
	nlsAns <- nls( y ~ blendChromatograms( x, weights), 
				start=starts, control=controlList, algorithm="port", lower=lowerBounds,
				upper=upperBounds)
	
	# extract the results
	nlsAns2 <- summary( nlsAns)
	coefs <- coefficients( nlsAns2)
	blendEsts <- rawEsts <- coefs[ ,1]
	names( blendEsts) <- names( seqs)
	pvals <- coefs[,4]
	pvalText <- ifelse( pvals < 0.001, formatC( pvals, format="e", digits=2), as.character( round(pvals, digits=3)))
	
	# turn these estimates into percentages
	pcts <- round( blendEsts * 100 / sum( blendEsts), digits=3)
	out1 <- data.frame( "Construct"=names(seqs), "Percentage"=pcts, "P.value"=pvals, stringsAsFactors=F)
	
	# order them, and the blend estimates too
	seqOrd <- order( out1$Percentage, decreasing=T)
	seqOrd <- order( out1$P.value, decreasing=F)
	out1 <- out1[ seqOrd, ]
	blendEsts <- blendEsts[ seqOrd]
	rownames(out1) <- 1:nrow(out1)
	
	# quantify how well the model fit the observed data
	obsInten <- sum( y)
	totResid <- sum( abs( nlsAns2$residuals))
	resPct <- round( totResid * 100 / obsInten, digits=2)
	modelPct <- 100 - resPct
	
	out <- list( "Model.Estimates"=out1, "Model.Fit.Percentage"=modelPct)
	
	# do we want to draw what we did?
	if (plot.chromatograms) {
		# create the residual by subtracting all non-zero components
		residChromo <- observedChromo
		residM <- residChromo$TraceM
		toDraw <- which( pcts >= min.pct.plot | pvals < max.pval.plot)
		NtoDraw <- length(toDraw)
		cex <- 1 - (0.05*NtoDraw)
		savMAI <- par( "mai")
		par( mfrow=c(NtoDraw+2, 1))
		par( mai=c(0.5,0.9,0.3,0.4))
		plotChromatogram( observedChromo, forceYmax=maxObsHeight, label=paste( "Observed Data:  ", label), cex=cex, cex.main=1.5)
		for ( j in seqOrd) {
			myEst <- rawEsts[j]
			myChromo <- synthChromos[[j]]
			myM <- myChromo$TraceM * myEst
			myChromo$TraceM <- myM
			residM <- residM - myM
			residM[ residM < 0] <- 0
			if( j %in% toDraw) {
				myYheight <- maxObsHeight * sqrt(myEst)
				plotChromatogram( myChromo, forceYmax=myYheight, cex=cex, cex.main=1.5,
							label=paste( "Model Element:  ", names(seqs)[j], " = ", round(pcts[j],digits=1), "%    P.value = ", pvalText[j], sep=""))
				dev.flush()
			}
		}
		# lastly, show the final residual
		residChromo$TraceM <- residM
		plotChromatogram( residChromo, forceYmax=maxObsHeight, cex=cex, cex.main=1.5, 
					label=paste( "Model Residual   (model explains ", modelPct, "% of raw intensity)", sep=""))	
		dev.flush()
		par( mfrow=c(1,1))
		par( mai=c(1,1,0.8,0.4))

	}
	
	return( out)
}
	


`deconvoluteChromatogram` <- function( obsChromo, seq=NULL, range=NULL, plot.chromatograms=T, doBlast=TRUE,
									blastTarget=c("protein","nucleotide"),
									label="", path=NULL, max.constructs=3, max.plot=2, verbose=T) {

	# given an observed chromatogram, try to automatically and iteratively deduce what constructs it
	# is composed of, by repeatedly modelling and subtracting the dominant signal
	
	# This is akin to model blending, but no sequences need be given a priori

	# allow being given a filename of a chromatogram
	if ( is.character(obsChromo) && file.exists( obsChromo[1])) {
		obsChromo <- loadChromatogram( obsChromo)
	}
	
	# decide where/what to call the files we make
	chromoFile <- obsChromo$Filename
	if ( is.null(path)) path <- dirname( chromoFile)
	baseFile <- sub( ".ab1$", "", basename(chromoFile))
	
	# default label is the file's name
	if ( label == "" && "Filename" %in% names(obsChromo)) {
		label <- baseFile
	}

	# were we given a request for a smaller region?
	if ( ! is.null( seq)) {
		obsChromo <- subsetChromatogramBySequence(obsChromo, seq=seq)
		# allow the subsequence request to return silently with no plot if not good sequence found
		if ( is.null( obsChromo)) return( NULL)
	} else if ( ! is.null( range)) {
		obsChromo <- subsetChromatogramByRange(obsChromo, range=range)
	}
	
	# accumulate the answers we will return
	dnaOut <- vector()
	intenOut <- vector()
	nameOut <- vector()
	nIter <- 0
	totalModelInten <- 0
	chromoToPlot <- vector( mode="list")
	residToPlot <- vector( mode="list")
	
	# get the starting trace matrix, to know when we should stop iterating
	# but first crop off any decayed tail and standarize the peak spacing
	obsChromo <- cropChromatogramLowSignalTail( obsChromo)
	obsChromo <- standardizeChromatogram( obsChromo, constant.height=TRUE)
	obsTraceM <- obsChromo$TraceM
	obsInten <- sum( obsTraceM, na.rm=T)
	stopInten <- obsInten * 0.05

	# repeat until we stall out or succeed
	curChromo <- obsChromo
	repeat {
	
		# watch for any reason to bail out
		if (nIter >= max.constructs) break
		curInten <- sum( curChromo$TraceM, na.rm=T)
		if ( curInten <= stopInten) break
	
		# take the DNA sequence as currently defined
		curDNA <- curChromo$DNA_Calls[1]
		
		# perhaps trim to exactly fit the sequence model will return
		curChromo <- subsetChromatogram( curChromo, seq=curDNA)
		
		# make the best fit model for this
		cat( "\nModel..")
		modelChromo <- modelFitChromatogram( curChromo, seq=curDNA, fixedPeaks=TRUE, effort=1, 
									doStandardize=TRUE, doSubset=FALSE, algorithm="GenSA")
		modelInten <- sum( modelChromo$TraceM, na.rm=T)
		
		# accumulate these answers
		nIter <- nIter + 1
		dnaOut[ nIter] <- curDNA
		intenOut[ nIter] <- modelInten
		nameOut[nIter] <- paste( "Construct", nIter, sep="_")
		cat( "  Answer=  ", nIter, "  Intensity= ", modelInten, "  Percent = ", round( modelInten * 100 / obsInten, digits=1))
		if (plot.chromatograms) chromoToPlot[[nIter]] <- modelChromo
		
		# subtract this model from the current chromatogram
		cat( "  Subtract..")
		deltaChromo <- subtractChromatogram( curChromo, modelChromo)
		if ( is.null( deltaChromo)) break
		if (plot.chromatograms) residToPlot[[nIter]] <- deltaChromo
		
		# re-pick where the peaks are now, watching for any errors
		cat( "  Re-PeakPick..")
		newChromo <- peakpickChromatogram( deltaChromo)
		if ( is.null( newChromo)) break
		# and restandardize it
		newChromo <- standardizeChromatogram( newChromo, constant.height=TRUE)
		
		# there is a chance that the chromatogram has more peaks than originally, due to how the peak picker operates
		# and the nature of the raw noisy data.  Try to prevent peak count creep.
		firstModel.size <- nchar( dnaOut[1])
		this.size <- nchar( newChromo$DNA_Calls[1])
		if ( this.size > firstModel.size) {
			croppedDNA <- substr( newChromo$DNA_Calls[1], 1, firstModel.size)
			newChromo <- subsetChromatogramBySequence( newChromo, croppedDNA)
		}
			
		# use this new shorter chromatogram with its own new peak locations/heights/calls, as the new current
		curChromo <- newChromo
	}

	# when we fall out, summarize and return the results
	residInten <- sum( curChromo$TraceM, na.rm=T)	
	intenPcts <- round( intenOut * 100 / obsInten, digits=1)
	out <- data.frame( "Name"=nameOut, "Percentage"=intenPcts, "DNA_Sequence"=dnaOut, stringsAsFactors=F)
	
	# submit the constructs to BLAST?
	fastaFile <- file.path( path, paste( baseFile, "Deconvolution.Constructs.fasta", sep="."))
	blastFile <- file.path( path, paste( baseFile, "Deconvolution.BlastOutput.xml", sep="."))
	resultsFile <- file.path( path, paste( baseFile, "Deconvolution.Results.csv", sep="."))
	plotFile <- file.path( path, paste( baseFile, "Deconvolution.Results.pdf", sep="."))
	if ( doBlast || !file.exists(blastFile) || (file.exists(blastFile) && file.info(blastFile)$size < 1000)) {
		writeFasta( as.Fasta( nameOut[1:nIter], dnaOut[1:nIter]), fastaFile)
		blastTarget <- match.arg( blastTarget)
		if ( blastTarget == "protein") {
			callBlastx( fastaFile, outfile=blastFile, evalue=10, wordsize=5, threads=6, outfmt=5, maxhits=1, verbose=verbose)
		} else {
			callBlastn( fastaFile, outfile=blastFile, evalue=10, wordsize=9, threads=6, outfmt=5, maxhits=1, verbose=verbose)
		}
	}
	cat( "\nExtracting results from XML..")
	blastAns <- extractBlastXMLdetails( blastFile, IDs=nameOut[1:nIter], IDprefix="")
	# supplement the output with what we found
	acc <- unlist( sapply( blastAns, `[[`, "accession"))
	def <- unlist( sapply( blastAns, `[[`, "definition"))
	eval <- unlist( sapply( blastAns, `[[`, "evalue"))
	scor <- unlist( sapply( blastAns, `[[`, "score"))
	matchStr <- unlist( sapply( blastAns, `[[`, "match"))
	# small chance that not all constructs gave valid hits...
	if ( length(acc) < nIter) acc <- c( acc, rep.int( "NA", nIter-length(acc)))
	if ( length(def) < nIter) def <- c( def, rep.int( "No good scoring hits found", nIter-length(def)))
	if ( length(eval) < nIter) eval <- c( eval, rep.int( 1, nIter-length(eval)))
	if ( length(scor) < nIter) scor <- c( scor, rep.int( 0, nIter-length(scor)))
	if ( length(matchStr) < nIter) matchStr <- c( matchStr, rep.int( "", nIter-length(matchStr)))
	# and the residual gets no details
	out$Accession <- acc
	out$Definition <- gsub( "&gt;", "  ", def, fixed=T)
	out$Evalue <- eval
	out$Score <- round( as.numeric(scor), digits=2)
	out$MatchString <- matchStr
	write.table( out, resultsFile, sep=",", quote=T, row.names=F)

	# do we want to draw what we did?
	if (plot.chromatograms) {
		if ( length(chromoToPlot) > max.plot) length(chromoToPlot) <- max.plot
		if ( length(residToPlot) > max.plot) length(residToPlot) <- max.plot
		NtoDraw <- 1 + length(chromoToPlot) + length(residToPlot)
		cex <- 1 - (0.05*NtoDraw)
		savMAI <- par( "mai")
		par( mfrow=c(NtoDraw, 1))
		par( mai=c(0.5,0.9,0.3,0.4))
		maxObsHeight <- quantile( apply( obsTraceM,1,max,na.rm=T),0.99)
		plotChromatogram( obsChromo, forceYmax=maxObsHeight, label=paste( "Observed Data:  ", label), cex=cex, cex.main=1.8)
		for ( j in 1:length(chromoToPlot)) {
			myChromo <- chromoToPlot[[j]]
			plotChromatogram( myChromo, forceYmax=maxObsHeight, cex=cex, cex.main=1.8,
							label=paste( "Deconvolution ", out$Name[j], " = ", intenPcts[j], "%    Score = ", out$Score[j], 
							"   Evalue = ", out$Evalue[j], sep=""))
			if ( "Definition" %in% colnames(out)) {
				textToShow <- out$Definition[j]
				if ( ! is.na(textToShow)) {
					textToShow <- gsub( "&gt;", "  ", textToShow, fixed=T)
					if ( nchar(textToShow) > 200) textToShow <- paste( substr( textToShow, 1, 200), "...", sep="")
					text( nrow(myChromo$TraceM)/2, maxObsHeight*0.925,  textToShow, cex=2)
				}
			}
			if ( "MatchString" %in% colnames(out)) {
				if ( ! is.na(out$Score[j]) && out$Score[j] > 10) {
					require( plotrix)
					textToShow <- out$MatchString[j]
						if ( ! is.na(textToShow) && nchar(textToShow) > 10) {
						textToShow <- strsplit( textToShow, split="\n")[[1]]
						textToShow <- sub( "^    ", "", textToShow)
						textToShow <- paste( c( "Construct  ", "           ", "Blast Hit  "), textToShow, sep="")
						textCharWidth <- nrow(myChromo$TraceM) / 400
						textX <- (nrow(myChromo$TraceM)/2) + (nchar(textToShow[1]) * textCharWidth * c(-1,1))
						textY <- maxObsHeight * c(0.15, 0.51)
						rect( textX[1], textY[1], textX[2], textY[2], border='black', col='white')
						text( nrow(myChromo$TraceM)/2, maxObsHeight*.43, textToShow[1], cex=1, family="mono")
						text( nrow(myChromo$TraceM)/2, maxObsHeight*.33, textToShow[2], cex=1, family="mono")
						text( nrow(myChromo$TraceM)/2, maxObsHeight*.23, textToShow[3], cex=1, family="mono")
					}
				}
			}
			dev.flush()
			myChromo <- residToPlot[[j]]
			myIntenPct <- round( sum( myChromo$TraceM, na.rm=T) * 100 / obsInten, digits=1)
			plotChromatogram( myChromo, forceYmax=maxObsHeight, cex=cex, cex.main=1.8,
							label=paste( "Residual after ", out$Name[j], " = ", myIntenPct, "%", sep=""))
			dev.flush()
		}
		dev.print( pdf, plotFile, width=20)
		par( mfrow=c(1,1))
		par( mai=c(1,1,0.8,0.4))
	}
	return( out)
}
