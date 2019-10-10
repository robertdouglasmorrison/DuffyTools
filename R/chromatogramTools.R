# chromatogramTools.R -- pieces to work directly from ABI chromatogram files
#			and to manipulate the underlaying A/C/G/T spectra


`loadChromatogram` <- function( chromoFile) {

	# read up an AB1 file from disk... into our standard format of trace, peaks, & sequence
	require( "sangerseqR")

	if ( ! file.exists( chromoFile)) {
		cat( "\nChromatogram file not found:  ", chromoFile)
		return( NULL)
	}
	ab1 <- sangerseq( read.abif( chromoFile))

	# the raw data is 3 parts:   a DNA sequence,  a 4x matrix of base intensities,  
	# and a matrix of peak locations
	primary.seq <- as.character( primarySeq(ab1))
	traceM <- traceMatrix( ab1)
	peakPosM <- peakPosMatrix( ab1)

	# we will always keep both strands, as both DNA and AA
	DNA.primary.set <- c( primary.seq, myReverseComplement(primary.seq))
	names(DNA.primary.set) <- c( "DNA", "rc(DNA)")
	AA.primary.set <- DNAtoAA( primary.seq, clipAtStop=F, readingFrame=1:6)
	names(AA.primary.set) <- paste( "AA_Frame", 1:6, sep="")

	# for now, just take the primary peak positions as ABI called them
	peakPos <- as.numeric( peakPosM[ ,1])
	names(peakPos) <- strsplit( primary.seq, split="")[[1]]

	# use those peak top locations to estimate where the trace stops being useful
	avgPeakSeparation <- median( diff( peakPos))
	rightTail <- max( peakPos) + avgPeakSeparation
	if ( rightTail < nrow(traceM)) {
		traceM <- traceM[ 1:rightTail, ]
	}
	colnames(traceM) <- c( "A", "C", "G", "T")
	rownames(traceM) <- 1:nrow(traceM)

	out <- list( "TraceM"=traceM, "PeakPosition"=peakPos, "DNA_Calls"=DNA.primary.set,
			"AA_Calls"=AA.primary.set)
	out
}


`revCompChromatogram` <- function( chromoObj) {

	traceM <- revCompTraceMatrix( chromoObj$TraceM)
	tmpPeakPos <- chromoObj$PeakPosition
	# turn into a matrix for the RevComp, then take back just column 1
	tmpPeakPosM <- matrix( tmpPeakPos, nrow=length(tmpPeakPos), ncol=1)
	#cat( "\nDebug:  typeof:  ", typeof( tmpPeakPos), "  dim: ", dim( tmpPeakPos))
	peakPosM <- revCompPeakPosMatrix( tmpPeakPosM, nrow.TM=nrow(traceM))
	peakPos <- peakPosM[ ,1]

	dna <- rev( chromoObj$DNA_Calls)
	aa <- DNAtoAA( dna[1], clipAtStop=F, readingFrame=1:6)
	
	names(peakPos) <- strsplit( dna[1], split="")[[1]]

	out <- list( "TraceM"=traceM, "PeakPosition"=peakPos, "DNA_Calls"=dna,
			"AA_Calls"=aa)
	out
}


`revCompTraceMatrix` <- function( traceM) {

	N <- nrow(traceM)
	# first do the reverse
	out <- tmp <- traceM[ N:1, ]
	# next do the  complement (A/C/G/T)
	out[ , 1] <- tmp[ ,4]
	out[ , 2] <- tmp[ ,3]
	out[ , 3] <- tmp[ ,2]
	out[ , 4] <- tmp[ ,1]
	out
}


`revCompPeakPosMatrix` <- function( peakPosM, nrow.TM) {

	N <- nrow(peakPosM)
	# first do the reverse
	out <- peakPosM[ N:1, ,drop=FALSE]
	# next update all the location pointers given the trace matrix extent
	for ( i in 1:ncol( peakPosM)) {
		out[ ,i] <- (nrow.TM - out[ ,i] + 1)
	}
	out
}


`baseCallOnePeak` <- function( peakSite, traceM, min.pct=0.10) {

	# use the data in a tiny window, not just the 1 point
	intenV <- apply( traceM[ (peakSite-1):(peakSite+1), ], MARGIN=2, sum)
	intenPcts <- intenV / sum(intenV)
	# good to call if at least 15% of total intensity
	good <- which( intenPcts >= min.pct)
	outBase <- colnames(traceM)[ good]
	if ( ! length(outBase)) {
		outBase <- "N"
		outPcts <- 0
		outOrd <- 1
	} else {
		outPcts <- round( intenPcts[ good], digits=2)
		outOrd <- order( outPcts, decreasing=T)
	}
	out <- outPcts[outOrd]
	names(out) <- outBase
	return( out)
}


`baseCallProfile` <- function( chromoObj, min.pct=0.10) {

	traceM <- chromoObj$TraceM
	peaks <- chromoObj$PeakPosition
	NP <- length( peaks)

	out <- lapply( peaks, baseCallOnePeak, traceM=traceM, min.pct=min.pct)
	names(out) <- paste( 1:NP, names(peaks), sep="_")

	return( out)
}


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
			dna <- c( newSeq, myReverseComplement(newSeq))
			names(dna) <- c( "DNA", "rc(DNA)")
			aa <- DNAtoAA( newSeq, clipAtStop=F, readingFrame=1:6)
			names(aa) <- paste( "AA_Frame", 1:6, sep="")
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

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "DNA_Calls"=dna, "AA_Calls"=aa)
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
		dna <- c( seq, myReverseComplement(seq))
		names(dna) <- c( "DNA", "rc(DNA)")
		aa <- DNAtoAA( seq, clipAtStop=F, readingFrame=1:6)
		names(aa) <- paste( "AA_Frame", 1:6, sep="")
	}

	# since we did interpolation on the trace matrix, clean it's significant digits some
	traceM <- round( traceM, digits=2)

	out <- list( "TraceM"=traceM, "PeakPosition"=peaks, "DNA_Calls"=dna, "AA_Calls"=aa)
	out
}


`synthChromatogram` <- function( x, height=1000, width=2.5, center=0) {

	# bares bones wrapper around synthetic chromatogram, given a DNA sequence and some constants
	#cat( "\nDebug Synth:  ", gaussian.height, gaussian.width)
	ans <- syntheticChromatogram( x, height=height, width=width, center=center, traceOnly=TRUE)
	return( ans$TraceM)
}


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

	ansModel <- synthChromatogram( seq, height=ht, width=wid, center=ctr)
	diff <- ( obs - ansModel)
	resid <- sqrt( sum( diff * diff))
	return( resid)
}


`subsetChromatogram` <- function( chromoObj, substring=NULL, range=NULL) {

	# get the data we need 
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	baseCall <- names( peakPos)
	NT <- nrow( traceM)
	NB <- length( peakPos)
	halfPeak <- floor( median( diff( peakPos)) / 2)

	# allow using just a subset of the full length
	firstTracePoint <- 1
	lastTracePoint <- NT

	if ( is.null(substring) && is.null(range)) stop( "Must specifiy a subset by substring or range..")

	if ( ! is.null( substring)) {
		# we will find DNA or AA in any reading frame
		subSeq <- as.character( substring)
		chromoDNA <- chromoObj$DNA_Calls
		chromoAA <- chromoObj$AA_Calls

		# look in both DNA and AA, (since it's 2 different scoring matrices, compensate a bit
		require( Biostrings)
		dnaScores <- pairwiseAlignment( chromoDNA, subSeq, type="local-global", scoreOnly=T) * 3
		aaScores <- pairwiseAlignment( chromoAA, subSeq, type="local-global", scoreOnly=T)
	
		if ( max( dnaScores) > max(aaScores)) {
			best <- which.max( dnaScores)
			bestDNA <- chromoDNA[ which.max( dnaScores)]
			pa <- pairwiseAlignment( bestDNA, subSeq, type="local-global", scoreOnly=F)
			startDNA <- start( pattern( pa))
			stopDNA <- startDNA + width( pattern( pa)) - 1
			needChromoRevComp <- (best > 1)
			AAoffset <- 1
		} else {
			best <- which.max( aaScores)
			bestAA <- chromoAA[ which.max( aaScores)]
			pa <- pairwiseAlignment( bestAA, subSeq, type="local-global", scoreOnly=F)
			startAA <- start( pattern( pa))
			stopAA <- startAA + width( pattern( pa)) - 1
			stopDNA <- (stopAA*3)
			startDNA <- startAA*3 - 2
			needChromoRevComp <- (best > 3)
			AAoffset <- if (best < 4) best else best - 3
		}

		# did we need to do a RevComp?
		if (needChromoRevComp) {
			tmpObj <- revCompChromatogram( chromoObj)
			traceM <- tmpObj$TraceM
			peakPos <- tmpObj$PeakPosition
			baseCall <- names( peakPos)
			NT <- nrow( traceM)
			NB <- length( peakPos)
		}

		# now turn those DNA spots into traceM locations, accounting for any reading frame shift
		startDNA <- startDNA + AAoffset - 1
		stopDNA <- stopDNA + AAoffset - 1
		firstTracePoint <- max( 1, peakPos[ startDNA] - halfPeak)
		lastTracePoint <- min( NT, peakPos[ stopDNA] + halfPeak)
	}

	if ( ! is.null(range)) {
		firstPeak <- peakPos[ min(range)]
		lastPeak <- peakPos[ max(range)]
		firstTracePoint <- max( 1, firstPeak - halfPeak)
		lastTracePoint <- min( NT, lastPeak + halfPeak)
	}

	# now trim all data to that subset
	xLimits <- c( firstTracePoint, lastTracePoint)
	traceOut <- traceM[ firstTracePoint:lastTracePoint, ]
	peaksOut <- peakPos[ peakPos > firstTracePoint & peakPos < lastTracePoint]
	peaksOut <- peaksOut - firstTracePoint + 1

	# make the other fields the chromatogram has
	seqOut <- paste( names(peaksOut), collapse="")
	dna <- c( seqOut, myReverseComplement(seqOut))
	names(dna) <- c( "DNA", "rc(DNA)")
	aa <- DNAtoAA( seqOut, clipAtStop=F, readingFrame=1:6)
	names(aa) <- paste( "AA_Frame", 1:6, sep="")

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "DNA_Calls"=dna, "AA_Calls"=aa)
	out
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
	peaksOut <- chromoObj1$PeakPosition
	dna <- chromoObj1$DNA_Calls
	aa <- chromoObj1$AA_Calls

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "DNA_Calls"=dna, "AA_Calls"=aa)
	out
}


`plotChromatogram` <- function( chromoObj, label="", subset=NULL, lwd=2, lty=1, cex=1, font=2, 
				add=FALSE, forceYmax=NULL, showAA=TRUE) {

	# given a chromatogram object, show the whole thing
	acgtBases <- c('A','C','G','T','N','-')
	acgtColors <- c('red','blue','orange','green','brown','black')

	neededNames <- c( "TraceM", "PeakPosition")
	if ( ! all( neededNames %in% names(chromoObj))) {
		cat( "Chromatogram object does not have needed fields..")
		return()
	}

	# get the data we need 
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	baseCall <- names( peakPos)
	NT <- nrow( traceM)
	NB <- length( peakPos)

	# allow using just a subset of the full length
	firstTracePoint <- 1
	lastTracePoint <- NT
	xLimits <- c( 1, NT)
	DO_SUBSET <- FALSE
	AAtoShow <- ""
	AAoffset <- 1

	# we will find DNA or AA in any reading frame
	if ( ! is.null( subset)) {
		subSeq <- as.character( subset)
		chromoDNA <- chromoObj$DNA_Calls
		chromoAA <- chromoObj$AA_Calls
		# look in both DNA and AA, (since it's 2 different scoring matrices, compensate a bit
		dnaScores <- pairwiseAlignment( chromoDNA, subSeq, type="local-global", scoreOnly=T) * 3
		aaScores <- pairwiseAlignment( chromoAA, subSeq, type="local-global", scoreOnly=T)

		if ( max( dnaScores) > max(aaScores)) {
			best <- which.max( dnaScores)
			bestDNA <- chromoDNA[ which.max( dnaScores)]
			pa <- pairwiseAlignment( bestDNA, subSeq, type="local-global", scoreOnly=F)
			startDNA <- start( pattern( pa))
			stopDNA <- startDNA + width( pattern( pa)) - 1
			needChromoRevComp <- (best > 1)
		} else {
			best <- which.max( aaScores)
			bestAA <- chromoAA[ which.max( aaScores)]
			pa <- pairwiseAlignment( bestAA, subSeq, type="local-global", scoreOnly=F)
			startAA <- start( pattern( pa))
			stopAA <- startAA + width( pattern( pa)) - 1
			stopDNA <- (stopAA*3)
			startDNA <- startAA*3 - 2
			needChromoRevComp <- (best > 3)
			AAtoShow <- chromoAA[ best]
			AAoffset <- if (best < 4) best else best - 3
			#cat( "\nDebug: ", best, needChromoRevComp, "|", startAA, stopAA, startDNA, stopDNA)
		}
		# did we need to do a RevComp?
		if (needChromoRevComp) {
			tmpObj <- revCompChromatogram( chromoObj)
			traceM <- tmpObj$TraceM
			peakPos <- tmpObj$PeakPosition
			baseCall <- names( peakPos)
			NT <- nrow( traceM)
			NB <- length( peakPos)
			#newStartDNA <- NB - stopDNA + 1
			#newStopDNA <- NB - startDNA + 1
			#startDNA <- newStartDNA
			#stopDNA <- newStopDNA
			if ( AAtoShow != "") {
				AAtoShow <- tmpObj$AA_Calls[AAoffset]
			}
		}

		# now turn those DNA spots into traceM locations, accounting for any reading frame shift
		startDNA <- startDNA + AAoffset - 1
		stopDNA <- stopDNA + AAoffset - 1
		halfPeak <- round( median( diff( peakPos)) / 2)
		firstTracePoint <- max( 1, peakPos[ startDNA] - halfPeak)
		lastTracePoint <- min( NT, peakPos[ stopDNA] + halfPeak)
		# and trim it to that subset
		xLimits <- c( firstTracePoint, lastTracePoint)
		DO_SUBSET <- TRUE
	} else {
		if ( is.logical(showAA) && showAA) AAtoShow <- chromoObj$AA_Calls[1]
		if ( is.numeric(showAA)) AAtoShow <- chromoObj$AA_Calls[showAA]
	}

	mainText <- paste( "Chromatogram:  ", label)
	x <- 1 : NT
	yLimits <- c( 0, max( traceM[ firstTracePoint:lastTracePoint, ], na.rm=T) * 1.1)
	if ( ! is.null(forceYmax)) yLimits[2] <- as.numeric( forceYmax[1])

	if ( ! add) plot( 1,1, type="n", main=mainText, xlim=xLimits, ylim=yLimits, ylab="Intensity", xlab=NA,
				xaxt="n", xaxs="i", las=2)

	for ( j in 1:4) {
		y <- traceM[ ,j]
		lines( x, y, col=acgtColors[j], lwd=lwd, lty=lty)
	}
	baseColor <- acgtColors[ match( baseCall, acgtBases)]
	# flag indel cleaning with extra color
	baseColor[ is.na(baseColor)] <- 'purple'

	if ( ! DO_SUBSET) {
		for (k in 1:NB) axis( side=1, at=peakPos[k], label=baseCall[k], col.axis=baseColor[k], 
				col.ticks=baseColor[k], font=font, cex.axis=cex)
	} else {
		for (k in 1:NB) {
			thisX <- peakPos[k]
			if (thisX < firstTracePoint || thisX > lastTracePoint) next
			axis( side=1, at=peakPos[k], label=baseCall[k], col.axis=baseColor[k], 
				col.ticks=baseColor[k], font=font, cex.axis=cex)
		}
	}

	if ( AAtoShow != "") {
		aaCall <- strsplit( AAtoShow, split="")[[1]]
		NAA <- min( length(aaCall), round(NB/3))
		if ( ! DO_SUBSET) {
			for (k in 1:NAA) {
				kk <- (k-1) * 3 + 1 + AAoffset
				axis( side=1, at=peakPos[kk], label=aaCall[k], line=1, col.axis='black', col.ticks=NA, 
					font=2, lwd.ticks=0, cex.axis=cex*1.4)
			}
		} else {
			for (k in 1:NAA) {
				kk <- (k-1) * 3 + 1 + AAoffset
				thisX <- peakPos[kk]
				if (thisX < firstTracePoint || thisX > lastTracePoint) next
				axis( side=1, at=peakPos[kk], label=aaCall[k], line=1, col.axis='black', col.ticks=NA, 
					font=2, lwd.ticks=0, cex.axis=cex*1.4)
			}
		}
	}
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



chromatogramInspector <- function( isolateName, seq, type=c("DNA","AA"), allowPartialMatch=FALSE,
							min.scorePerAA=1, finger=c("DBL4","CLAG2"), ...) {

	require( "sangerseqR")
	type <- match.arg( type)
	finger <- match.arg( finger)

	chromoSet <- findSeqInChromatograms( isolateName, seq, type=type, 
						allowPartialMatch=allowPartialMatch, min.scorePerAA=min.scorePerAA, 
						finger=finger)
	ans <- plotChromatogramSet( chromoSet, isolateName, seq, type=type, ...)
	ans
}


findSeqInChromatograms <- function( isolateName, seq, type=c("DNA","AA"), allowPartialMatch=FALSE,
							min.scorePerAA=1, finger=c("DBL4","CLAG2")) {

	# given a (small) sequence, find and show that region in all the chromatograms that saw it
	isolateID <- sub( "_DBL4$", "", isolateName)
	isolateID <- sub( "_CLAG2$", "", isolateID)
	fragFile <- file.path( RESULTS.PATH, paste( isolateID, ".", finger, ".Fragment", type, ".fasta", sep=""))
	fa <- loadFasta( fragFile, verbose=T)
	nFrag <- length( fa$desc)

	# we will find this sequence in any fragments, and save up what we need to render them all afterward
	chromoList <- vector( mode="list")
	nChromo <- 0

	require( Biostrings)
	cat( "\nSearching..")
	for ( i in 1:nFrag) {
		thisFragID <- fa$desc[i]
		thisSeq <- fa$seq[i]
		# find where/if the given seq is in this fragment
		if (allowPartialMatch) {
			pa.type <- "local"
		} else {
			pa.type <- "global-local"
		}
		pa <- pairwiseAlignment( seq, thisSeq, type=pa.type, scoreOnly=F)
		if ( score( pa) < (nchar(seq)*min.scorePerAA)) next
		# yes, it's there, extract the location
		fragStart <- start( subject( pa))
		fragStop <- width( subject( pa)) + fragStart - 1

		if ( ! allowPartialMatch) {
			thisLen <- fragStop - fragStart + 1
			seqLen <- nchar( seq)
			if ( thisLen < seqLen) {
				cat( "\nFound, but too short:  ", thisLen, seqLen, thisFragID)
				next
			}
		}

		# grab that chromatogram
		# strip out any suffixes we added to the fragments name
		thisFragID <- sub( "_Fwd_F.+", "", thisFragID)
		thisFragID <- sub( "_Rev_F.+", "", thisFragID)
		chromoFile <- file.path( SEQUENCE.PATH, isolateName, paste( thisFragID, ".ab1", sep=""))
		if ( ! file.exists( chromoFile)) {
			cat( "\nCan not find/open chromatogram file:  ", chromoFile)
			next
		}
		ab1 <- sangerseq( read.abif( chromoFile))
		primary.seq <- as.character( primarySeq(ab1))
		traceM <- traceMatrix( ab1)
		peakPosM <- peakPosMatrix( ab1)
		peak.width <- round( nrow(traceM) / nrow(peakPosM))

		# turn this DNA string into whatever it may represent
		if ( type == "DNA") { 
			primary.set <- c( primary.seq, myReverseComplement(primary.seq))
			names(primary.set) <- c( "DNA", "rc(DNA)")
		} else {
			primary.set <- DNAtoAA( primary.seq, clipAtStop=F, readingFrame=1:6)
			names(primary.set) <- paste( "AA_Frame", 1:6, sep="")
		}

		# find this small seq in the full length chromatograms
		pa <- pairwiseAlignment( primary.set, seq, type="local-global", scoreOnly=T)
		best <- which.max( pa)
		bestScore <- pa[best]
		scorePerAA <- round( bestScore / nchar(seq), digits=2)
		pa <- pairwiseAlignment( primary.set[best], seq, type="local-global", scoreOnly=F)
		fragStart <- start( pattern( pa))
		fragStop <- width( pattern( pa)) + fragStart - 1
		fragString <- substr( primary.set[best], fragStart, fragStop)

		# given this best place, resolve any DNA & RevComp issues
		needRevComp <- FALSE
		if (type == "DNA" && best == 2) needRevComp <- TRUE
		if (type == "AA" && best > 3) needRevComp <- TRUE
		if (needRevComp) {
			traceM <- revCompTraceMatrix( traceM)
			peakPosM <- revCompPeakPosMatrix( peakPosM, nrow(traceM))
			thisFragID <- paste( "RevComp(", thisFragID, ")")
		}
		fragStartDNA <- fragStartAA <- fragStart
		fragStopDNA <- fragStopAA <- fragStop
		if (type == "DNA") fragStartAA <- fragStopAA <- NA
		if (type == "AA") {
			fragStartDNA <- (fragStartAA-1) * 3 + 1
			fragStopDNA <- fragStopAA * 3
			# fine tune the locations for reading frame
			frame.shift <- 0
			if ( best %in% c(2,5)) frame.shift <- 1
			if ( best %in% c(3,6)) frame.shift <- 2
			fragStartDNA <- fragStartDNA + frame.shift
			fragStopDNA <- fragStopDNA + frame.shift
		}

		# OK, we have all the details we need to draw this chromatogram fragment.
		# save it up so we can draw them all at once
		tmp <- list( "ID"=thisFragID, "AA_START"=fragStartAA, "AA_STOP"=fragStopAA, "DNA_START"=fragStartDNA, 
				"DNA_STOP"=fragStopDNA, "SEQ"=fragString, "TRACE_M"=traceM, "PEAK_M"=peakPosM, 
				"PEAK_WIDTH"=peak.width)
		nChromo <- nChromo + 1
		chromoList[[ nChromo]] <- tmp

		cat( "\nDebug:  ", i, thisFragID, "RevComp=", needRevComp, "best=", best, "Score=", bestScore, "PerAA=", scorePerAA,
				fragStartAA, fragStopAA)
	}
	cat( "  N_Found: ", nChromo)
	return( chromoList)
}


plotChromatogramSet <- function( chromoSet, isolateName, seq, type=type, lwd=2, cex=1, font=2) {

	# given a set of sub-regions in chromatograms, show them all in one plot
	nChromo <- length( chromoSet)
	if ( nChromo < 1) return()

	par( mfrow=c( (nChromo+1), 1))
	par( mai=c(0.5, 0.5, 0.4, 0.2))

	acgtBases <- c('A','C','G','T')
	acgtColors <- c('red','blue','orange','green')

	cat( "  Sizing..")
	consensusID <- consensusSize <- 1
	for ( i in 1:nChromo) {
		tmp <- chromoSet[[i]]
		thisFragID <- tmp$ID
		startAA <- tmp$AA_START
		stopAA <- tmp$AA_STOP
		startDNA <- tmp$DNA_START
		stopDNA <- tmp$DNA_STOP
		fragString <- tmp$SEQ
		traceM <- tmp$TRACE_M
		peakPosM <- tmp$PEAK_M
		half.peak.width <- round( tmp$PEAK_WIDTH/2)

		# the region we want is always DNA, plus a half peak width
		leftPos <- max( peakPosM[ startDNA, 1]-half.peak.width, 1)
		rightPos <- min( peakPosM[ stopDNA, 1]+half.peak.width, nrow(traceM))
		smallTraceM <- traceM[ leftPos:rightPos, ]
		smallPeakM <- peakPosM[ startDNA:stopDNA, ]

		# find the first one that is "full size", to set our limits, etc.
		if ( nrow(smallTraceM) > consensusSize) {
			# let's try to build a consensus chromatogram
			consensusTraceM <- smallTraceM
			rownames(consensusTraceM) <- leftPos : rightPos
			colnames(consensusTraceM) <- c("A","C","G","T")
			consensusSize <- nrow(smallTraceM)
			consensusPeakM <- smallPeakM
			consensusOffset <- leftPos
			consensusFirstPeak <- smallPeakM[1,1]
			consensusID <- i
			consensusLeftPos <- leftPos
			consensusRightPos <- rightPos
		}
	}

	cat( "  Plotting..")
	for ( i in 1:nChromo) {
		tmp <- chromoSet[[i]]
		thisFragID <- tmp$ID
		startAA <- tmp$AA_START
		stopAA <- tmp$AA_STOP
		startDNA <- tmp$DNA_START
		stopDNA <- tmp$DNA_STOP
		fragString <- tmp$SEQ
		traceM <- tmp$TRACE_M
		peakPosM <- tmp$PEAK_M
		half.peak.width <- round( tmp$PEAK_WIDTH/2)

		# the region we want is always DNA, plus a half peak width
		leftPos <- max( peakPosM[ startDNA, 1]-half.peak.width, 1)
		rightPos <- min( peakPosM[ stopDNA, 1]+half.peak.width, nrow(traceM))
		smallTraceM <- traceM[ leftPos:rightPos, ]
		smallPeakM <- peakPosM[ startDNA:stopDNA, ]

		# let's try to build a consensus chromatogram
		padleft <- padright <- FALSE
		peakCenterOffset <- 0
		if ( i != consensusID) {
			thisSize <- nrow(smallTraceM)

			# there is a chance that this chunk is not full length
			if ( nchar(fragString) < nchar(seq)) {
				pa <- pairwiseAlignment( fragString, seq, type="global-local")
				pattStart <- start( pattern( pa))
				subjStart <- start( subject( pa))
				if ( pattStart == 1 && subjStart > 1) {
					# we need to pad-fill at the left
					padleft <- TRUE
				} else {
					padright <- TRUE
				}
				nExtra <- consensusSize - nrow( smallTraceM)
				extraTraceM <- matrix( 0, nrow=nExtra, ncol=4)
				if (padleft) {
					smallTraceM <- rbind( extraTraceM, smallTraceM)
					leftPos <- leftPos
					rightPos <- rightPos + nExtra
					peakCenterOffset <- nrow( extraTraceM)
				} else {
					smallTraceM <- rbind( smallTraceM, extraTraceM)
					rightPos <- rightPos + nExtra
				}
			}

			tmpTraceM <- matrix( NA, nrow=consensusSize, ncol=4)
			scaleFactor <- nrow(smallTraceM) / consensusSize
			linearShift <- (smallPeakM[1,1] - leftPos) - (consensusFirstPeak - consensusOffset)
			cat( "\nDebug: scale, shift: ", scaleFactor, linearShift)
			tmpRowPtr <- round( 1:consensusSize * scaleFactor)
			tmpRowPtr <- tmpRowPtr + linearShift
			tmpRowPtr[ tmpRowPtr < 1] <- NA
			tmpRowPtr[ tmpRowPtr > consensusSize] <- NA
			for ( j in 1:4) tmpTraceM[ , j] <- smallTraceM[ tmpRowPtr, j]
			tmpTraceM[ is.na(tmpTraceM)] <- 0
			consensusTraceM <- consensusTraceM + tmpTraceM
		}

		plot( 1,1, type="n", main=thisFragID, xlim=c(leftPos,rightPos), ylim=c(0,max(smallTraceM)), ylab=NA, xlab=NA,
				xaxt="n", xaxs="i", las=2)
		for ( j in 1:4) {
			x <- leftPos : rightPos
			y <- smallTraceM[ ,j]
			lines( x, y, col=acgtColors[j], lwd=lwd)
		}
		peakCenters <- smallPeakM[ , 1]
		peakCentersOut <- smallPeakM[ , 1] + peakCenterOffset
		callBase <- acgtBases[ apply( traceM[ peakCenters, ], MARGIN=1, which.max)]
		colorBase <- acgtColors[ match( callBase, acgtBases)]
		for (k in 1:length(callBase)) axis( side=1, at=peakCentersOut[k], label=callBase[k], col.axis=colorBase[k], 
				col.ticks=colorBase[k], font=font, cex.axis=cex)

		callAA <- DNAVtoAAV( callBase)
		aaCenters <- peakCentersOut[ seq(2, length(peakCenters),by=3)]
		for ( k in 1:length(callAA)) axis( side=1, at=aaCenters[k], label=callAA[k], col=1, cex=1.4*cex, 
				font=font, line=1, tick=F)

	}

	# lastly show that consensus
	plot( 1,1, type="n", main="Consensus Chromatogram", xlim=c(1,consensusSize), ylim=c(0,max(consensusTraceM, na.rm=T)), 
			ylab=NA, xlab=NA, xaxt="n", xaxs="i", las=2)
	for ( j in 1:4) {
		x <- 1 : consensusSize
		y <- consensusTraceM[ ,j]
		lines( x, y, col=acgtColors[j], lwd=lwd)
	}
	peakCenters <- consensusPeakM[ , 1] - consensusOffset + 1
	callBase <- acgtBases[ apply( consensusTraceM[ peakCenters, ], MARGIN=1, which.max)]
	colorBase <- acgtColors[ match( callBase, acgtBases)]
	for (k in 1:length(callBase)) axis( side=1, at=peakCenters[k], label=callBase[k], col.axis=colorBase[k], 
			col.ticks=colorBase[k], font=font, cex.axis=cex)

	callAA <- DNAVtoAAV( callBase)
	aaCenters <- peakCenters[ seq(2, length(peakCenters),by=3)]
	for ( k in 1:length(callAA)) axis( side=1, at=aaCenters[k], label=callAA[k], col=1, cex=1.4*cex, 
			font=font, line=1, tick=F)

	# package up the consensus
	consensusTracePct <- consensusTraceM
	rowSums <- apply( consensusTraceM, 1, sum, na.rm=T)
	for (j in 1:4) consensusTracePct[ ,j] <- round( consensusTracePct[ ,j] * 100 / rowSums, digits=0)
	traceAns <- as.data.frame( consensusTracePct)
	traceAns$PeakCall <- ""
	traceAns$PeakCall[ peakCenters] <- callBase

	out <- list( "TracePcts"=traceAns, "DNA"=paste( callBase, collapse=""), "AA"=paste( callAA, collapse=""))
	out
}


`cleanChromatogramDNA` <- function( chromoObj, referenceDNA, referenceAA=NULL, nSkip=20, verbose=TRUE) {

	# we find cases where Sanger sequencing is inserting duplicate bases, which throw off reading frame
	# try to find and correct, using a reference DNA and perhaps AA sequence.
	require( Biostrings)

	# first efforts, quite subjective, currently only changing the DNA & AA fields
	#		not yet modifying the trace matrix or peak calls...
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	chromoDNA <- chromoObj$DNA_Calls

	# step 1: see which strand is better match to reference
	dnaScores <- pairwiseAlignment( chromoDNA, referenceDNA, type="global-local", scoreOnly=T)
	dna <- chromoDNA[ which.max(dnaScores)]

	# step 2: try to re-call "N" bases manually, but ignore the edges
	hasNcall <- which( names(peakPos) == "N")
	NP <- length(peakPos)
	hasNcall <- setdiff( hasNcall, c( 1:nSkip, (NP-nSkip+1):NP))
	if ( length( hasNcall)) {
		for (k in hasNcall) {
			baseAns <- baseCallOnePeak( peakPos[k], traceM)
			newBase <- names(baseAns)[1]
			if ( newBase != "N") {
				substr( dna, k, k) <- newBase
				names(peakPos)[k] <- newBase
			}
		}
	}

	# step 3: repeatedly see if we find a 1-2bp gap
	ans <- DNAtoCleanChromatogramDNA( dna, referenceDNA=referenceDNA, referenceAA=referenceAA, verbose=verbose)
	dnaNow <- ans$DNA
	details <- ans$CleaningDetails
	if ( ! is.null( details)) {
		# the details show where bases got added/removed...
		# we can update the names of the 'peakPos' data so future plots see the change
		for ( k in 1:nrow(details)) {
			changeLoc <- details$Location[k]
			baseWas <- details$PreviousBase[k]
			baseNew <- details$CleanedBase[k]
			chromoBase <- names( peakPos)[ changeLoc]
			# verify they say what they should
			if ( chromoBase != baseWas) cat( "\nDebug assertion failed: cleaning base error: ", chromoBase, baseWas)
			names(peakPos)[changeLoc] <- baseNew
		}
	}

	# repackage and return
	dnaOut <- c( dnaNow, myReverseComplement(dnaNow))
	names(dnaOut) <- c( "DNA", "rc(DNA)")
	aaOut <- DNAtoAA( dnaNow, clipAtStop=F, readingFrame=1:6)
	names(aaOut) <- paste( "AA_Frame", 1:6, sep="")

	out <- list( "TraceM"=traceM, "PeakPosition"=peakPos, "DNA_Calls"=dnaOut, "AA_Calls"=aaOut)
	out
}


# some routines from trying to quantify codon calls at specific mutation sites

# get the motif region directly from the ABI file
motifChromatogram <- function( f, motif, gene="", plot=TRUE, min.score=0, verbose=TRUE, referenceDNA=NULL) {

	motifName <- names(motif)[1]
	motif <- toupper( motif[1])
	if (is.null(motifName) || is.na(motifName) || motifName == "") motifName <- motif

	# given one ABI file, extract the subset that we want
	ch <- loadChromatogram( f)
	ch2 <- subsetChromatogram( ch, substring=motif)
	if ( ! is.null( referenceDNA)) ch2 <- cleanChromatogramDNA( ch2, referenceDNA=referenceDNA, nSkip=1, verbose=verbose)

	if (plot) {
		mainText <- paste( "  Motif:  ", motifName, "    Gene:  ", gene, 
					"\nFile: ", sub( ".ab1$","",basename(f)))
		plotChromatogram( ch2, label=mainText)
	}

	# sanity check to see if we got that motif
	aaSeqs <- ch2$AA_Calls[1:3]
	pa <- pairwiseAlignment( aaSeqs, motif, type="local-global", scoreOnly=T)
	best <- which.max( pa)
	if ( pa[best] < min.score) {
		if( verbose) {
			cat( "\nFailed to find motif in chromatogram.\nMotif=", motif, 
				"  File", basename(f), "  Score=", pa[best])
		}
		if (plot) {
			txt <- paste( "Failed to find motif: ", motif)
			usr <- par('usr')
			text( mean(usr[1:2]), mean(usr[3:4]), txt, pos=3, cex=1, font=2)
		}
		return( NULL)
	}
	return( ch2)
}


# zoom in on just the center codon
motifCodonChromatogram <- function( ch, motif, plot=TRUE) {

	# given the chromatogram subset that contains the motif, 
	# extract the exact center as the codon of interest
	if ( is.null(ch)) return(NULL)

	motifName <- names(motif)[1]
	motif <- toupper( motif[1])
	if (is.null(motifName) || is.na(motifName) || motifName == "") motifName <- motif

	aaSeqs <- ch$AA_Calls[1:3]
	pa <- pairwiseAlignment( aaSeqs, motif, type="local-global", scoreOnly=T)
	best <- which.max( pa)
	pa <- pairwiseAlignment( aaSeqs[best], motif, type="local-global", scoreOnly=F)
	subjStart <- start( subject( pa))
	pattStart <- start( pattern( pa))
	# use exactly the center
	nFlankAA <- floor( nchar(motif)/2)
	alleleLocAA <- pattStart + nFlankAA
	alleleLocDNA <- alleleLocAA * 3 - 2
	traceMstart <- ch$PeakPosition[alleleLocDNA]
	traceMstop <- ch$PeakPosition[alleleLocDNA+2]
	traceHalfWidth <- round( (traceMstop - traceMstart) / 4)
	traceMstart <- traceMstart - traceHalfWidth
	traceMstop <- traceMstop + traceHalfWidth
	smallTM <- ch$TraceM[ traceMstart:traceMstop, ]
	rownames(smallTM) <- 1:nrow(smallTM)

	# draw what/where we found it
	if (plot) {
		bigY <- max( smallTM, na.rm=T)
		rect( traceMstart, -bigY/30, traceMstop, bigY, border='black', lwd=4)
	}

	# extract just that region to send back
	smallPeaks <- ch$PeakPosition[ alleleLocDNA:(alleleLocDNA+2)]
	smallDNA <- substr( ch$DNA_Calls[1], alleleLocDNA, alleleLocDNA+2)
	smallAA <- substr( ch$AA_Calls[best], alleleLocAA, alleleLocAA)

	# adjust the peak centers to these new local units
	smallPeaks <- smallPeaks - traceMstart + 1

	# send this tiny chromatogram back
	out <- list( "TraceM"=smallTM, "PeakPosition"=smallPeaks, "DNA_Calls"=smallDNA, "AA_Calls"=smallAA)
	return( out)
}


measureDominantCodons <- function( ch, nCodons=1, plot=T) {

	# start with what is given
	peaks <- ch$PeakPosition
	tm <- ch$TraceM
	calledAA <- substr( ch$AA_Calls[1], 1, 1)

	# use just the center top region to ignore the tails
	# mask out the tails
	tmIn <- tm
	halfWidth <- 2
	maskM <- matrix( 0, nrow=nrow(tm), ncol=ncol(tm))
	for ( p in peaks) maskM[ (p-halfWidth):(p+halfWidth), ] <- 1
	tm <- tm * maskM
	tmFull <- tm
	total <- apply( tmFull, MARGIN=1, sum, na.rm=T)
	
	# we may be asked to call more than one, so set up to loop around
	aaOut <- dnaOut <- rep.int( "", nCodons)
	rankOut <- pctOut <- rep.int( NA, nCodons)
	myPct <- vector( length=3)
	for ( i in 1:nCodons) {
	
		# grab the called base data and the overall sums
		big <- apply( tm, MARGIN=1, max, na.rm=T)
		calledBase <- apply( tm[ peaks, ], MARGIN=1, function(x) colnames(tm)[ which.max(x)])
		calledAA <- DNAtoAA( paste(calledBase,collapse=""), readingFrame=1)
		for ( ip in 1:3) {
			p <- peaks[ip]
			px <- (p-halfWidth):(p+halfWidth)
			myPct[ip] <- sum(big[px]) * 100 / sum(total[px])
		}
		# the percentage of the dominant allele is from the lowest of the 3 bases!
		whoLow <- which.min( myPct)
		pctAns <- round( myPct[whoLow])
		lowPeakHeight <- tm[ peaks[whoLow], calledBase[whoLow]]

		# store the answer for this codon call, with a lower limit on how much
		# we can trust what we found
		if (pctAns < 10) break
		aaOut[i] <- calledAA
		dnaOut[i] <- paste( calledBase, collapse="")
		pctOut[i] <- pctAns
		rankOut[i] <- i

		# if we need more than one call, subtract out what we called and repeat
		if ( i < nCodons) {
			for ( j in 1:3) {
				myBase <- calledBase[j]
				p <- peaks[j]
				px <- (p-halfWidth):(p+halfWidth)
				newValues <- tm[ px, myBase] - lowPeakHeight
				newValues <- pmax( newValues, 0)
				tm[ px, myBase] <- newValues
			}
		}
	}

	# only keep/show what we filled in
	nFound <- sum( aaOut != "")
	length(aaOut) <- length(dnaOut) <- length(pctOut) <- length(rankOut) <- nFound
	if ( plot) {
		usr <- par( "usr")
		myX <- mean(usr[1:2])
		myY <- max(tmFull, na.rm=T)
		txt <- paste( aaOut, "=", pctOut, "%", sep="", collapse="; ")
		text( myX, myY, txt, font=2, cex=1.4, pos=3)
	}
	return( data.frame( "Status"="Pass", "AA_Call"=aaOut, "Codon"=dnaOut, "Percent"=pctOut, "Rank"=rankOut,
					stringsAsFactors=F))
}


