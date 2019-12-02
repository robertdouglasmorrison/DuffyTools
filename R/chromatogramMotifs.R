# chromatogramMotifs.R -- pieces to find known mutation site regions
#			and extract/quantify the mutaions calls

# some routines from trying to quantify codon calls at specific mutation sites
# we use the term "motif" to mean a ~11aa short sequence of context that is 
# centered on a single AA of interest, and we use the term "codon" as the center AA of interest

# get the motif region directly from an ABI chromatogram file
`motifSubsetChromatogram` <- function( f, motif, gene="", plot=TRUE, min.score=0, verbose=TRUE, referenceDNA=NULL) {

	motifName <- names(motif)[1]
	motif <- toupper( motif[1])
	if (is.null(motifName) || is.na(motifName) || motifName == "") motifName <- motif

	# given one ABI file, extract the subset that we want
	ch <- loadChromatogram( f)
	ch2 <- subsetChromatogram( ch, seq=motif)
	if ( ! is.null( referenceDNA)) {
		ch2 <- cleanChromatogramDNA( ch2, referenceDNA=referenceDNA, nSkip=1, verbose=verbose)
	}

	if (plot) {
		mainText <- paste( "  Motif:  ", motifName, "    Gene:  ", gene, 
					"\nFile: ", sub( ".ab1$","",basename(f)))
		plotChromatogram( ch2, label=mainText, showConfidence=T)
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


# zoom in on just the center codon of a motif, and return it as a chromatogram data object
`motifCodonChromatogram` <- function( chromoObj, motif, plot=TRUE) {

	# given the chromatogram subset that contains the motif, 
	# extract the exact center as the codon of interest
	if ( is.null(chromoObj)) return(NULL)

	motifName <- names(motif)[1]
	motif <- toupper( motif[1])
	if (is.null(motifName) || is.na(motifName) || motifName == "") motifName <- motif

	require( Biostrings)
	data( BLOSUM62)
	aaSeqs <- chromoObj$AA_Calls[1:3]
	motif <- gsub( "?", "X", motif, fixed=T)
	pa <- pairwiseAlignment( aaSeqs, motif, type="local-global", substitutionMatrix=BLOSUM62, scoreOnly=T)
	best <- which.max( pa)
	pa <- pairwiseAlignment( aaSeqs[best], motif, type="local-global", substitutionMatrix=BLOSUM62, scoreOnly=F)
	subjStart <- start( subject( pa))
	pattStart <- start( pattern( pa))
	# use exactly the center
	nFlankAA <- floor( nchar(motif)/2)
	alleleLocAA <- pattStart + nFlankAA
	alleleLocDNA <- alleleLocAA * 3 - 2
	traceMstart <- chromoObj$PeakPosition[alleleLocDNA]
	traceMstop <- chromoObj$PeakPosition[alleleLocDNA+2]
	traceHalfWidth <- round( (traceMstop - traceMstart) / 4)
	traceMstart <- traceMstart - traceHalfWidth
	traceMstop <- traceMstop + traceHalfWidth
	smallTM <- chromoObj$TraceM[ traceMstart:traceMstop, ]
	rownames(smallTM) <- 1:nrow(smallTM)

	# draw what/where we found it
	if (plot) {
		bigY <- max( smallTM, na.rm=T)
		rect( traceMstart, -bigY/30, traceMstop, bigY, border='black', lwd=4)
	}

	# extract just that region to send back
	smallPeaks <- chromoObj$PeakPosition[ alleleLocDNA:(alleleLocDNA+2)]
	smallConf <- chromoObj$PeakConfidence[ alleleLocDNA:(alleleLocDNA+2)]
	smallDNA <- substr( chromoObj$DNA_Calls[1], alleleLocDNA, alleleLocDNA+2)
	smallAA <- substr( chromoObj$AA_Calls[best], alleleLocAA, alleleLocAA)

	# adjust the peak centers to these new local units
	smallPeaks <- smallPeaks - traceMstart + 1

	# send this tiny chromatogram back
	out <- list( "TraceM"=smallTM, "PeakPosition"=smallPeaks, "PeakConfidence"=smallConf,
				"DNA_Calls"=smallDNA, "AA_Calls"=smallAA, "Filename"=chromoObj$Filename)
	return( out)
}


`measureDominantCodons` <- function( chromoObj, nCodons=1, plot=T) {

	# start with what is given
	peaks <- chromoObj$PeakPosition
	tm <- chromoObj$TraceM
	calledAA <- substr( chromoObj$AA_Calls[1], 1, 1)

	# use just the center top region to ignore the tails
	# mask out the tails
	tmIn <- tm
	halfWidth <- 2
	maskM <- matrix( 0, nrow=nrow(tm), ncol=ncol(tm))
	for ( p in peaks) maskM[ (p-halfWidth):(p+halfWidth), ] <- 1
	tm <- tm * maskM
	tmFull <- tm
	total <- apply( tmFull, MARGIN=1, sum, na.rm=T)
	
	# we may be asked to call more than one codon at this location, so set up to loop around
	# repeatedly: make the dominant call, subtract that out, and repeat on the residual
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


# quantify how confident the dominant call is for a region of peaks
# note that confidence is now an explicit element of the chromatogram object
# returns a value in 0 to 100 range
`motifCodonConfidence` <- function( chromoObj, centerPeak, nNeighbors=0) {

	# get the raw chromatogram data and the peak list
	if ( is.null(chromoObj)) return( NULL)
	peaks <- chromoObj$PeakPosition
	confs <- chromoObj$PeakConfidence

	# the peak we were told to center on must be a real peak
	myPeakPtr <- match( centerPeak, peaks, nomatch=0)
	if ( !myPeakPtr) {
		cat( "\nGiven Peak location is not in the chromatogram peak set: ", centerPeak)
		cat( "\nNot one of: ", peaks)
		return( NULL)
	}

	# the overall confidence for the interval is the average of all the peak confidence calls
	leftPtr <- max( myPeakPtr - nNeighbors, 1)
	rightPtr <- min( myPeakPtr + nNeighbors, length(peaks))
	# show the region of the trace matrix that this covers
	leftPt <- peaks[ leftPtr] - 5
	rightPt <- peaks[ rightPtr] + 5
	avgConf <- mean( confs[ leftPtr:rightPtr], na.rm=T)
	# scale to 0..100
	avgConf <- round( avgConf * 100, digits=2)
	return( list( "Confidence"=avgConf, "StartPoint"=leftPt, "StopPoint"=rightPt))
}


# higher level call that gets the confidence and shows it on the current chromatogram plot
`measureMotifConfidence` <- function( chromoObj, nNeighbors=7, plot=T) {

	# start with what is given
	peaks <- chromoObj$PeakPosition

	# use just the center top region to ignore the tails
	centerPeak <- peaks[ floor( (length(peaks)+1) / 2)]
	ans <- motifCodonConfidence( chromoObj, centerPeak=centerPeak, nNeighbors=nNeighbors)
	if ( is.null(ans)) return(0)

	conf <- ans$Confidence

	if ( plot) {
		# draw something?...
		tm <- chromoObj$TraceM
		xl <- ans$StartPoint
		xh <- ans$StopPoint
		# use the exact same scaling as the plotChromatogram tool
		# then scale by the confidence
		bigY <- max( tm, na.rm=T) * 1.1
		scaledY <- bigY / 100 * conf
		rect( xl, 0, xh, scaledY, border=1, lwd=1)
		text( xh, scaledY*0.98, paste( "Confidence = ", conf, "%", sep=""), font=2, cex=0.95, pos=if(conf > 90) 1 else 3)
	}
	return( conf)
}
