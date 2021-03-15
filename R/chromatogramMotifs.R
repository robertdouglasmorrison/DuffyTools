# chromatogramMotifs.R -- pieces to find known mutation site regions
#			and extract/quantify the mutaions calls

# some routines from trying to quantify codon calls at specific mutation sites
# we use the term "motif" to mean a ~11aa short sequence of context that is 
# centered on a single AA of interest, and we use the term "codon" as the center AA of interest

# get the motif region directly from an ABI chromatogram file
`motifSubsetChromatogram` <- function( f, motif, gene="", plot=TRUE, min.score=0, verbose=TRUE, 
					referenceDNA=NULL, ...) {

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
		plotChromatogram( ch2, label=mainText, showConfidence=T, ...)
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
	cat( "\nDebug old codon measure tool: ", nFound, aaOut, pctOut)
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



`modelFitOneCodon` <- function( chromoObj, ref.AA, alt.AA, min.percent=1, plot=T, visualize.fit=F, verbose=T) {

	# given a tiny chromatogram that spans just a single codon, fit this to the
	# expected amino acid calls being searched for.  The fit tool will search for any codons
	# that code for the reference AA, the expected alternate AA (like a resistance mutation)
	# and the codon AA as called by the chromatogram itself

	failAns <- data.frame( "Status"="FAIL", "Best_AA_Call"="", "Best_Codon"="", "Best_Percent"=NA, 
				"Ref_Name"="", "Ref_Percent"=NA, "Mutant_Name"="", 
				"Mutant_Percent"=NA, "Confidence"=0, stringsAsFactors=F)
	if ( is.null( chromoObj)) return( failAns)
	
	# start with what is given
	peaks <- chromoObj$PeakPosition
	conf <- chromoObj$PeakConfidence
	tm <- chromoObj$TraceM
	calledAA <- substr( chromoObj$AA_Calls[1], 1, 1)
	
	# do a simple fit to get the optimal peak width
	fit0 <- modelFitChromatogram( chromoObj, fixedPeaks=TRUE, 
				effort=1, doStandardize=T, doSubset=F, algorithm="GenSA")
	bestPeakWidth <- fit0$FitPeakWidth

	# we will use the model fit tools that expect the model sequences to be explicitly given.
	# use the codon map and trim down to just the non-degenerate AA of interest
	cmap <- getCodonMap()
	cmap <- cmap[ grep( "^[ACGT]+$", cmap$DNA), ]
	myCmap <- subset( cmap, AA %in% c( ref.AA, alt.AA, calledAA))
	
	# name these as being 'ref' or 'mutant'
	modelSeqs <- myCmap$DNA
	prefix <- c( "Ref", "Mutant", "Other")[ match( myCmap$AA, c( ref.AA, alt.AA, calledAA))]
	names(modelSeqs) <- paste( prefix, myCmap$AA, myCmap$DNA, sep="_")
	
	# there is a small chance that the given pool of Ref/Mutant/Called codons do not fully capture
	# the raw data.  Check each of the 3 bases in the codon for minor peaks, and perhaps extend our
	# set of models
	for (k in 1:3) {
		givenCodon <- chromoObj$DNA_Calls[1]
		tmPoint <- peaks[k]
		tmValues <- apply( tm[ (tmPoint-1):(tmPoint+1), ], MARGIN=2, sum)
		baseOrd <- order( tmValues, decreasing=T)
		base2 <- colnames(tm)[ baseOrd[2]]
		heights <- tmValues[ baseOrd[1:2]]
		if (heights[2] >= heights[1]/2) {
			# got a second peak at this base
			extraCodon <- givenCodon
			substr( extraCodon, k, k) <- substr(base2,1,1)
			# have we already yet seen this codon?
			if ( extraCodon %in% myCmap$DNA) next
			moreCmap <- subset( cmap, DNA == extraCodon)
			if (nrow(moreCmap)) {
				names(extraCodon) <- paste( "Other", moreCmap$AA[1], moreCmap$DNA[1], sep="_")
				modelSeqs <- c( modelSeqs, extraCodon)
			}
		}
	}
	if (verbose) cat( "\nN_Model_Seq: ", length(modelSeqs), "\n  ", names(modelSeqs))
	
	# perhaps generate images of the codon fit inputs
	if ( visualize.fit) plotCodonFitInputs( chromoObj, modelSeqs, width=bestPeakWidth)
	
	# call the modeler
	ans <- modelBlendChromatogram( chromoObj, modelSeqs, plot.chromatograms=F, 
				synthetic.width=bestPeakWidth, noise.seqs=FALSE, verbose=verbose)
	if ( is.null( ans)) return( NULL)
	fitAns <- ans$Model.Estimates
	if ( is.null( fitAns)) return( failAns)
	
	# perhaps generate images of the codon fit inputs
	if ( visualize.fit) plotCodonFitResults( chromoObj, modelSeqs, fitAns, width=bestPeakWidth)
	
	# recombine what we got back into a format that best answers the question about what we saw
	# the P-values are useless since the data is so small we have no power.  
	# Just look at % calls of those that are above zero
	fitAns <- subset( fitAns, Percentage > 0)
	ord <- order( fitAns$Percentage, decreasing=T)
	fitAns <- fitAns[ ord, ]
	bestConstruct <- fitAns$Construct[1]
	# these are <prefix>_<AA>_<codon>; so keep the prefix and AA as the 'Construct' ID
	constructIDs <- sub( "_[ACGT]+$", "", fitAns$Construct)
	bestID <- constructIDs[1]
	bestCodon <- sub( ".+_","", bestConstruct)
	bestAA <- sub( "(^.+_)([A-Z])(_.+$)", "\\2", bestConstruct)
	
	# tally up how much percentage went with each ID, and renormalize to 100%
	amounts <- tapply( fitAns$Percentage, factor( constructIDs), sum, na.rm=T)
	pcts <- round( amounts * 100 / sum(amounts), digits=2)
	# enforce a minimum percentage to prevent small noise from acting
	# like a true minor variant.  Is so, renormalize again.
	drops <- which( pcts < min.percent)
	if ( length( drops)) {
		pcts <- pcts[ -drops]
		pcts <- round( pcts * 100 / sum(pcts), digits=2)
	}
	nCalls <- length(pcts)
	pcts <- round( pcts, digits=0)
	
	# how to decipher depends on how many came back
	if ( nCalls > 1) {
		# put them is pct sorted order
		pcts <- sort( pcts, decreasing=T)
		# we have 2+, so find the best of both groups: Ref & Not Ref
		isREF <- grep( "^Ref", names(pcts))
		isMUT <- grep( "^Mut", names(pcts))
		isOTH <- grep( "^Oth", names(pcts))
		if ( length(isREF)) {
			bestRefPtr <- isREF[ which.max( pcts[isREF])]
			refPct <- pcts[ bestRefPtr]
			refName <- names(refPct)
		} else {
			refPct <- 0
			refName <- ""
		}
		if ( length(isMUT)) {
			bestMutPtr <- isMUT[ which.max( pcts[isMUT])]
			mutPct <- pcts[ bestMutPtr]
			mutName <- names(mutPct)
		} else {
			mutPct <- 0
			mutName <- ""
		}
		if ( length(isOTH)) {
			bestOthPtr <- isOTH[ which.max( pcts[isOTH])]
			othPct <- pcts[ bestOthPtr]
			othName <- names(othPct)
		} else {
			othPct <- 0
			othName <- ""
		}
		# lastly with the top 2 known, set the best percentage
		bestPct <- if ( bestID == refName) refPct else mutPct
		if (othPct > bestPct) {
			bestPct <- othPct
		}
	} else {
		# with just the one, which is it...
		if ( bestAA == ref.AA) {
			refPct <- pcts[1]
			refName <- names(pcts)[1]
			mutPct <- 0
			mutName <- ""
			othPct <- 0
			othName <- ""
		} else {
			refPct <- 0
			refName <- ""
			mutPct <- othPct <- 0
			mutName <- othName <- ""
			if (grepl("^Oth", names(pcts[1]))) {
				othPct <- pcts[1]
				othName <- names(pcts)[1]
			} else {
				mutPct <- pcts[1]
				mutName <- names(pcts)[1]
			}
		}
		bestPct <- pcts[1]
	}
	
	# the final confidence is the average confidence of the codon bases
	confOut <- round( mean(conf,na.rm=T) * 100, digits=2)
	
	out <- data.frame( "Status"="Pass", "Best_AA_Call"=bestAA, "Best_Codon"=bestCodon, "Best_Percent"=bestPct, 
				"Ref_Name"=refName, "Ref_Percent"=refPct, "Mutant_Name"=mutName, 
				"Mutant_Percent"=mutPct, "Confidence"=confOut, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	
	if ( plot) {
		usr <- par( "usr")
		myX <- mean(usr[1:2])
		myY <- max(tm, na.rm=T)
		aaOut <- c( refName, mutName, othName)
		pctOut <- as.numeric( c( refPct, mutPct, othPct))
		# we only report those above zero
		keep <- which( pctOut >= 1)
		aaOut <- aaOut[keep]
		pctOut <- pctOut[keep]
		txt <- paste( aaOut, "=", pctOut, "%", sep="", collapse="; ")
		text( myX, myY, txt, font=2, cex=1.1, pos=3, offset=0.25)
	}
	
	return( out)
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


`plotCodonFitInputs` <- function( chromoObj, modelSeqs, width=2.5) {
	
	# given a tiny chromatogram of a single codon, and all the candidate codons that will go into a fit
	# try to make one image that conveys it all
	Nseq <- length( modelSeqs)
	Nplots <- Nseq + 1 + 1 # leave a space for the final fit residual
	sqrtN <- ceiling( sqrt(Nplots))
	
	# use a new window
	X11( type="dbcairo", bg='white', width=12, height=8)
	
	# set up to draw many small plots
	savMAI <- par('mai')
	on.exit( par(mai=savMAI))
	par( mai=c(0.6,0.7,0.5,0.2))
	mf <- c( sqrtN, sqrtN)
	if ( Nplots < 21) mf <- c(4,5)
	if ( Nplots < 17) mf <- c(4,4)
	if ( Nplots < 13) mf <- c(3,4)
	if ( Nplots < 10) mf <- c(3,3)
	if ( Nplots < 7) mf <- c(2,3)
	if ( Nplots < 5) mf <- c(2,2)
	par( mfcol=mf)
	
	# start with the observed data
	plotChromatogram( chromoObj, label="Observed_Codon", showAA=T, cex=1.4, cex.main=1.25, lwd=3, shiftAA=2,
			main.prefix="")

	# now make and draw every model
	bigHeight <- max( chromoObj$TraceM)
	for ( i in 1:Nseq) {
		thisSeq <- modelSeqs[i]
		thisName <- names( modelSeqs)[i]
		thisObj <- syntheticChromatogram( thisSeq, height=bigHeight, width=width)
		plotChromatogram( thisObj, label=thisName, showAA=T, cex=1.4, cex.main=1.25, lwd=3, shiftAA=2,
				main.prefix="Model:  ")
	}
	dev.flush()
	Sys.sleep(1)
	dev.print( pdf, "Chromatogram.Codon.Fit.Inputs.pdf", width=12)
	dev.off()
}


`plotCodonFitResults` <- function( chromoObj, modelSeqs, fitAns, width=2.5) {
	
	# given a tiny chromatogram of a single codon, and all the candidate codons that will go into a fit
	# try to make one image that conveys it all
	Nseq <- length( modelSeqs)
	Nplots <- Nseq + 1 + 1  # leava space for the final fit residual
	sqrtN <- ceiling( sqrt(Nplots))
	
	# use a new window
	X11( type="dbcairo", bg='white', width=12, height=8)
	
	# set up to draw many small plots
	savMAI <- par('mai')
	on.exit( par(mai=savMAI))
	par( mai=c(0.6,0.7,0.5,0.2))
	mf <- c( sqrtN, sqrtN)
	if ( Nplots < 21) mf <- c(4,5)
	if ( Nplots < 17) mf <- c(4,4)
	if ( Nplots < 13) mf <- c(3,4)
	if ( Nplots < 10) mf <- c(3,3)
	if ( Nplots < 7) mf <- c(2,3)
	if ( Nplots < 5) mf <- c(2,2)
	par( mfcol=mf)

	bigHeight <- max( chromoObj$TraceM)
	yMax <- bigHeight * 1.1
	
	# start with the observed data
	plotChromatogram( chromoObj, label="Observed_Codon", showAA=T, cex=1.4, cex.main=1.25, lwd=3, shiftAA=2,
			forceYmax=yMax, main.prefix="")
			
	# save a copy for showing the final residual
	residObj <- standardizeChromatogram( chromoObj, call.Ns=F)
	residTM <- obsTM <- residObj$TraceM

	# grab the model estimates, we will use the percentages.  Put them back into input order
	ord <- match( names(modelSeqs), fitAns$Construct)
	fitAns <- fitAns[ ord, ]
	
	# now make and draw every model
	for ( i in 1:Nseq) {
		thisSeq <- modelSeqs[i]
		thisPct <- fitAns$Percentage[i]
		thisName <- paste( names( modelSeqs)[i], "  ", round(thisPct,digits=1), "%", sep="")
		thisObj <- syntheticChromatogram( thisSeq, height=(bigHeight*thisPct/100), width=width)
		plotChromatogram( thisObj, label=thisName, showAA=(thisPct>=0.5), cex=1.4, cex.main=1.25, lwd=3, shiftAA=2,
				forceYmax=yMax, main.prefix="Fit:  ")
		if ( thisPct > 0) {
			thisTM <- thisObj$TraceM
			residTM <- residTM - thisTM
			residTM[ residTM < 0] <- 0
		}
	}
	
	# lastly show the final residual
	residObj$TraceM <- residTM
	residPct <- sum(residTM) * 100 / sum(obsTM)
	residName <- paste( "Residual after Fit  ", round(residPct,digits=1), "%", sep="")
	plotChromatogram( residObj, label=residName, showAA=F, cex=1.4, cex.main=1.25, lwd=3, shiftAA=2,
			forceYmax=yMax, main.prefix="")
	# done
	dev.flush()
	Sys.sleep(1)
	dev.print( pdf, "Chromatogram.Codon.Fit.Results.pdf", width=12)
	dev.off()
}

