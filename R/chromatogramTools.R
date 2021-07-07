# chromatogramTools.R -- pieces to work directly from ABI chromatogram files
#			and to manipulate the underlaying A/C/G/T spectra


`loadChromatogram` <- function( chromoFile, curated=TRUE) {

	if ( curated) {
		curatedFile <- sub( "\\.ab1$", ".rda", chromoFile)
		if ( file.exists( curatedFile)) {
			return( loadCuratedChromatogram( curatedFile))
		}
	}
	return( loadABiChromatogram( chromoFile))
}


`loadCuratedChromatogram` <- function( chromoFile) {

	curatedFile <- sub( "\\.ab1$", ".rda", chromoFile)
	if ( !file.exists( curatedFile)) {
		cat( "\nCurated Chromatogram file not found:  ", curatedFile)
		return( NULL)
	}
	
	# load that R data object into local environment
	chromoObj <- NULL
	load( curatedFile, envir=environment())
	if ( is.null( chromoObj)) {
		cat( "\nLoading Curated Chromatogram object failed:  ", curatedFile)
		return( NULL)
	}
	return( chromoObj)
}


`loadABiChromatogram` <- function( chromoFile) {

	# read up an AB1 file from disk... into our standard format of trace, peaks, & sequence
	require( "sangerseqR")

	if ( ! file.exists( chromoFile)) {
		cat( "\nChromatogram file not found:  ", chromoFile)
		return( NULL)
	}
	ab1 <- sangerseq( read.abif( chromoFile))

	# the raw data is 3 parts:   a DNA sequence,  a 4x matrix of base intensities,  
	# and a matrix of peak locations
	dna.seq <- as.character( primarySeq(ab1))
	traceM <- traceMatrix( ab1)
	peakPosM <- peakPosMatrix( ab1)

	# we will always keep both strands, as both DNA and AA
	seqData <- chromatogramSequences( dna.seq)
	
	# for now, just take the primary peak positions as ABI called them
	peakPos <- as.numeric( peakPosM[ ,1])
	names(peakPos) <- strsplit( dna.seq, split="")[[1]]

	# use those peak top locations to estimate where the trace stops being useful
	avgPeakSeparation <- median( diff( peakPos))
	rightTail <- max( peakPos) + avgPeakSeparation
	if ( rightTail < nrow(traceM)) {
		traceM <- traceM[ 1:rightTail, ]
	}

	# there could also be rare cases where the trace matrix does not contain all the peaks.
	# this should never happen, but check and fix if it does
	NTR <- nrow( traceM)
	maxPeakLoc <- NTR - floor( avgPeakSeparation/2)
	if ( any( peakPos > maxPeakLoc)) {
		# throw away those peaks
		drops <- which( peakPos > maxPeakLoc)
		peakPos <- peakPos[ -drops]
		# this means we need to truncate the DNA calls too
		if ( nchar(dna.seq) > length(peakPos)) {
			dna.seq <- substr( dna.seq, 1, length(peakPos))
			seqData <- chromatogramSequences( dna.seq)
		}
	}

	colnames(traceM) <- c( "A", "C", "G", "T")
	rownames(traceM) <- 1:nrow(traceM)
	
	# lastly, calculate confidence scores for each peak
	peakConfidence <- calcChromatogramPeakConfidence( traceM, peakPos)

	out <- list( "TraceM"=traceM, "PeakPosition"=peakPos, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoFile)
	out
}


`loadMultipleChromatograms` <- function( chromoFiles, curated=TRUE) {

	out <- list()
	n <- 0
	for ( f in chromoFiles) {
		ans <- loadChromatogram( f, curated=curated)
		if (is.null(ans)) next
		n <- n + 1
		out[[n]] <- ans
		# make a provisional name from the filename
		chromoName <- sub( "\\.ab1$", "", basename(f))
		names(out)[n] <- chromoName
	}
	if ( n == 0) return( NULL)
	return( out)
}


`writeCuratedChromatogram` <- function( chromoObj, chromoFile=NULL) {

	# filename is already in the object
	curatedFile <- sub( "\\.ab1$", ".rda", chromoObj$Filename)
	csvFile <- sub( "\\.ab1$", ".csv", chromoObj$Filename)
	if ( ! is.null( chromoFile)) {
		curatedFile <- sub( "\\.ab1$", ".rda", chromoFile)
		csvFile <- sub( "\\.ab1$", ".csv", chromoFile)
		chromoObj$Filename <- chromoFile
	}

	# save both as R object and as .CSV table
	save( chromoObj, file=curatedFile)
	write.table( chromatogramToTable( chromoObj), csvFile, sep=",", quote=T, row.names=F, na="")
}

	
`chromatogramSequences` <- function( dna.seq) {

	# given one DNA seq as called by the ABI tools, expand it to cover all
	# ways we may want to see that sequence
	DNA.set <- c( dna.seq, myReverseComplement(dna.seq))
	names(DNA.set) <- c( "DNA", "RevComp(DNA)")
	AA.set <- DNAtoAA( dna.seq, clipAtStop=F, readingFrame=1:6)
	names(AA.set) <- paste( "AA_Frame", 1:6, sep="")
	
	# clean the sequences to remove any unexpected characters that may break
	# sequence search/compare/score tools later
	DNA.set <- gsub( "?", "N", DNA.set, fixed=T)
	AA.set <- gsub( "?", "X", AA.set, fixed=T)

	return( list( "DNA"=DNA.set, "AA"=AA.set))
}


`revCompChromatogram` <- function( chromoObj) {

	traceM <- revCompTraceMatrix( chromoObj$TraceM)
	peakPos <- revCompPeakPosition( chromoObj$PeakPosition, nrow(traceM))
	peakConfidence <- rev( chromoObj$PeakConfidence)

	# since we already have the RevComp of the DNA, just use it
	dna.seq <- chromoObj$DNA_Calls[2]
	seqData <- chromatogramSequences( dna.seq)
	
	names(peakPos) <- strsplit( dna.seq, split="")[[1]]

	out <- list( "TraceM"=traceM, "PeakPosition"=peakPos, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename)
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
	rownames(out) <- 1:N
	out
}


`revCompPeakPosition` <- function( peakPos, nrow.TM) {

	# do the reverse of where the peaks are
	out <- rev( nrow.TM - peakPos + 1)
	out
}


`baseCallOnePeak` <- function( peakSite, traceM, min.pct=0.10, extra.width=1) {

	# use the data in a tiny window, not just the 1 point
	intenV <- apply( traceM[ (peakSite-extra.width):(peakSite+extra.width), , drop=FALSE], MARGIN=2, sum)
	intenPcts <- intenV / sum(intenV)
	# good to call if at least some % of total intensity
	good <- which( intenPcts >= min.pct)
	outBase <- colnames(traceM)[ good]
	if ( length(outBase)) {
		outPcts <- round( intenPcts[ good], digits=2)
		names(outPcts) <- outBase
		ord <- order( outPcts, decreasing=T)
	} else {
		outPcts <- 0
		names(outPcts) <- "N"
		ord <- 1
	}
	out <- outPcts[ord]
	return( out)
}


`baseCallProfile` <- function( chromoObj, min.pct=0.10, extra.width=1) {

	traceM <- chromoObj$TraceM
	peaks <- chromoObj$PeakPosition
	NP <- length( peaks)

	out <- lapply( peaks, baseCallOnePeak, traceM=traceM, min.pct=min.pct, extra.width=extra.width)
	names(out) <- paste( 1:NP, names(peaks), sep="_")

	return( out)
}


`fixChromatogramNcalls` <- function( chromoObj, verbose=FALSE) {

	# brute force fix any 'N' calls to be whatever is highest
	traceM <- chromoObj$TraceM
	peaks <- chromoObj$PeakPosition
	dna <- chromoObj$DNA_Calls[1]
	NP <- length( peaks)

	# re-call each N site, and update if we get a valid call
	isN <- which( names(peaks) == "N")
	if ( length( isN)) {
		nFixed <- 0
		for ( i in isN) {
			thisPeakPos <- peaks[i]
			thisBase <- names(peaks)[i]
			ans <- baseCallOnePeak( thisPeakPos, traceM=traceM, min.pct=0.1, extra.width=0)
			newBase <- names(ans)[1]
			if ( newBase != thisBase) {
				names(peaks)[i] <- newBase
				substr( dna, i, i) <- newBase
				nFixed <- nFixed + 1
			}
		}
		if (verbose) cat( "\nFixed 'N' base calls: ", basename(chromoObj$Filename), " \tWas:", 
				length(isN), " \tNow:", (length(isN)-nFixed))
	}

	# use the new data to rebuild what we need
	seqData <- chromatogramSequences( dna)

	# with better base calls, re-calculate confidence scores for each peak
	peakConfidence <- calcChromatogramPeakConfidence( traceM, peaks)

	out <- list( "TraceM"=traceM, "PeakPosition"=peaks, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename)
	out
}


`baseCallVerify` <- function( chromoObj) {

	# the peak calls 'SHOULD' match the Trace Matrix...
	# verify that explicitly and report the percentage agreement
	traceM <- chromoObj$TraceM
	peaks <- chromoObj$PeakPosition
	NP <- length( peaks)

	# see which base is the highest at the peak centerpoints
	centerM <- traceM[ as.numeric(peaks), ]
	centerBase <- colnames(traceM)[ apply( centerM, MARGIN=1, which.max)]
	calledBase <- as.character( names(peaks))

	nMatch <- sum( centerBase == calledBase)
	nNotMatch <- sum( centerBase != calledBase)
	pctMatch <- round( nMatch * 100 / NP, digits=2)
	pctNotMatch <- round( nNotMatch * 100 / NP, digits=2)
	matchTbl <- table( centerBase, calledBase)
	return( list( "PctMatch"=pctMatch, "PctNotMatch"=pctNotMatch, "MatchTable"=matchTbl))
}


`chromatogramToTable` <- function( chromoObj) {

	# turn all the pieces into what we need to make one table
	# trace matrix is trivial
	traceM <- chromoObj$TraceM
	ntr <- nrow( traceM)
	row.names <- rownames( traceM)

	# put the peak locations where they go
	peakVec <- peakCount <- peakConf <- rep.int( NA, ntr)
	peakLocs <- as.integer( chromoObj$PeakPosition)
	peakVec[ peakLocs] <- peakLocs
	peakCount[ peakLocs] <- 1:length(peakLocs)
	peakConf[ peakLocs] <- as.numeric( chromoObj$PeakConfidence)

	# show just the 'fwd' strand calls at the peaks
	dna <- chromoObj$DNA_Calls
	dnaV <- rep.int( "", ntr)
	dnaCalls <- strsplit( dna[1], split="")
	dnaV[ peakLocs] <- dnaCalls[[1]]

	# and show all the AA reading frames too
	aa <- chromoObj$AA_Calls
	aaM <- matrix( "", nrow=ntr, ncol=6)
	colnames(aaM) <- names(aa)
	aaCalls <- strsplit( aa, split="")
	read1.4 <- seq( 2, length(peakLocs), by=3)
	read2.5 <- seq( 3, length(peakLocs), by=3)
	read3.6 <- seq( 4, length(peakLocs), by=3)

	# makw sure the lengths are all the same
	if ( length(aaCalls[[1]]) > length( read1.4)) {
		length(aaCalls[[1]]) <- length(aaCalls[[4]]) <- length(read1.4)
	}
	if ( length(aaCalls[[1]]) < length( read1.4)) length(read1.4) <- length(aaCalls[[1]])
	if ( length(aaCalls[[2]]) > length( read2.5)) {
		length(aaCalls[[2]]) <- length(aaCalls[[5]]) <- length(read2.5)
	}
	if ( length(aaCalls[[2]]) < length( read2.5)) length(read2.5) <- length(aaCalls[[2]])
	if ( length(aaCalls[[3]]) > length( read3.6)) {
		length(aaCalls[[3]]) <- length(aaCalls[[6]]) <- length(read3.6)
	}
	if ( length(aaCalls[[3]]) < length( read3.6)) length(read3.6) <- length(aaCalls[[3]])

	# now stuff those AA calls in
	aaM[ peakLocs[read1.4], 1] <- aaCalls[[1]]
	aaM[ peakLocs[read1.4], 4] <- aaCalls[[4]]
	aaM[ peakLocs[read2.5], 2] <- aaCalls[[2]]
	aaM[ peakLocs[read2.5], 5] <- aaCalls[[5]]
	aaM[ peakLocs[read3.6], 3] <- aaCalls[[3]]
	aaM[ peakLocs[read3.6], 6] <- aaCalls[[6]]

	out <- data.frame( "Row"=row.names, traceM, "PeakPosition"=peakVec, "PeakCount"=peakCount, 
			"PeakConfidence"=peakConf, "DNA_Call"=dnaV, aaM, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	out
}


`chromatogramFromTable` <- function( tbl) {

	expectedColumns1 <- c( "Row", "A", "C", "G", "T", "PeakPosition", "PeakConfidence", "DNA_Call")
	expectedColumns2 <- paste( "AA_Frame", 1:6, sep="")
	if ( ncol(tbl) != 14 || ! all( colnames(tbl)[1:8] == expectedColumns1)) {
		cat( "\nGiven 'tbl' object does not have chromatogram column names..")
		stop()
	}

	traceIn <- as.matrix( tbl[ , 2:5])
	rownames(traceIn) <- as.integer( tbl$Row)
	traceIn[ is.na(traceIn)] <- 0

	peaks <- as.integer( tbl$PeakPosition)
	good <- which( !is.na( peaks))
	peaksV <- peaks[ good]
	peaksConf <- as.numeric( tbl$PeakConfidence[good])

	dna <- tbl$DNA_Call
	dnaV <- dna[good]
	dnaStr <- paste( dnaV, collapse="")
	seqData <- chromatogramSequences(dnaStr)

	# at this point, we don't know the file's name, so just leave a placeholder
	out <- list( "TraceM"=traceIn, "PeakPosition"=peaksV, "PeakConfidence"=peaksConf, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"="")
	out
}


`subsetChromatogram` <- function( chromoObj, seq=NULL, range=NULL, min.unit.score=NULL, verbose=FALSE) {

	# allow being given a filename
	if ( is.character(chromoObj) && file.exists( chromoObj[1])) {
		chromoObj <- loadChromatogram( chromoObj)
	}
	
	if ( ! is.null( seq)) {
		return( subsetChromatogramBySequence( chromoObj, seq=seq, min.unit.score=min.unit.score, verbose=verbose))
	} else if ( ! is.null( range)) {
		return( subsetChromatogramByRange( chromoObj, range=range, verbose=verbose))
	} else {
		cat( "\nMust specifiy a chromatogram subset by sequence or range.")
	}
	return( NULL)
}


`subsetChromatogramBySequence` <- function( chromoObj, seq, min.unit.score=NULL, subset.type=NULL, verbose=FALSE) {

	# get the data we need 
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	peakConf <- chromoObj$PeakConfidence
	baseCall <- names( peakPos)
	NT <- nrow( traceM)
	NB <- length( peakPos)
	halfPeak <- floor( median( diff( peakPos)) / 2)

	# allow using just a subset of the full length
	firstTracePoint <- 1
	lastTracePoint <- NT

	# we will find DNA or AA in any reading frame
	subSeq <- toupper( as.character( seq)[1])
	chromoDNA <- chromoObj$DNA_Calls
	chromoAA <- chromoObj$AA_Calls

	# look in both DNA and AA, (since it's 2 different scoring matrices, compensate a bit)
	require( Biostrings)
	data( BLOSUM62)
	DNA_MATRIX <- nucleotideSubstitutionMatrix()
	notDNA <- grepl( "Q|E|I|L|F|P|J|Z|X|\\*", subSeq)
	if (notDNA) {
		subSeq <- gsub( "?", "X", subSeq, fixed=T)
	} else {
		subSeq <- gsub( "?", "N", subSeq, fixed=T)
		subSeq <- gsub( "-", "N", subSeq, fixed=T)
	}
	
	# by default, we assume the subset sequence we are finding is smaller than the full chromatogram,
	# and we use this to find/extract the smaller construct, but this is user settable even if it breaks some
	# assumptions below...
	if ( is.null(subset.type)) subset.type <- "local-global"
	
	# find how well the given sequence matches in both worlds
	dnaScores <- if ( notDNA) -99999 else pairwiseAlignment( chromoDNA, subSeq, type=subset.type, substitutionMatrix=DNA_MATRIX, scoreOnly=T) * 3
	aaScores <- pairwiseAlignment( chromoAA, subSeq, type=subset.type, substitutionMatrix=BLOSUM62, scoreOnly=T)

	# get the exact limits on both DNA and AA, regardless of which one we were given
	if ( max(dnaScores) > max(aaScores)) {
		best <- which.max( dnaScores)
		bestScore <- dnaScores[ best]
		bestDNA <- chromoDNA[ best]
		pa <- pairwiseAlignment( bestDNA, subSeq, type="local-global", substitutionMatrix=DNA_MATRIX, scoreOnly=F)
		startDNA <- start( pattern( pa))
		stopDNA <- startDNA + width( pattern( pa)) - 1
		needChromoRevComp <- (best > 1)
		AAoffset <- 1
		if (verbose) cat( "\nBest Subset was DNA: ", best, bestScore, substr(bestDNA,startDNA,stopDNA))
	} else {
		best <- which.max( aaScores)
		bestScore <- aaScores[ best]
		bestAA <- chromoAA[ best]
		pa <- pairwiseAlignment( bestAA, subSeq, type="local-global", substitutionMatrix=BLOSUM62, scoreOnly=F)
		startAA <- start( pattern( pa))
		stopAA <- startAA + width( pattern( pa)) - 1
		stopDNA <- (stopAA*3)
		startDNA <- startAA*3 - 2
		needChromoRevComp <- (best > 3)
		AAoffset <- if (best < 4) best else best - 3
		if (verbose) cat( "\nBest Subset was AA:  ", best, bestScore, substr(bestAA,startAA,stopAA))
	}

	# see if we should return 'not found'
	unitScore <- bestScore / nchar( subSeq)
	if ( ! is.null( min.unit.score)) {
		if ( unitScore < min.unit.score) return( NULL)
	}

	# did we need to do a RevComp?
	if (needChromoRevComp) {
		tmpObj <- revCompChromatogram( chromoObj)
		traceM <- tmpObj$TraceM
		peakPos <- tmpObj$PeakPosition
		peakConf <- tmpObj$PeakConfidence
		baseCall <- names( peakPos)
		NT <- nrow( traceM)
		NB <- length( peakPos)
	}

	# now turn those DNA spots into traceM locations, accounting for any reading frame shift
	startDNA <- startDNA + AAoffset - 1
	stopDNA <- stopDNA + AAoffset - 1
	firstTracePoint <- max( 1, peakPos[ startDNA] - halfPeak)
	lastTracePoint <- min( NT, peakPos[ stopDNA] + halfPeak)

	# now trim all data to that subset
	xLimits <- c( firstTracePoint, lastTracePoint)
	traceOut <- traceM[ firstTracePoint:lastTracePoint, ]
	peaksOut <- peakPos[ peakPos > firstTracePoint & peakPos < lastTracePoint]
	confOut <- peakConf[ peakPos > firstTracePoint & peakPos < lastTracePoint]
	
	# make the peaks be in this local reference system now
	peaksOut <- peaksOut - firstTracePoint + 1

	# make the other fields the chromatogram has
	seqOut <- paste( names(peaksOut), collapse="")
	seqData <- chromatogramSequences( seqOut)

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "PeakConfidence"=confOut, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename,
				"Offset"=firstTracePoint, "BestAAframe"=AAoffset, "UnitScore"=unitScore)
	out
}


`subsetChromatogramByRange` <- function( chromoObj, range, verbose=FALSE) {

	# get the data we need 
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	peakConf <- chromoObj$PeakConfidence
	baseCall <- names( peakPos)
	NT <- nrow( traceM)
	NB <- length( peakPos)
	halfPeak <- floor( median( diff( peakPos)) / 2)

	# allow using just a subset of the full length
	firstTracePoint <- 1
	lastTracePoint <- NT

	peak1 <- max( 1, min( range, na.rm=T))
	peakN <- min( length(peakPos), max( range, na.rm=T))
	firstPeak <- peakPos[ peak1]
	lastPeak <- peakPos[ peakN]
	firstTracePoint <- max( 1, firstPeak - halfPeak)
	lastTracePoint <- min( NT, lastPeak + halfPeak)

	# now trim all data to that subset
	xLimits <- c( firstTracePoint, lastTracePoint)
	traceOut <- traceM[ firstTracePoint:lastTracePoint, ]
	peaksOut <- peakPos[ peakPos > firstTracePoint & peakPos < lastTracePoint]
	confOut <- peakConf[ peakPos > firstTracePoint & peakPos < lastTracePoint]
	
	# make the peaks be in this local reference system now
	peaksOut <- peaksOut - firstTracePoint + 1

	# make the other fields the chromatogram has
	seqOut <- paste( names(peaksOut), collapse="")
	seqData <- chromatogramSequences( seqOut)

	out <- list( "TraceM"=traceOut, "PeakPosition"=peaksOut, "PeakConfidence"=confOut, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename,
				"Offset"=firstTracePoint)
	out
}



`peakpickChromatogram` <- function( chromoObj) {

	# get the data we need 
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	peakConf <- chromoObj$PeakConfidence

	# ignore all the existing peak calls, and just call new ones directly from the trace matrix
	peakPos <- peakpickTraceM( traceM)
	
	# with this new set of peak calls, remake all the other parts
	baseCalls <- names( peakPos)
	dna.seq <- paste( baseCalls, collapse="")
		
	# we will always keep both strands, as both DNA and AA
	seqData <- chromatogramSequences( dna.seq)
	
	# use those peak top locations to estimate where the trace stops being useful
	avgPeakSeparation <- median( diff( peakPos))
	rightTail <- max( peakPos) + avgPeakSeparation
	if ( rightTail < nrow(traceM)) {
		traceM <- traceM[ 1:rightTail, ]
	}

	# there could also be rare cases where the trace matrix does not contain all the peaks.
	# this should never happen, but check and fix if it does
	NTR <- nrow( traceM)
	maxPeakLoc <- NTR - floor( avgPeakSeparation/2)
	if ( any( peakPos > maxPeakLoc)) {
		# throw away those peaks
		drops <- which( peakPos > maxPeakLoc)
		peakPos <- peakPos[ -drops]
		# this means we need to truncate the DNA calls too
		if ( nchar(dna.seq) > length(peakPos)) {
			dna.seq <- substr( dna.seq, 1, length(peakPos))
			seqData <- chromatogramSequences( dna.seq)
		}
	}

	colnames(traceM) <- c( "A", "C", "G", "T")
	rownames(traceM) <- 1:nrow(traceM)
	
	# lastly, calculate confidence scores for each peak
	peakConfidence <- calcChromatogramPeakConfidence( traceM, peakPos)

	out <- list( "TraceM"=traceM, "PeakPosition"=peakPos, "PeakConfidence"=peakConfidence, 
				"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename)
	out
}


`peakpickTraceM` <- function( traceM) {

	# given a matrix of A/C/G/T trace lines, find the peaks from first principles
	
	# start with total intensity in all channels
	N <- nrow(traceM)
	inten <- apply( traceM, 1, sum, na.rm=T)
	intenM1 <- c( 0, inten[1:(N-1)])
	intenP1 <- c( inten[2:N], 0)
	min.inten <- max( inten) * 0.01
	
	# find all extrema as simple high points above background
	isExtrema <- which( inten > min.inten & inten > intenM1 & inten >= intenP1 )
	extremaHts <- inten[ isExtrema]
	Nextrema <- length( isExtrema)
	if ( Nextrema < 4) {
		cat( "\nError: not enough peaks found.  Invalid chromatogram")
		return( NULL)
	}
	
	# we can get a guess of typical peak spacing
	avg.peak.separation <- median( diff( isExtrema))
	peak.tail <- round( avg.peak.separation * 0.6)
	
	# now visit every extrema in order by height, and mask out the region under that peak
	usedMask <- rep.int( FALSE, N)
	visitOrd <- order( extremaHts[ isExtrema], decreasing=T)
	peakLoc <- peakBase <- vector()
	Npeaks <- 0
	for ( i in 1:Nextrema) {
		thisPeak <- isExtrema[ visitOrd[i]]
		# not a real peak if already covered by someone taller
		if ( usedMask[thisPeak]) next
		
		# OK, this is a good peak to keep
		Npeaks <- Npeaks + 1
		peakLoc[Npeaks] <- thisPeak
		peakBase[Npeaks] <- colnames(traceM)[ which.max( traceM[ thisPeak, ])]
		
		# now mask its region out
		loEdge <- max( 1, thisPeak - peak.tail)
		hiEdge <- min( N, thisPeak + peak.tail)
		usedMask[ loEdge:hiEdge] <- TRUE
	}
	
	# done visiting the obvious peaks
	
	# perhaps search for dead spots...??
	
	# package up the results
	ord <- order( peakLoc)
	peaksOut <- peakLoc[ ord]
	names( peaksOut) <- peakBase[ord]
	return( peaksOut)
}


`plotChromatogram` <- function( chromoObj, label="", seq=NULL, range=NULL, 
				lwd=2, lty=1, cex=1, font=2, add=FALSE, forceYmax=NULL, 
				showAA=TRUE, showTraceRowNumbers=FALSE, showConfidence=FALSE,
				min.unit.score=NULL, xlim=NULL, shiftAA=0, min.intensity.plot=0.1,
				main.prefix="Chromatogram:  ", ...) {

	# allow being given a filename of a chromatogram
	if ( is.character(chromoObj) && file.exists( chromoObj[1])) {
		chromoObj <- loadChromatogram( chromoObj)
	}
	
	# given a chromatogram object, show the whole thing
	acgtBases <- c('A','C','G','T','N','-')
	acgtColors <- c('red','blue','orange','green','brown','black')

	neededNames <- c( "TraceM", "PeakPosition")
	if ( ! all( neededNames %in% names(chromoObj))) {
		cat( "Chromatogram object does not have needed fields..")
		return()
	}

	# were we given a request for a smaller region?
	if ( ! is.null( seq)) {
		chromoObj <- subsetChromatogramBySequence(chromoObj, seq=seq, min.unit.score=min.unit.score)
		# allow the subsequence request to return silently with no plot if not good sequence found
		if ( is.null( chromoObj)) return( NULL)
	} else if ( ! is.null( range)) {
		chromoObj <- subsetChromatogramByRange(chromoObj, range=range)
	}
	
	# get the data we need 
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	peakConf <- chromoObj$PeakConfidence
	baseCall <- names( peakPos)
	NT <- nrow( traceM)
	NB <- length( peakPos)

	# allow using just a subset of the full length
	firstTracePoint <- 1
	lastTracePoint <- NT
	xLimits <- c( 1, NT)
	# allow the caller to shift the plot along X when doing 2+ chromatotrams
	if ( ! is.null( xlim)) {
		xLimits[1] <- min( xlim)
		xLimits[2] <- max( xlim)
	}
	AAtoShow <- ""
	AAoffset <- 1
	showAllAA <- ( is.character(showAA) && showAA =="all")
	if ( showAllAA || (is.logical(showAA) && showAA)) AAtoShow <- chromoObj$AA_Calls[AAoffset]
	if ( "BestAAframe" %in% names(chromoObj)) {
		AAoffset <- as.integer( chromoObj$BestAAframe)
 		AAtoShow <- chromoObj$AA_Calls[AAoffset]
	}
	if ( is.numeric(showAA)) {
		AAtoShow <- chromoObj$AA_Calls[showAA]
		AAoffset <- showAA
	}

	mainText <- paste( main.prefix, label)
	x <- 1 : NT
	yLimits <- c( 0, max( traceM[ firstTracePoint:lastTracePoint, ], na.rm=T) * 1.1)
	if ( ! is.null(forceYmax)) yLimits[2] <- as.numeric( forceYmax[1])

	if ( ! add) plot( 1,1, type="n", main=mainText, xlim=xLimits, ylim=yLimits, ylab="Intensity", xlab=NA,
				xaxt="n", xaxs="i", las=2, ...)

	for ( j in 1:4) {
		y <- traceM[ ,j]
		if ( all( y <= min.intensity.plot)) next
		lines( x, y, col=acgtColors[j], lwd=lwd, lty=lty)
	}
	baseColor <- acgtColors[ match( baseCall, acgtBases)]
	# flag indel cleaning with extra color
	baseColor[ is.na(baseColor)] <- 'purple'

	for (k in 1:NB) axis( side=1, at=peakPos[k], label=baseCall[k], col.axis=baseColor[k], 
				col.ticks=baseColor[k], font=font, cex.axis=cex)

	if ( AAtoShow != "" && nchar(AAtoShow) >= 1) {
		aaCall <- strsplit( AAtoShow, split="")[[1]]
		NAA <- min( length(aaCall), round(NB/3))
		for (k in 1:NAA) {
			kk <- (k-1) * 3 + 1 + AAoffset
			if (kk > 0 & kk <= length(peakPos)) {
				axis( side=1, at=peakPos[kk]+shiftAA, label=aaCall[k], line=1, col.axis='black', col.ticks=NA, 
					tick=FALSE, font=2, lwd.ticks=0, cex.axis=cex*1.4)
			}
		}
	}
	if ( showAllAA) {
		frameToDo <- setdiff( 1:3, AAoffset)
		lineNow <- 2.0
		for (frame in frameToDo) {
			AAtoShow <- chromoObj$AA_Calls[frame]
			aaCall <- strsplit( AAtoShow, split="")[[1]]
			NAA <- min( length(aaCall), round(NB/3))
			if ( NAA > 1) for (k in 1:NAA) {
				kk <- (k-1) * 3 + 1 + frame
				if (kk > 0 & kk <= length(peakPos)) {
					axis( side=1, at=peakPos[kk]+shiftAA, label=aaCall[k], line=lineNow, col.axis='black', col.ticks=NA, 
							tick=FALSE, font=1, lwd.ticks=0, cex.axis=cex*1.0)
				}
			}
			lineNow <- lineNow + 0.75
		}
	}

	if ( showTraceRowNumbers) {
		traceRowValues <- as.numeric( rownames(traceM)[ firstTracePoint:lastTracePoint])
		myAtValues <- pretty( traceRowValues, n=9)
		myAts <- myAtValues - traceRowValues[1] + 1
		axis( side=3, at=myAts, labels=myAtValues, padj=1)
	}

	if ( showConfidence) {
		# put zero to a hundred on the right axis, and scale the confidence data
		atConf <- seq( 0, 100, by=20)
		yScale <- yLimits[2] / 100
		atY <- atConf * yScale
		axis( side=4, at=atY, labels=atConf, las=0, tck=-0.01, mgp=c(3,0.2,0))
		mtext( side=4, "Confidence", line=1)
		confX <- peakPos
		# raw confidence is 0..1, scale to 0..100
		confY <- peakConf * 100 * yScale
		lines( confX, confY, lty=3, lwd=1, col='black')
		# and the moving average too
		if ( length(confY) >= 15) {
			names(confY) <- confX
			moveAvg <- movingAverage( confY, window=11)
			lines( confX, moveAvg, lty=1, lwd=1, col='black')
		}
	}
	dev.flush()
}


`setChromatogramBestAAframe` <- function( chromoObj, referenceAAseq) {

	# given a expected protein sequence, set the internal object pointer to 
	# which AA frame is best match	
	require( Biostrings)
	data( BLOSUM62)
	mySeqs <- chromoObj$AA_Calls
	scores <- pairwiseAlignment( mySeqs, referenceAAseq, type="local", scoreOnly=T, substitutionMatrix=BLOSUM62)
	chromoObj$BestAAframe <- which.max( scores)
	return( chromoObj)
}



`plotMultipleChromatograms` <- function( chromoSet, label="", seq=NULL, showAA=TRUE, 
					showTraceRowNumbers=FALSE, showConfidence=FALSE, 
					min.unit.score=NULL, max.to.show=4, xlim=NULL, ...) {

	# set up to draw more than one
	nChromo <- length( chromoSet)
	nShow <- min( nChromo, max.to.show)
	par( mfrow=c( nShow, 1))
	par( mai=c(0.8, 0.8, 0.7, 0.4))
	
	# allow being given a list of X limits
	myXlimits <- xlim

	for ( i in 1:nChromo) {
		chromoObj <- chromoSet[[i]]
		chromoLabel <- paste( label, names(chromoSet)[i], sep="  ")
		if ( ! is.null( xlim)) {
			if ( is.list( xlim)) myXlimits <- xlim[[ min( i, length(xlim))]]
		}
		plotChromatogram( chromoObj, label=chromoLabel, seq=seq, showAA=showAA,
				showTraceRowNumbers=showTraceRowNumbers, showConfidence=showConfidence, 
				min.unit.score=min.unit.score, xlim=myXlimits, ...)
	}
}


`plotFullChromatogramAsPDF` <- function( chromoFile, label=basename(chromoFile), reference.seq=NULL, 
				cex=1, forceYmax=NULL, showAA="all", showTraceRowNumbers=TRUE, 
				showConfidence=TRUE, plot.path=dirname(chromoFile), ...) {

	# given a chromatogram object, show the whole thing
	chromoObj <- loadChromatogram( chromoFile)
	if ( is.null(chromoObj)) return(NULL)

	# PDF will be put where the raw data is
	pdfFile <- sub( ".ab1$", ".pdf", basename(chromoFile))
	pdfFile <- file.path( plot.path, pdfFile)

	# use the length of the raw data to set the device
	nTraceRows <-nrow( chromoObj$TraceM)
	plotWidth <- round( nTraceRows/95)
	plotWidth <- max( plotWidth, 6)
	plotHeight <- 6

	# since the raw data usually has big spikes, clip Y a bit
	if ( is.null( forceYmax)) {
		peakPts <- chromoObj$PeakPosition
		hts <- apply( chromoObj$TraceM[ peakPts, ], MARGIN=1, max, na.rm=T)
		forceYmax <- quantile( hts, 0.96)
	}
	
	# build a informative label
	label <- paste( label, "     ", "DNA_Len=", nchar(chromoObj$DNA_Calls[1]), "   AA_Len=", 
			nchar(chromoObj$AA_Calls[1]), "   Avg_Confidence=", 
			round(mean(chromoObj$PeakConfidence)*100,digits=1))

	pdf( pdfFile, width=plotWidth, height=plotHeight)

	plotChromatogram( chromoObj, label=label, seq=reference.seq, showAA=showAA,
			showTraceRowNumbers=showTraceRowNumbers, showConfidence=showConfidence, 
			forceYmax=forceYmax, ...)

	dev.off()
}


`cleanChromatogramDNA` <- function( chromoObj, referenceDNA, referenceAA=NULL, nSkip=1, verbose=TRUE) {

	# we find cases where Sanger sequencing is inserting duplicate bases, which throw off reading frame
	# try to find and correct, using a reference DNA and perhaps AA sequence.
	require( Biostrings)
	DNA_MATRIX <- nucleotideSubstitutionMatrix()

	# first efforts, quite subjective, currently only changing the DNA & AA fields
	#		not yet modifying the trace matrix or peak calls...
	chromoDNA <- chromoObj$DNA_Calls

	# set up to know what to return
	outObj <- chromoObj
	fixedNs <- FALSE
	fixedDNA <- FALSE
	isRevComp <- FALSE

	# step 1: see which strand is better match to reference
	dnaScores <- pairwiseAlignment( chromoDNA, referenceDNA, type="global-local", substitutionMatrix=DNA_MATRIX, scoreOnly=T)
	bestStrand <- which.max( dnaScores)
	if (bestStrand == 2) isRevComp <- TRUE

	# if it is the Reverse, then we need to RevComp everything before we proceed
	if (isRevComp) {
		outObj <- revCompChromatogram( outObj)
	}
	traceM <- outObj$TraceM
	peakPos <- outObj$PeakPosition
	peakConf <- outObj$PeakConfidence
	dna <- outObj$DNA_Calls[1]

	# step 2: try to re-call "N" bases manually, but ignore the edges
	hasNcall <- which( names(peakPos) == "N")
	NP <- length(peakPos)
	hasNcall <- setdiff( hasNcall, c( 1:nSkip, (NP-nSkip+1):NP))
	if ( NhasN <- length( hasNcall)) {
		hasNnew <- hasNstr <- vector( length=NhasN)
		for (ik in 1:NhasN) {
			k <- hasNcall[ ik]
			baseAns <- baseCallOnePeak( peakPos[k], traceM)
			hasNnew[ik] <- newBase <- names(baseAns)[1]
			hasNstr[ik] <- paste( names(baseAns), as.numeric(baseAns), sep=":", collapse="; ")
			if ( newBase != "N") {
				substr( dna, k, k) <- newBase
				names(peakPos)[k] <- newBase
				names(peakConf)[k] <- newBase
				fixedNs <- TRUE
			}
		}
		hasN.details <- data.frame( "Peak"=hasNcall, "OldBase"="N", "NewBase"=hasNnew, "Base.Percents"=hasNstr, 
					stringsAsFactors=F)
	} else {
		hasN.details <- NULL
	}

	# step 3: repeatedly see if we find a 1-2bp gap
	ans <- DNAtoCleanChromatogramDNA( dna, referenceDNA=referenceDNA, referenceAA=referenceAA, verbose=verbose)
	dnaNow <- ans$DNA
	details <- ans$CleaningDetails
	outPeaks <- peakPos
	outConf <- peakConf
	if ( ! is.null( details)) {
		# the details show where bases got added/removed...
		# we can push these changes into the chromatogram
		for ( k in 1:nrow(details)) {
			changeLoc <- details$Location[k]
			changeType <- details$CleaningType[k]
			baseWas <- details$OriginalBase[k]
			baseNew <- details$CleanedBase[k]
			chromoBase <- names( peakPos)[ changeLoc]
			# verify they say what they should
			if ( chromoBase != baseWas) cat( "\nDebug assertion failed: cleaning base error: ", chromoBase, baseWas)
			#names(peakPos)[changeLoc] <- baseNew
			# instead of just noting this in the name field, incorporate it
			if (changeType == "ExtraBase") {
				# the raw chromatogram seems to have an extra base that messes up reading frame
				# force remove it
				outPeaks[ changeLoc] <- NA
				outConf[ changeLoc] <- NA
			}
			if (changeType == "MissingBase") {
				# the raw chromatogram seems to be missing a base that messes up reading frame
				# force add an extra, by shifting this one back a hair and adding another ahead a hair
				oldCenter <- outPeaks[ changeLoc]
				oldConf <- outConf[ changeLoc]
				outPeaks[ changeLoc] <- oldCenter - 1
				newPeak <- oldCenter + 2
				newConf <- oldConf
				names(newPeak) <- names(newConf) <- substr( baseNew,2,2)
				outPeaks <- c( outPeaks, newPeak)
				outConf <- c( outConf, newConf)
			}
		}
		# after all the changes, update the peak details
		isNA <- which( is.na( outPeaks))
		if (length(isNA)) {
			outPeaks <- outPeaks[ -isNA]
			outConf <- outConf[ -isNA]
		}
		ord <- order( outPeaks)
		outPeaks <- outPeaks[ord]
		outConf <- outConf[ord]
		fixedDNA <- TRUE
	}

	# do we need to repackage?
	if ( !( fixedNs || fixedDNA)) return( chromoObj)

	# OK, the cleaned DNA string becomes the truth now
	dnaOut <- paste( names( outPeaks), collapse="")
	seqData <- chromatogramSequences( dnaOut)
	out <- list( "TraceM"=outObj$TraceM, "PeakPosition"=outPeaks, "PeakConfidence"=outConf,
			"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=outObj$Filename)
	# and undo the Rev Comp if we did it
	if ( isRevComp) {
		out <- revCompChromatogram( out)
		# also, the cleaning details need a bit of rev comp cleaning...
		if ( ! is.null( details)) {
			for ( k in 1:nrow(details)) {
				details$Location[k] <- length(outPeaks) - details$Location[k] + 1
				details$OriginalBase[k] <- myReverseComplement( details$OriginalBase[k])
				details$CleanedBase[k] <- myReverseComplement( details$CleanedBase[k])
			}
		}
	}

	# append the cleaning details
	if ( !is.null( details) || !is.null( hasN.details)) {
		out$CleaningDetails <- list( "Fixed.N.Calls"=hasN.details, "Fixed.IndelGaps"=details)
	}
	out
}


`DNAtoCleanChromatogramDNA` <- function( dna, referenceDNA, referenceAA=NULL, verbose=T) {

	# lower level tool to do the cleaning...
	require( Biostrings)
	DNA_SUBM <- nucleotideSubstitutionMatrix()
	data( BLOSUM62)
	AA_SUBM <- BLOSUM62

	ORDER <- base::order
	PASTE <- base::paste
	SETDIFF <- base::setdiff
	STRSPLIT <- base::strsplit
	SUBSTR <- base::substr
	WHICH <- base::which
	WHICH.MAX <- base::which.max
	WHICH.MIN <- base::which.min

	# DNA from chromatograms can have repeats or gaps that destroy reading frame.
	# try to find and fix
	dnaNow <- dna[1]

	# accumulate a set of locations already visited
	maskBases <- vector()

	# and try to build a table of details if the caller wants to know exactly what changed
	fixLoc <- fixType <- fixWas <- fixNew <- vector()
	nFix <- 0

	repeat {
		pa <- pairwiseAlignment( dnaNow, referenceDNA, type="local", scoreOnly=F, substitutionMatrix=DNA_SUBM)
		dnaStart <- start( pattern( pa))
		dnaStr <- as.character( alignedPattern( pa))
		refStart <- start( subject( pa))
		refStr <- as.character( alignedSubject( pa))
		scoreNow <- score(pa)
		# if we don't have a good score, don't even try
		if ( scoreNow < nchar(dnaNow)*0.1) break

		# do we see any context that looks like a 1bp gap?
		# use a fixed size flank of 9bp
		N_FLANK <- 9
		#dnaGapAns <- as.numeric( gregexpr( "[ACGTN]{9}\\-[ACGTN]{9}", dnaStr)[[1]])
		#refGapAns <- as.numeric( gregexpr( "[ACGTN]{9}\\-[ACGTN]{9}", refStr)[[1]])
		# let's not try to fix anything too close to 'N' calls...
		dnaGapAns <- as.numeric( gregexpr( "[ACGT]{9}\\-[ACGT]{9}", dnaStr)[[1]])
		refGapAns <- as.numeric( gregexpr( "[ACGT]{9}\\-[ACGT]{9}", refStr)[[1]])

		# don't revisit any sites already done
		if ( length( maskBases)) {
			dnaGapAns <- SETDIFF( dnaGapAns, maskBases)
			if( ! length( dnaGapAns)) dnaGapAns <- -1
			refGapAns <- SETDIFF( refGapAns, maskBases)
			if( ! length( refGapAns)) refGapAns <- -1
		}

		# in the case there is more than one, pick the one nearest the middle
		if ( length(dnaGapAns) > 1) {
			midPt <- nchar( dnaStr) / 2
			dx <- abs( dnaGapAns + N_FLANK - midPt)
			best <- WHICH.MIN( dx)
			dnaGapAns <- dnaGapAns[ best]
		}
		if ( length(refGapAns) > 1) {
			midPt <- nchar( refStr) / 2
			dx <- abs( refGapAns + N_FLANK - midPt)
			best <- WHICH.MIN( dx)
			refGapAns <- refGapAns[ best]
		}
		if ( dnaGapAns < 0 && refGapAns < 0) break

		# which one? and does removing it make things better?
		if ( refGapAns > 0) {
			# gap in reference means extra base in chromatogram
			gapLocation <- refGapAns + N_FLANK
			nOtherGaps <- sum( gregexpr( "-", SUBSTR( dnaStr, 1, refGapAns), fixed=T)[[1]] > 0)
			dnaLocation <- dnaStart + gapLocation - 1 - nOtherGaps
			newLeft <- SUBSTR( dnaNow, 1, dnaLocation-1)
			newRight <- SUBSTR( dnaNow, dnaLocation+1, nchar(dnaNow))
			dnaNew <- PASTE( newLeft, newRight, sep="")

			# is this better?
			if ( is.null( referenceAA)) {
				scoreNew <- pairwiseAlignment( dnaNew, referenceDNA, type="local", scoreOnly=T, substitutionMatrix=DNA_SUBM)
			} else {
				aaNow <- DNAtoBestPeptide(dnaNow, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				aaNew <- DNAtoBestPeptide(dnaNew, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				scoreNow <- pairwiseAlignment( aaNow, referenceAA, type="local", scoreOnly=T, substitutionMatrix=AA_SUBM)
				scoreNew <- pairwiseAlignment( aaNew, referenceAA, type="local", scoreOnly=T, substitutionMatrix=AA_SUBM)
			}
			extraBase <- SUBSTR( dnaNow, dnaLocation, dnaLocation)
			oldContext <- SUBSTR( dnaNow, dnaLocation-1, dnaLocation+1)
			newContext <- PASTE( SUBSTR(oldContext,1,1), SUBSTR(oldContext,3,3), sep="")
			pctBetter <- round( (scoreNew-scoreNow) * 100 / scoreNow, digits=2)
			if (verbose) {
				cat( "\n--------------------------")
				cat( "\nCleaning Chromatogram DNA:  found extra base in Chromatogram.")
				cat( "\nDNA Context:", SUBSTR(dnaStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nRef Context:", SUBSTR(refStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nExtra Base: ", extraBase,"   Was: ", oldContext, "   Now: ", newContext)
			}
			if ( scoreNew > scoreNow) {
				if (verbose) cat( "\nScore:    Was: ", scoreNow, "   Now: ", scoreNew, "   PctImprove: ", pctBetter, "%")
				dnaNow <- dnaNew
				nFix <- nFix + 1
				fixLoc[nFix] <- dnaLocation
				fixWas[nFix] <- extraBase
				fixNew[nFix] <- ""
				fixType[nFix] <- "ExtraBase"
			} else {
				if (verbose) cat( "\nBut Score no better!    Was: ", scoreNow, "   Now: ", scoreNew)
			}
			if (verbose) cat( "\n--------------------------\n")
			maskBases <- sort( unique.default( c( maskBases, (refGapAns - N_FLANK %/% 2):(refGapAns + N_FLANK %/% 2))))
			next
		}
		if ( dnaGapAns > 0) {
			# gap in chromatogram means extra base in reference 
			gapLocation <- dnaGapAns + N_FLANK
			nOtherGaps <- sum( gregexpr( "-", SUBSTR( dnaStr, 1, dnaGapAns), fixed=T)[[1]] > 0)
			extraRefBase <- SUBSTR( refStr, gapLocation, gapLocation)
			dnaLocation <- dnaStart + gapLocation - 1 - nOtherGaps
			newLeft <- SUBSTR( dnaNow, 1, dnaLocation-1)
			# since the gap is in my string, don't step past it
			#newRight <- SUBSTR( dnaNow, dnaLocation+1, nchar(dnaNow))
			newRight <- SUBSTR( dnaNow, dnaLocation, nchar(dnaNow))
			dnaNew <- PASTE( newLeft, extraRefBase, newRight, sep="")

			# is this better?
			if ( is.null( referenceAA)) {
				scoreNew <- pairwiseAlignment( dnaNew, referenceDNA, type="local", scoreOnly=T, substitutionMatrix=DNA_SUBM)
			} else {
				aaNow <- DNAtoBestPeptide(dnaNow, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				aaNew <- DNAtoBestPeptide(dnaNew, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				scoreNow <- pairwiseAlignment( aaNow, referenceAA, type="local", scoreOnly=T, substitutionMatrix=AA_SUBM)
				scoreNew <- pairwiseAlignment( aaNew, referenceAA, type="local", scoreOnly=T, substitutionMatrix=AA_SUBM)
			}
			oldContext <- SUBSTR( dnaNow, dnaLocation-1, dnaLocation)
			newContext <- PASTE( SUBSTR(oldContext,1,1), extraRefBase, SUBSTR(oldContext,2,2), sep="")
			pctBetter <- round( (scoreNew-scoreNow) * 100 / scoreNow, digits=2)
			if (verbose) {
				cat( "\n--------------------------")
				cat( "\nCleaning Chromatogram DNA:  found missing base in Chromatogram.")
				cat( "\nDNA Context:", SUBSTR(dnaStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nRef Context:", SUBSTR(refStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nMissing Base: ", extraRefBase,"   Was: ", oldContext, "   Now: ", newContext)
			}
			if ( scoreNew > scoreNow) {
				if (verbose) cat( "\nScore:    Was: ", scoreNow, "   Now: ", scoreNew, "   PctImprove: ", pctBetter, "%")
				dnaNow <- dnaNew
				nFix <- nFix + 1
				# let's call it as being appended to previous instead of prepended to the current base
				fixLoc[nFix] <- dnaLocation - 1
				fixWas[nFix] <- SUBSTR(oldContext,1,1)
				fixNew[nFix] <- SUBSTR(newContext,1,2)
				fixType[nFix] <- "MissingBase"
			} else {
				if (verbose) cat( "\nBut Score no better!    Was: ", scoreNow, "   Now: ", scoreNew)
			}
			if (verbose) cat( "\n--------------------------\n")
			maskBases <- sort( unique.default( c( maskBases, (dnaGapAns - N_FLANK %/% 2):(dnaGapAns + N_FLANK %/% 2))))
			next
		}
	}
	
	dnaOut <- dnaNow
	detailsOut <- NULL
	if (nFix) {
		detailsOut <- data.frame( "Location"=fixLoc, "OriginalBase"=fixWas, "CleanedBase"=fixNew,
					"CleaningType"=fixType, stringsAsFactors=F)
		# note that every cleaning modification we did throws off the locations w.r.t. the original DNA
		# we want the 'Location' data to be in caller's original units
		if ( nFix > 1) {
			ord <- ORDER( detailsOut$Location)
			detailsOut <- detailsOut[ ord, ]
			rownames(detailsOut) <- 1:nFix
			curDelta <- nchar(detailsOut$OriginalBase[1]) - nchar(detailsOut$CleanedBase[1])
			for ( i in 2:nFix) {
				detailsOut$Location[i] <- detailsOut$Location[i] + curDelta
				curDelta <- curDelta + nchar(detailsOut$OriginalBase[i]) - nchar(detailsOut$CleanedBase[i])
			}
		}
	}

	return( list( "DNA"=dnaOut, "CleaningDetails"=detailsOut))
}


# calculate the confidence of all peak calls
# by looking at the proportion of intensity that the top base explains
# returns a value in 0 to 1 range
`calcChromatogramPeakConfidence` <- function( traceM, peaks, windowSize=11) {

	# deduce how far apart peaks are, to estimate how many points around the center top to include
	avgPeakSep <- median( diff( peaks))
	half.width <- max( 1, floor( avgPeakSep / 3))

	# given how many peaks to look at, set up storage for our answers
	NTR <- nrow( traceM)
	NP <- length( peaks)
	peakDepth <- peakCount <- rep.int( NA, NP)
	
	# visit each peak, and look at that fraction of the trace matrix
	for (i in 1:NP) {
		myPeak <- peaks[i]
		leftPt <- max( 1, myPeak - half.width)
		rightPt <- min( NTR, myPeak + half.width)
		smlTM <- traceM[ leftPt:rightPt, ]

		# see how deep and how consistent
		cm <- apply( smlTM, MARGIN=2, sum, na.rm=T)
		peakDepth[i] <- sum( cm)

		# the calls are not "always" the biggext value, although they should be.
		# instead use the base call for the peak to get the one correct column
		myBase <- names(peaks)[i]
		peakCount[i] <- cm[ match(myBase, names(cm))]
		# 'N' calls would break this, so catch
		# but since 'N' is a terrible base call, reflect that in the peak's count
		if ( is.na( peakCount[i])) peakCount[i] <- ( max(cm) * 0.25)   # max( cm)
	}
	
	# prevent divide by zero
	peakDepth[ peakDepth < 1] <- 1
	peakPct <- peakCount / peakDepth

	# certainty is a function of Counts depth
	certainty <- 1.0 - (2 ^ -peakDepth)

	# we want/need to add another metric:   when a peak's height is way out of line with it's neighbors,
	# that is a sign of local trouble.  Use a smoothing window to assess what's expected, and punish if too
	# far outside that range.   We do need a minimum size to calc this...
	if ( NP > windowSize * 2) {
		tempCount <- peakCount
		names(tempCount) <- 1:NP
		smoothCount <- movingAverage( tempCount, window=windowSize)
		# turn the actual peak size into a ratio vs the average
		smoothCount[ smoothCount < 1] <- 1
		peakRatio <- peakCount / smoothCount
		# we expect as 'normal' peaks that are within 2x of average, like from mixed populations
		peakRatio[ peakRatio >= 0.5 & peakRatio <= 2.0] <- 1.0
		# for all other's, thurn this net ratio into a penalty on zero to one
		isBig <- which( peakRatio > 1.0)
		if ( length(isBig)) peakRatio[isBig] <- ( 1.0 / peakRatio[isBig])
		# lastly turn those zero to 0.5 values scaled up to one
		isSmall <- which( peakRatio <= 0.5)
		if ( length(isSmall)) peakRatio[isSmall] <- peakRatio[isSmall] * 2
	} else {
		peakRatio <- rep.int( 1, NP)
	}

	# confidence is the product of the dominant percentage times the certainty, 
	# with the new peak height ratio term too

	# always on the interval 0..1
	conf <- peakPct * certainty * peakRatio
	names(conf) <- names( peaks)
	conf
}


`selectBestChromatogramSequence` <- function( chromoObj, reference, type=c("DNA","AA")) {

	# given one chromatogram and reference sequence,
	# find the one chromatogram sequence (DNA or AA) that best matches that reference
	type <- match.arg( type)
	require( Biostrings)
	if ( type == "DNA") {
		subM <- nucleotideSubstitutionMatrix()
		mySeqs <- chromoObj$DNA_Calls
	} else {
		data( BLOSUM62)
		subM <- BLOSUM62
		mySeqs <- chromoObj$AA_Calls
	}
	
	# do the similarity scoring test
	paScores <- pairwiseAlignment( mySeqs, reference, type="local", scoreOnly=T, substitutionMatrix=subM)
	best <- which.max( paScores)
	
	# gather what we send back
	bestSeq <- mySeqs[ best]
	bestScore <- paScores[ best]
	bestLen <- nchar(bestSeq)
	unitScore <- round( bestScore / bestLen, digits=3)
	
	# also see exactly where this best sequence lands in the reference
	pa2 <- pairwiseAlignment( bestSeq, reference, type="local", scoreOnly=F, substitutionMatrix=subM)
	pattStart <- start( pattern( pa2))
	refStart <- start( subject( pa2))
	refStop <- width( subject( pa2)) + refStart - 1
	if ( pattStart > 1) {
		refStart <- refStart - pattStart + 1
		refStop <- refStop - pattStart + 1
	}
	
	out <- list( "SequenceName"=names(bestSeq), "ReferenceName"=names(reference), "Score"=bestScore, 
			"UnitScore"=unitScore, "RefStart"=refStart, "RefStop"=refStop, "Sequence"=bestSeq)

	# lastly, put the type into some of the column naems
	names(out)[3:7] <- paste( type, names(out)[3:7], sep=".")
	out
}


`getChromatogramSequenceConfidence` <- function( chromoObj, seq) {

	# given one chromatogram and one of it's DNA or AA sequences,
	# return the vector of confidence values for the elements of that seq string
	peakConf <- chromoObj$PeakConfidence
	
	# the forward DNA is trivial, as that is exact the confidence that is stored
	if ( seq == chromoObj$DNA_Calls[1]) return( peakConf)
	
	# the reverse DNA is almost as easy
	if ( seq == chromoObj$DNA_Calls[2]) {
		out <- rev( peakConf)
		names(out) <- strsplit(chromoObj$DNA_Calls[2], split="")[[1]]
		return( out)
	}

	# for the AA, we need to do merge/mean of the 3 bases per codon
	who <- match( seq, chromoObj$AA_Calls, nomatch=0)
	if  ( who == 0) return(NULL)
	
	# build the pair of indices that define the range for each codon
	N <- length( peakConf)
	if ( who %in% 1:3) {
		myFrom <- seq( who, N, by=3)
	} else {
		peakConf <- rev( peakConf)
		myFrom <- seq( (who-3), N, by=3)
	}
	myTo <- myFrom + 2
	
	# now we can step along each to calc the confidence for each aa
	aaConf <- mapply( myFrom, myTo, FUN=function(x,y) mean( peakConf[ x:y]))
	names( aaConf) <- strsplit( chromoObj$AA_Calls[who], split="")[[1]]
	return( aaConf)
}



`cropChromatogramLowSignalTail` <- function( chromoObj, min.signal.percent=20, windowSize=11, verbose=TRUE) {

	# grab the raw data that we start from
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	peakConf <- chromoObj$PeakConfidence
	NP <- length( peakPos)
	if ( NP <= windowSize) {
		if (verbose) cat( "  Crop Warning:  sequence too short, cropped to 'empty'.")
		return( NULL)
	}

	# get the intensity heights of the peaks
	peakHts <- apply( traceM[ peakPos, ], MARGIN=1, FUN=max)
	
	# apply a smoothing window to the peak heights
	smoothHts <- peakHts
	if ( NP >= windowSize * 2) {
		names(peakHts) <- peakPos
		smoothHts <- movingAverage( peakHts, window=windowSize)
	}

	# the cutoff may be given on 0 to one, or zero to 100
	# but the raw data is always zero to one
	if ( min.signal.percent > 1) min.signal.percent <- min.signal.percent / 100

	# make sure the front center portion -- which should be high height -- is
	centerPts <- round( c( NP*0.2, NP*0.5))
	meanHt <- mean( smoothHts[ centerPts[1] : centerPts[2]], na.rm=T)
	cropHtCutoff <- round( meanHt * min.signal.percent)

	# keep all the front, but
	# find the first place in the back half that is below the minimum height
	badRight <- which( smoothHts[ centerPts[2] : NP] < cropHtCutoff)
	keepLeft <- 1
	if ( length( badRight)) {
		keepRight <- min( badRight) + centerPts[2] - 1
	} else {
		keepRight <- NP
	}
	keep <- keepLeft : keepRight

	# grap those parts for each piece
	outPeaks <- peakPos[ keep]
	outConfs <- peakConf[ keep]
	NP2 <- length( outPeaks)
	halfPeak <- floor( median( diff( peakPos)) / 2)
	firstTrace <- 1
	lastTrace <- min( outPeaks[NP2] + halfPeak, nrow(traceM))
	outTrace <- traceM[ firstTrace:lastTrace, ]

	# now reset all the peak pointers
	rownames(outTrace) <- 1:nrow(outTrace)
	outPeaks <- outPeaks - firstTrace + 1

	# OK, the cleaned DNA string becomes the truth now
	dnaOut <- paste( names( outPeaks), collapse="")
	seqData <- chromatogramSequences( dnaOut)
	out <- list( "TraceM"=outTrace, "PeakPosition"=outPeaks, "PeakConfidence"=outConfs,
			"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename)
	out
}


`cropChromatogramByConfidence` <- function( chromoObj, min.confidence=60, windowSize=11, verbose=TRUE) {

	# grab the raw data that we start from
	traceM <- chromoObj$TraceM
	peakPos <- chromoObj$PeakPosition
	peakConf <- chromoObj$PeakConfidence
	NP <- length( peakConf)
	if ( NP <= windowSize) {
		if (verbose) cat( "  Crop Warning:  sequence too short, cropped to 'empty'.")
		return( NULL)
	}

	# apply a smoothing window to the raw confidence
	smoothConf <- confV <- as.numeric( peakConf)
	if ( length(confV) >= windowSize * 2) {
		names(confV) <- peakPos
		smoothConf <- movingAverage( confV, window=windowSize)
	}

	# the cutoff may be given on 0 to one, or zero to 100
	# but the raw data is always zero to one
	if ( min.confidence > 1) min.confidence <- min.confidence / 100

	# make sure the center -- which should be high confidence -- is
	centerPts <- round( c( NP*0.45, NP*0.55))
	if ( diff( centerPts) < 30) centerPts <- round( c( NP*0.4, NP*0.6))
	if ( diff( centerPts) < 30) centerPts <- round( c( NP*0.35, NP*0.65))
	middleConf <- mean( smoothConf[ centerPts[1] : centerPts[2]], na.rm=T)
	if ( middleConf <= min.confidence) {
		if (verbose) cat( "  Crop Warning:  center confidence too low, cropped to 'empty'.")
		return( NULL)
	}

	# find the last place in the front half that is below the confidence cutoff
	# and the first place in the back half
	badLeft <- which( smoothConf[ 1 : centerPts[1]] < min.confidence)
	badRight <- which( smoothConf[ centerPts[2] : NP] < min.confidence)
	if ( length( badLeft)) {
		keepLeft <- max( badLeft)
	} else {
		keepLeft <- 1
	}
	if ( length( badRight)) {
		keepRight <- min( badRight) + centerPts[2] - 1
	} else {
		keepRight <- NP
	}
	keep <- keepLeft : keepRight

	# grap those parts for each piece
	outPeaks <- peakPos[ keep]
	outConfs <- peakConf[ keep]
	NP2 <- length( outPeaks)
	halfPeak <- floor( median( diff( peakPos)) / 2)
	firstTrace <- max( outPeaks[1] - halfPeak, 1)
	lastTrace <- min( outPeaks[NP2] + halfPeak, nrow(traceM))
	outTrace <- traceM[ firstTrace:lastTrace, ]

	# now reset all the peak pointers
	rownames(outTrace) <- 1:nrow(outTrace)
	outPeaks <- outPeaks - firstTrace + 1

	# OK, the cleaned DNA string becomes the truth now
	dnaOut <- paste( names( outPeaks), collapse="")
	seqData <- chromatogramSequences( dnaOut)
	out <- list( "TraceM"=outTrace, "PeakPosition"=outPeaks, "PeakConfidence"=outConfs,
			"DNA_Calls"=seqData$DNA, "AA_Calls"=seqData$AA, "Filename"=chromoObj$Filename)
	out
}

