# chromatogramAssembly.R -- tools to build a single consensus assembly from sets of chromatogram files.  
#				Intended for cases where multiple Sanger sequencing primers are used
# 				to try to reconstruct one larger piece of gene


# top level tool designed to do the entire process of turning many chromatograms 
# covering one gene region into a single best consensus assembly

`chromatogramAssembly` <- function( chromoFiles, referenceFastaFile, doCleaning=TRUE, 
				cdsBases=NULL, verbose=TRUE) {

	# For conserved genes, 'Cleaning' can catch fix chromatogram errors.
	# But not advised for highly variant genes where the isolate likely has many indels

	# Step 1:  load the set of chromatograms into memory
	chromoSet <- loadMultipleChromatograms( chromoFiles, curated=T)
	nChromo <- length( chromoSet)
	if ( is.null( chromoSet)) return( NULL)
	
	# Step 2: load the reference DNA sequence.  We could be given 2+ different references, 
	# so we may need to find the one best reference for this isolate sample
	refFA <- loadFasta( referenceFastaFile, verbose=FALSE)
	refAns <- selectBestChromatogramReference( chromoSet, refFA)
	refDNA <- refAns$ReferenceSequence
	refName <- names(refDNA)[1]
	refAA <- translateChromatogramSequence( refDNA, mode="string", badAA="")
	refStats <- refAns$ScoreDetails
	if (verbose && !is.null( refStats)) {
		cat( "\nSelection Scoring of Best Reference:\n")
		print( refStats)
	}

	# Step 3:  now with the best reference DNA known, we may want to clean the chromatograms for 
	# typical sequence errors and artifacts
	if ( doCleaning) {
		for ( ichr in 1:nChromo) {
			if ( verbose) cat( "\nCleaning: ", names(chromoSet)[ichr])
			chromoObj <- chromoSet[[ichr]]
			chromoObj <- cleanChromatogramDNA( chromoObj, referenceDNA=refDNA, verbose=verbose)
			if ( "CleaningDetails" %in% names( chromoObj)) {
				writeCuratedChromatogram( chromoObj)
				chromoSet[[ichr]] <- chromoObj
			}
		}
	}
	
	# Step 4:  now extract the one best (correct strand) sequence for each chromatogram and it's 
	# confidence scores. we will build a data frame of details to pass back, 
	# and some info to help construct the matrix of base calls
	seqDF <- data.frame()
	seqInfo <- vector()
	confInfo <- list()
	for ( ichr in 1:nChromo) {
		chromoObj <- chromoSet[[ichr]]
		chromoName <- names(chromoSet)[ichr]
		smlAns <- selectBestChromatogramSequence( chromoObj, reference=refDNA, type="DNA")
		seqDF <- rbind( seqDF, data.frame( "Chromatogram"=chromoName, smlAns, stringsAsFactors=F))
		seqInfo[ichr] <- mySeq <- smlAns$DNA.Sequence
		myConf <- getChromatogramSequenceConfidence( chromoObj, seq=mySeq)
		confInfo[[ichr]] <- myConf
	}
	rownames(seqDF) <- 1:nrow(seqDF)
	seqDF$DNA.AvgConfidence <- round( sapply( confInfo, mean, na.rm=T) * 100)
	
	# Step 5:  use all these sequences and confidence scores to build one large matrix of base calls at all locations
	names(seqInfo) <- names(confInfo) <- seqDF$Chromatogram
	matrixAns <- constructChromatogramBaseMatrix( seqInfo, confInfo, refDNA)
	baseM <- matrixAns$BaseMatrix
	confM <- matrixAns$ConfidenceMatrix

	# Step 6:  now, use the confidence scores as the voting criteria, to extact the one final consensus
	consensusSeqAns <- extractChromatogramConsensusSequence( baseM, confM)
	finalDNAseq <- consensusSeqAns$sequence
	finalDNAstr <- paste( finalDNAseq, collapse="")
	finalDNAconf <- consensusSeqAns$confidence

	# Step 7:  we can now turn the DNA back to AA, for both the final consensus and all the fragments
	finalAA <- translateChromatogramSequence( finalDNAseq, mode="vector", badAA="", cdsBases=cdsBases)
	finalAAconf <- translateChromatogramConfidence( finalDNAconf, mode="vector", badAA="", cdsBases=cdsBases)
	names(finalAAconf) <- strsplit( finalAA, split="")[[1]]
	fragAA <- translateChromatogramSequence( baseM, mode="matrix", badAA=" ", cdsBases=cdsBases)
	fragAAconf <- translateChromatogramConfidence( confM, mode="matrix", badAA="", cdsBases=cdsBases)


	# we can make the AA scores and locations too, to complete the info
	trimmedAA <- gsub( " ", "", fragAA)
	data( BLOSUM62)
	pa <- pairwiseAlignment( trimmedAA, refAA, type="local", scoreOnly=F, substitutionMatrix=BLOSUM62)
	aaScores <- score(pa)
	aaStarts <- start( pattern( pa))
	aaStops <- width( pattern( pa)) + aaStarts - 1
	seqDF$AA.Score <- aaScores
	seqDF$AA.UnitScore <- round( aaScores / nchar( trimmedAA), digits=3)
	seqDF$AA.RefStart <- aaStarts
	seqDF$AA.RefStop <- aaStops
	seqDF$AA.AvgConfidence <- round( sapply( fragAAconf, mean, na.rm=T))
	seqDF$AA.Sequence <- fragAA
	# lastly, move the DNA seqeunce to the very end
	moveColumn <- which( colnames(seqDF) == "DNA.Sequence")
	detailsDF <- cbind( seqDF[ , -moveColumn], "DNA.Sequence"=seqDF$DNA.Sequence, stringsAsFactors=F)

	out <- list( "ReferenceName"=refName, "ReferenceAA"=refAA, "ReferenceDNA"=as.character(refDNA), "ConsensusAA"=finalAA,
			"ConfidenceAA"=finalAAconf, "ConsensusDNA"=finalDNAstr, "ConfidenceDNA"=finalDNAconf, 
			"FragmentDetails"=detailsDF, "BaseMatrix"=baseM, "ConfidenceMatrix"=confM)
	return( out)
}


`selectBestChromatogramReference` <- function( chromoSet, refFA) {

	# given a list of chromatograms and a FASTA object of DNA references
	# find the one best reference for this entire set of chromatograms
	nRef <- length( refFA$desc)
	
	# if we only have one, nothing to do
	if ( nRef < 2) {
		outSeq <- refFA$seq[1]
		names(outSeq) <- refFA$desc[1]
		out <- list( "ReferenceSequence"=outSeq, "ScoreDetails"=NULL)
		return( out)
	}
	
	# we expect a list of chromatograms, but catch the case if just one
	if ( all( c( "TraceM", "PeakPosition", "DNA_Calls", "AA_Calls") %in% names(chromoSet))) {
		tmp <- list()
		tmp[[1]] <- chromoSet
		chromoSet <- tmp
	}
	nChromo <- length( chromoSet)
	
	# if there are many chromatograms, only do a random subsample of them
	toDo <- if (nChromo <= 10) 1:nChromo else sample( nChromo, size=10)
	
	# set up storage to evaluate how well they match
	scoreM <- matrix( 0, nrow=length(toDo), ncol=nRef)
	colnames(scoreM) <- refFA$desc
	require( Biostrings)
	subM <- nucleotideSubstitutionMatrix()
	
	# do the similarity scoring over all pairs, looking at both Fwd and RevComp at the same time
	for (ichr in 1:length(toDo)) {
		ch <- chromoSet[[ichr]]
		dnaSet <- ch$DNA_Calls
		for ( iref in 1:nRef) {
			refSeq <- refFA$seq[ iref]
			paScores <- pairwiseAlignment( dnaSet, refSeq, type="local", scoreOnly=T, substitutionMatrix=subM)
			scoreM[ ichr, iref] <- max( paScores)
		}
	}
	
	# average over all the chromatograms tested, to give which reference is the best
	avgScore <- round( apply( scoreM, MARGIN=2, mean), digits=2)
	best <- which.max( avgScore)

	# package up the answer
	outSeq <- refFA$seq[ best]
	names(outSeq) <- refFA$desc[ best]
	outScores <- sort( avgScore, decreasing=T)
	out <- list( "ReferenceSequence"=outSeq, "ScoreDetails"=outScores)
	return(out)
}


`constructChromatogramBaseMatrix` <- function( seqInfo, confInfo, refDNA) {

	# given all the sequences and confidence scores that cover one reference construct.
	# populate a matrix that holds all the calls at every location
	refBases <- base::strsplit( refDNA, split="")[[1]]
	nBases <- length( refBases)

	nSeqs <- length( seqInfo)
	seqBases <- base::strsplit( seqInfo, split="")

	# build a matrix for the bases, and one for the confidence scores
	baseM <- matrix( "", nrow=nSeqs, ncol=nBases)
	confM <- matrix( 0, nrow=nSeqs, ncol=nBases)
	colnames(baseM) <- colnames(confM) <- refBases
	rownames(baseM) <- rownames(confM) <- names(seqInfo)

	# we will do pairwise alignment to place each sequence in the best location
	require( Biostrings)
	subM <- nucleotideSubstitutionMatrix()
	for ( i in 1:nSeqs) {
		mySeq <- seqInfo[i]
		myBases <- seqBases[[i]]
		myConfs <- confInfo[[i]]

		# we expect the chromatogram sequences to be smaller than the full length reference,
		# but explicitly see if any sequence extends too far.  If so, truncate now before the main
		# consensus building steps...
		pa <- pairwiseAlignment( mySeq, refDNA, type="local", scoreOnly=F, substitutionMatrix=subM)
		myStart <- start( pattern(pa))
		myStop <- myStart + width( pattern(pa)) - 1
		refStart <- start( subject(pa))
		refStop <- refStart + width( subject(pa)) - 1
		leftTail <- myStart - refStart
		if (leftTail > 0) {
			# crop from head
			mySeq <- substr( mySeq, (leftTail+1), nchar(mySeq))
			myBases <- myBases[ -c(1:leftTail)]
			myConfs <- myConfs[ -c(1:leftTail)]
			myStart <- myStart - leftTail
			myStop <- myStop - leftTail
		}
		myExtra <- nchar(mySeq) - myStop
		refExtra <- nchar(refDNA) - refStop
		rightTail <- myExtra - refExtra
		if (rightTail > 0) {
			# crop from tail
			newLen <- nchar(mySeq) - rightTail
			mySeq <- substr( mySeq, 1, newLen)
			myBases <- myBases[ 1:newLen]
			myConfs <- myConfs[ 1:newLen]
			if (myStop > newLen) myStop <- newLen
		}

		#  now do it for reals
		pa <- pairwiseAlignment( mySeq, refDNA, type="local", scoreOnly=F, substitutionMatrix=subM)
		myStart <- start( pattern(pa))
		myStop <- myStart + width( pattern(pa)) - 1
		refStart <- start( subject(pa))
		refStop <- refStart + width( subject(pa)) - 1
		myAlignSeq <- as.character( alignedPattern(pa))
		refAlignSeq <- as.character( alignedSubject(pa))
		myStartInRef <- refStart - myStart + 1

		# indels in the sequence are fine, but indels in the reference can not be allowed.
		# check for those, and shift/concatenate those bases into previous cell
		hasIndel <- grepl( "-", refAlignSeq, fixed=T)
		if (hasIndel) {
			indelSites <- gregexpr( "-", refAlignSeq, fixed=T)[[1]]
			nGapSoFar <- 0
			for (k in indelSites) {
				# because we shift every time we deal with a gap, the locations keep moving
				kNow <- k - nGapsSoFar
				baseNow <- myBases[kNow]
				kStore <- kNow + 1
				myBases[kStore] <- paste( myBases[kStore], baseNow, sep="")
				myBases <- myBases[ -kNow]
				myConfs <- myConfs[ -kNow]
				nGapSoFar <- nGapSoFar + 1
			}
		}
		myLenNow <- length( myBases)
		myStopInRef <- myStartInRef + myLenNow - 1

		# place the sequence data where it goes
		baseM[ i, myStartInRef:myStopInRef] <- myBases
		confM[ i, myStartInRef:myStopInRef] <- myConfs
	}

	# lastly, turn all the confidence into intergers on 0..100
	if ( all( confM <= 1)) confM <- round( confM * 100)

	# Done!...
	out <- list( "BaseMatrix"=baseM, "ConfidenceMatrix"=confM)
	out
}


`extractChromatogramConsensusSequence` <- function( baseM, confM) {
	
	# given the matrix of all base calls and the matrix of confidence scores, 
	# extract the best consensus
	refBases <- colnames(baseM)
	refLen <- length(refBases)
	outBases <- rep.int( " ", refLen)
	outConfs <- rep.int( 0, refLen)

	# step along, and count up how many votes for each
	for ( i in 1:refLen) {
		myBases <- baseM[ , i]
		myConfs <- as.integer( confM[ , i])
		myVotes <- rep( myBases, times=myConfs)
		# no votes at all is empty, skip
		if ( ! length(myVotes)) next
		bestBase <- names( sort( table( myVotes), decreasing=T))[1]
		bestConf <- round( mean( myConfs[ myBases == bestBase]))
		outBases[i] <- bestBase
		outConfs[i] <- bestConf
	}
	names(outBases) <- names(outConfs) <- refBases

	out <- list( "sequence"=outBases, "confidence"=outConfs)
	out
}


`translateChromatogramSequence` <- function( seqData, cdsBases=NULL, mode=c("vector","string","matrix"),
						badAA="X") {

	# this function translates DNA from these chromatogram tools into AA

	# figure out how many translations we will be doing
	mode <- match.arg( mode)
	if ( mode == "string") {
		seqData <- strsplit( seqData, split="")
		nToDo <- length(seqData)
		N <- length( seqData[[1]])
	} else if ( mode == "matrix") {
		nToDo <- nrow(seqData)
		N <- ncol(seqData)
	} else {
		nToDo <- 1
		N <- length(seqData)
	}

	# were we given the CDS bases? default is full length as one exon
	if ( is.null( cdsBases)) {
		cdsBases <- 1:N
	} else {
		cdsBases <- intersect( 1:N, cdsBases)
	}

	# ready to translate
	out <- rep.int( "", nToDo)
	for ( i in 1:nToDo) {
		if ( mode == "string") {
			myData <- seqData[[i]]
		} else if ( mode == "matrix") {
			myData <- seqData[ i, ]
		} else {
			myData <- seqData
		}
		dnaV <- myData[ cdsBases]
		dnaStr <- paste( dnaV, collapse="")
		aaStr <- DNAtoAA( dnaStr, clipAtStop=F, readingFrame=1)
		aaStr <- gsub( "?", badAA, aaStr, fixed=T)
		out[i] <- aaStr
	}
	return( out)
}


`translateChromatogramConfidence` <- function( confData, cdsBases=NULL, mode=c("vector","matrix"), badAA="X") {

	# this function translates DNA confidence calls from these chromatogram tools into AA confidence calls

	# figure out how many translations we will be doing
	mode <- match.arg( mode)
	if ( mode == "matrix") {
		nToDo <- nrow(confData)
		N <- ncol(confData)
	} else {
		nToDo <- 1
		N <- length(confData)
	}

	# were we given the CDS bases? default is full length as one exon
	if ( is.null( cdsBases)) {
		cdsBases <- 1:N
	} else {
		cdsBases <- intersect( 1:N, cdsBases)
	}

	# ready to translate:  combine triplets by mean
	out <- vector( mode="list", length=nToDo)
	for ( i in 1:nToDo) {
		if ( mode == "matrix") {
			myData <- confData[ i, ]
		} else {
			myData <- confData
		}
		confV <- myData[ cdsBases]

		# take the mean of the triplets, where any zero values denote no data at that base
		lenNow <- length( confV)
		starts <- seq( 1, lenNow, by=3)
		stops <- starts + 2
		stops[ length(stops)] <- lenNow
		aaConf <- mapply( starts, stops, FUN=function(x,y) {
					confValues <- confV[ x : y]
					if ( any( confValues == 0)) return( 0)
					return( mean( confValues))
				})
		aaConf <- round( aaConf)
		# we want to force the confidence to exactly match the AA sequence, so drop any values that got truncated
		if ( badAA == "") {
			drops <- which( aaConf == 0)
			if ( length( drops)) aaConf <- aaConf[ -drops]
		}
		out[[i]] <- aaConf
	}

	if ( nToDo == 1) return( out[[1]])
	return( out)
}
