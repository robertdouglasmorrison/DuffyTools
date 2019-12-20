# chromatogramAssembly.R -- tools to build a single consensus assembly from sets of chromatogram files.  
#				Intended for cases where multiple Sanger sequencing primers are used
# 				to try to reconstruct one larger piece of gene


# top level tool designed to do the entire process of turning many chromatograms 
# covering one gene region into a single best consensus assembly

`chromatogramAssembly` <- function( chromoFiles, referenceFastaFile, doCleaning=TRUE, crop=TRUE,
				curated=TRUE, cdsBases=NULL, min.length=30, min.confidence=0.6, 
				min.DNA.unitscore=0.5, bridgeFastaFile=NULL, verbose=TRUE) {

	# For conserved genes, 'Cleaning' can catch fix chromatogram errors.
	# But not advised for highly variant genes where the isolate likely has many indels

	# Step 0:  load the set of chromatograms into memory
	chromoSet <- loadMultipleChromatograms( chromoFiles, curated=curated)
	nChromo <- length( chromoSet)
	if ( is.null( chromoSet)) return( NULL)
	use <- 1:nChromo
	if ( verbose) cat( "\nGiven as input", nChromo, "chromatograms.")
	
	# Step 1: crop off low confidence edges
	if ( crop) {
		for ( ichr in 1:nChromo) {
			if ( verbose) cat( "\nCropping: ", names(chromoSet)[ichr])
			chromoObj <- chromoSet[[ichr]]
			if ( verbose) cat( "  Was: ", length(chromoObj$PeakPosition))
			chromoObj <- cropChromatogramByConfidence( chromoObj, min.confidence=min.confidence, verbose=verbose)
			if ( is.null( chromoObj)) {
				use <- setdiff( use, ichr)
			} else {
				if ( verbose) cat( "  Now N_Base: ", length(chromoObj$PeakPosition))
				chromoSet[[ichr]] <- chromoObj
			}
		}
	}

	# Step 1.5:   don't let any very short sequences get used
	for ( ichr in use) {
		chromoObj <- chromoSet[[ichr]]
		dnaLen <- nchar( chromoObj$DNA_Calls[1])
		if ( dnaLen < min.length) {
			use <- setdiff( use, ichr)
		}
	}
	# a chance that no chromatograms are useable
	if ( ! length(use)) return( NULL)
	
	# Step 2: load the reference DNA sequence.  We could be given 2+ different references, 
	# so we may need to find the one best reference for this isolate sample
	if ( verbose) cat( "\nSelecting best Reference..")
	refFA <- loadFasta( referenceFastaFile, verbose=FALSE)
	refAns <- selectBestChromatogramReference( chromoSet[use], refFA)
	refDNA <- refAns$ReferenceSequence
	refName <- names(refDNA)[1]
	refAA <- translateChromatogramSequence( refDNA, mode="string", badAA="")
	refStats <- refAns$ScoreDetails
	if (verbose && !is.null( refStats)) {
		cat( "\nTop Scoring Best References:\n")
		print( head( refStats, 10))
	}

	# Step 3:  now with the best reference DNA known, we may want to clean the chromatograms for 
	# typical sequence errors and artifacts
	if ( doCleaning) {
		for ( ichr in use) {
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
	# do every chromatogram, even those flagged for exclusion
	if ( verbose) cat( "\nExtracting best DNA from ", length(use), "chromatograms..")
	for ( ichr in 1:nChromo) {
		chromoObj <- chromoSet[[ichr]]
		chromoName <- names(chromoSet)[ichr]
		smlAns <- selectBestChromatogramSequence( chromoObj, reference=refDNA, type="DNA")
		seqDF <- rbind( seqDF, data.frame( "Chromatogram"=chromoName, smlAns, stringsAsFactors=F))
		seqInfo[ichr] <- mySeq <- smlAns$DNA.Sequence
		myConf <- getChromatogramSequenceConfidence( chromoObj, seq=mySeq)
		confInfo[[ichr]] <- myConf
	}

	# Step 4.5:  Bridge constructs
	# we may have blind spots due to primer and sequencing limitations.   Allow being passed in some
	# non-chromatogram linkage bridging constructs
	nBridge <- 0
	if ( ! is.null( bridgeFastaFile)) {
		if ( verbose) cat( "\nAdding bridge linker sequences..")
		bridgeFA <- loadFasta( bridgeFastaFile, verbose=F)
		# use the REference name to find just those to use
		keep <- grep( refName, bridgeFA$desc, fixed=T)
		bridgeSeqs <- bridgeFA$seq[keep]
		bridgeNames <- paste( "Bridge", bridgeFA$desc[keep], sep="_")
		if ( nBridge <- length( bridgeSeqs)) {
			smlAns <- data.frame( "SequenceName"="Bridge", "ReferenceName"=refName, "DNA.Score"=nchar(bridgeSeqs),
						"DNA.UnitScore"=1.0, "DNA.RefStart"=NA, "DNA.RefStop"=NA, "DNA.Sequence"=bridgeSeqs,
						stringsAsFactors=F)
			seqDF <- rbind( seqDF, data.frame( "Chromatogram"=bridgeNames, smlAns, stringsAsFactors=F))
			for ( k in 1:nBridge) {
				seqInfo[nChromo+k] <- mySeq <- bridgeSeqs[k]
				confInfo[[nChromo+k]] <- rep.int( 1, nchar(mySeq))
			}
		}
		# so now every future step should see these linkers as just other chromatogram sequences
		use <- c( use, (1:nBridge) + nChromo)
		nChromo <- nChromo + nBridge
	}
	rownames(seqDF) <- 1:nrow(seqDF)
	seqDF$DNA.AvgConfidence <- round( sapply( confInfo, mean, na.rm=T) * 100)
	seqDF$AA.Score <- 0;  seqDF$AA.UnitScore <- 0; 
	seqDF$AA.RefStart <- NA;  seqDF$AA.RefStop <- NA; 
	seqDF$AA.AvgConfidence <- 0;  seqDF$AA.Sequence <- ""
	seqDF <- seqDF[ , c( 1:7, 9:15, 8)]

	# do a rough PA on the garbage sequences, just so they show up later in plots, etc.
	notUse <- union( which( seqDF$DNA.UnitScore < min.DNA.unitscore), setdiff( 1:nChromo, use))
	use <- setdiff( 1:nChromo, notUse)
	data( BLOSUM62)
	for ( j in notUse) {
		tmpDNA <- seqDF$DNA.Sequence[j]
		if ( is.na(tmpDNA) || nchar(tmpDNA) < min.length) next
		tmpAA <- DNAtoBestPeptide( tmpDNA)
		if ( is.na(tmpAA) || nchar(tmpAA) < min.length/3) next
		pa <- pairwiseAlignment( tmpAA, refAA, type="local", scoreOnly=F, substitutionMatrix=BLOSUM62)
		aaScore <- score(pa)
		aaStart <- start( subject( pa))
		aaStop <- width( subject( pa)) + aaStart - 1
		seqDF$AA.Score[j] <- aaScore
		seqDF$AA.UnitScore[j] <- round( aaScore / (nchar(tmpDNA)/3), digits=3)
		seqDF$AA.RefStart[j] <- aaStart
		seqDF$AA.RefStop[j]  <- aaStop
		seqDF$AA.AvgConfidence[j] <- round( mean( seqDF$DNA.AvgConfidence[j]))
		seqDF$AA.Sequence[j] <- tmpAA
	}

	# don't always use all, some may be terrible
	# and if absolutely none passed, send back an empty result now
	if ( length( use) <= nBridge) {
		out <- list( "ReferenceName"="FAIL", "ReferenceAA"=refAA, "ReferenceDNA"=as.character(refDNA), "ConsensusAA"="",
			"ConfidenceAA"=integer(0), "EditDistanceAA"=NA,
			"ConsensusDNA"="", "ConfidenceDNA"=integer(0), "EditDistanceDNA"=NA,
			"FragmentDetails"=seqDF, "BaseMatrix"=NULL, "ConfidenceMatrix"=NULL)
		return(out)
	}

	# Step 5:  use all these sequences and confidence scores to build one large matrix of base calls at all locations
	names(seqInfo) <- names(confInfo) <- seqDF$Chromatogram
	if ( verbose) cat( "\nConstructing Consensus DNA matrix using confidence scores..")
	matrixAns <- constructChromatogramBaseMatrix( seqInfo[use], confInfo[use], refDNA)
	baseM <- matrixAns$BaseMatrix
	confM <- matrixAns$ConfidenceMatrix

	# Step 6:  now, use the confidence scores as the voting criteria, to extact the one final consensus
	#          we can use the DNA scores as relative weights
	if ( verbose) cat( "\nExtract Final Consensus DNA..")
	fragWts <- round( seqDF$DNA.UnitScore[use] * 10)
	consensusSeqAns <- extractChromatogramConsensusSequence( baseM, confM, wts=fragWts)
	finalDNAseq <- consensusSeqAns$sequence
	finalDNAstr <- paste( finalDNAseq, collapse="")
	finalDNAconf <- consensusSeqAns$confidence

	# Step 7:  we can now turn the DNA back to AA, for both the final consensus and all the fragments
	if ( verbose) cat( "\nTranslate results to AA..")
	finalAA <- translateChromatogramSequence( finalDNAseq, mode="vector", badAA="X", cdsBases=cdsBases, referenceAA=refAA)
	finalAAconf <- translateChromatogramSeqConfidence( finalDNAseq, finalDNAconf, mode="vector", badAA="X", cdsBases=cdsBases)[[1]]
	# tiny chance of length mismatch, so trap by hand
	nbMin <- min( nchar(finalAA), length(finalAAconf))
	names(finalAAconf)[1:nbMin] <- strsplit( finalAA, split="")[[1]][1:nbMin]
	fragAA <- translateChromatogramSequence( baseM, mode="matrix", badAA=" ", cdsBases=cdsBases, referenceAA=refAA)
	fragAAconf <- translateChromatogramSeqConfidence( baseM, confM, mode="matrix", badAA=" ", cdsBases=cdsBases)

	# Step 8:  we can make the AA scores and locations for all the fragments
	trimmedAA <- gsub( " ", "", fragAA)
	trimmedAA[ nchar(trimmedAA) == 0] <- "X"
	data( BLOSUM62)
	pa <- pairwiseAlignment( trimmedAA, refAA, type="local", scoreOnly=F, substitutionMatrix=BLOSUM62)
	aaScores <- score(pa)
	aaStarts <- start( subject( pa))
	aaStops <- width( subject( pa)) + aaStarts - 1
	seqDF$AA.Score[use] <- aaScores
	seqDF$AA.UnitScore[use] <- round( aaScores / nchar( trimmedAA), digits=3)
	seqDF$AA.RefStart[use] <- aaStarts
	seqDF$AA.RefStop[use]  <- aaStops
	seqDF$AA.AvgConfidence[use] <- round( sapply( fragAAconf, function(x) if ( length(x)) mean(x,na.rm=T) else 0))
	seqDF$AA.Sequence[use] <- fragAA
	# put these all into sequence order
	aaMidpt <- (seqDF$AA.RefStart + seqDF$AA.RefStop) / 2
	ord <- order( aaMidpt)
	seqDF <- seqDF[ ord, ]
	rownames(seqDF) <- 1:nrow(seqDF)

	# also calculate the edit distance, but don't penalize the AA for 'no data' regions
	edDNA <- adist( as.character(refDNA), finalDNAstr)[1]
	edAA <- adist( refAA, finalAA)[1]
	nNdna <- sum( gregexpr( " ", finalDNAstr, fixed=T)[[1]] > 0)
	nXaa <- sum( gregexpr( "X", finalAA, fixed=T)[[1]] > 0)
	edDNA <- max( edDNA - nNdna, 0)
	edAA <- max( edAA - nXaa, 0)
	# there is a tiny chance that the final AA call is all 'no call' Xs.
	if ( grepl( "^X+$", finalAA)) edAA <- NA

	out <- list( "ReferenceName"=refName, "ReferenceAA"=refAA, "ReferenceDNA"=as.character(refDNA), 
			"ConsensusAA"=finalAA, "ConfidenceAA"=finalAAconf, "EditDistanceAA"=edAA, 
			"ConsensusDNA"=finalDNAstr, "ConfidenceDNA"=finalDNAconf, "EditDistanceDNA"=edDNA,
			"FragmentDetails"=seqDF, "BaseMatrix"=baseM, "ConfidenceMatrix"=confM)
	if ( verbose) cat( "\nDone.")
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
	if (nChromo > 10) {
		# see how big they all are, and use the bigger/longer ones
		nBaseEach <- sapply( chromoSet, function(x) return( length( x$PeakPosition)))
		toDo <- order( nBaseEach, decreasing=T)[1:10]
	} else {
		toDo <- 1:nChromo
	}
	
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
	outScores <- base::sort( avgScore, decreasing=T)
	out <- list( "ReferenceSequence"=outSeq, "ScoreDetails"=outScores)
	return(out)
}


`constructChromatogramBaseMatrix` <- function( seqInfo, confInfo, refDNA, wts=NULL) {

	# given all the sequences and confidence scores that cover one reference construct.
	# populate a matrix that holds all the calls at every location
	refBases <- base::strsplit( refDNA, split="")[[1]]
	nBases <- length( refBases)

	nSeqs <- length( seqInfo)
	seqBases <- base::strsplit( seqInfo, split="")

	# build a matrix for the bases, and one for the confidence scores
	baseM <- matrix( " ", nrow=nSeqs, ncol=nBases)
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

		#  now do the pairwise alignment for reals
		pa <- pairwiseAlignment( mySeq, refDNA, type="local", scoreOnly=F, substitutionMatrix=subM)
		myStart <- start( pattern(pa))
		myStop <- myStart + width( pattern(pa)) - 1
		refStart <- start( subject(pa))
		refStop <- refStart + width( subject(pa)) - 1
		myStartInRef <- refStart - myStart + 1
		# we also need the aligned versions of the 2 strings to see where the indel gaps are
		myAlignSeq <- as.character( alignedPattern(pa))
		refAlignSeq <- as.character( alignedSubject(pa))

		# Indels in the chromatogram are allowed.
		# Check for those, and re-assess what the sequence is as needed
		hasSeqIndel <- grepl( "-", myAlignSeq, fixed=T)
		if (hasSeqIndel) {
			# indels shift the bases around, which also dictates that confidence scores shift too
			myAlignBases <- strsplit( myAlignSeq, split="")[[1]]
			indelSites <- which( myAlignBases == "-")
			myAlignConfs <- rep.int( 0, length(myAlignBases))
			goodSites <- setdiff( 1:length(myAlignBases), indelSites)
			# slight chance that value locations may not align perfectly, catch and punish
			myAlignConfs[ goodSites] <- myConfs[ 1:length(goodSites)]
			myAlignConfs[ is.na(myAlignConfs)] <- 0.25  #  25% of a great score...
			#done with shift, put these new values back
			mySeq <- myAlignSeq
			myBases <- myAlignBases
			myConfs <- myAlignConfs
		}

		# Indels in the reference can not be allowed.
		# Check for those, and shift/concatenate those bases into previous cell
		hasRefIndel <- grepl( "-", refAlignSeq, fixed=T)
		if (hasRefIndel) {
			indelSites <- gregexpr( "-", refAlignSeq, fixed=T)[[1]]
			nGapSoFar <- 0
			for (k in indelSites) {
				# because we shift every time we deal with a gap, the locations keep moving
				kNow <- k - nGapSoFar
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


`extractChromatogramConsensusSequence` <- function( baseM, confM, wts=NULL) {
	
	# given the matrix of all base calls and the matrix of confidence scores, 
	# extract the best consensus
	refBases <- colnames(baseM)
	refLen <- length(refBases)
	outBases <- rep.int( " ", refLen)
	outConfs <- rep.int( 0, refLen)

	# we want to downplay any 'N' calls, regardless of how much confidence they had
	confM[ baseM == "N"] <- 1

	# force any weigtht to be integers
	if ( ! is.null( wts)) wts <- as.integer( wts)

	# step along, and count up how many votes for each
	SORT <- base::sort
	TABLE <- base::table
	for ( i in 1:refLen) {
		myBases <- baseM[ , i]
		myConfs <- as.integer( confM[ , i])
		if (any( is.na( myConfs))) {
			cat( "\nDebug: bases: ", myBases)
			cat( "\nDebug: confs: ", myConfs)
		}
		if (is.null(wts)) {
			myVotes <- rep( myBases, times=myConfs)
		} else {
			myVotes <- rep( myBases, times=(myConfs*wts))
		}
		# no votes at all is empty, skip
		if ( ! length(myVotes)) next
		bestBase <- names( SORT( TABLE( myVotes), decreasing=T))[1]
		bestConf <- round( mean.default( myConfs[ myBases == bestBase]))
		outBases[i] <- bestBase
		outConfs[i] <- bestConf
	}
	names(outBases) <- names(outConfs) <- refBases

	out <- list( "sequence"=outBases, "confidence"=outConfs)
	out
}


`translateChromatogramSequence` <- function( seqData, cdsBases=NULL, mode=c("vector","string","matrix"),
						badAA="X", referenceAA=NULL) {

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
	# keep in mind these are in the units of the reference
	if ( is.null( cdsBases)) {
		cdsBases <- 1:N
	} else {
		cdsBases <- intersect( 1:N, cdsBases)
	}

	# ready to translate
	out <- rep.int( "", nToDo)
	for ( i in 1:nToDo) {
		# we need to bear in mind that these DNA constructs may have indel gaps
		myCdsBases <- cdsBases
		if ( mode == "string") {
			myData <- seqData[[i]]
			indelSites <- as.numeric( gregexpr( "-", myData, fixed=T)[[1]])
			if ( sum( indelSites > 0)) myCdsBases <- setdiff( myCdsBases, indelSites)
		} else if ( mode == "matrix") {
			myData <- seqData[ i, ]
			indelSites <- which( myData == "-")
			if ( length(indelSites)) myCdsBases <- setdiff( myCdsBases, indelSites)
		} else {
			myData <- seqData
			indelSites <- which( myData == "-")
			if ( length(indelSites)) myCdsBases <- setdiff( myCdsBases, indelSites)
		}
		dnaV <- myData[ myCdsBases]
		dnaStr <- paste( dnaV, collapse="")

		# in spite of all best efforts, we can't garauntee that everything stays perfectly in frame
		# so use the most robust method
		# aaStr.plain <- DNAtoAA( dnaStr, clipAtStop=F, readingFrame=1:3)
		aaStr.best <- DNAtoFrameShiftingPeptides( dnaStr, referenceAA=referenceAA, details=F)
		# keep the full length that is most like what we expect
		# dm <- adist( aaStr.plain, aaStr.best)
		# aaStr <- aaStr.plain[ which.min( dm)]
		aaStr <- aaStr.best
		aaStr <- gsub( "?", badAA, aaStr, fixed=T)

		# the frame shifting tool uses 'X' for dead/missing data
		if ( badAA != "X") aaStr <- gsub( "X", badAA, aaStr, fixed=T)

		# done
		out[i] <- aaStr
	}
	return( out)
}


`translateChromatogramSeqConfidence` <- function( seqData, confData, cdsBases=NULL, mode=c("vector","matrix"), badAA="X") {

	# this function translates DNA confidence calls from these chromatogram tools into AA confidence calls
	# (but we need the sequence data too to find any indel gaps

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
		# we need to bear in mind that these DNA constructs may have indel gaps
		myCdsBases <- cdsBases
		if ( mode == "matrix") {
			myData <- confData[ i, ]
			mySeq <- seqData[ i, ]
			indelSites <- which( mySeq == "-")
			if ( length(indelSites)) myCdsBases <- setdiff( myCdsBases, indelSites)
		} else {
			myData <- confData
			mySeq <- seqData
			indelSites <- which( mySeq == "-")
			if ( length(indelSites)) myCdsBases <- setdiff( myCdsBases, indelSites)
		}
		confV <- myData[ myCdsBases]

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
	return( out)
}

