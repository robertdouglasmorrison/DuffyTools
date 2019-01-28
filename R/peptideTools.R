# peptideTools.R -  peptide based tools


`DNAtoBestPeptide` <- function( dnaSet, clipAtStop=FALSE, readingFrames=1:6,
				tieBreakMode=c("evalue","sample","reference"), reference=NULL) {

	tieBreakMode <- match.arg( tieBreakMode)
	if ( tieBreakMode == "reference") {
		require( Biostrings)
		data( BLOSUM62)
		refAA <- toupper( as.character( reference))
	}

	out <- sapply( as.character( dnaSet), function( dna) {

			pepsIn <- DNAtoAA( dna, clipAtStop=clipAtStop, readingFrames=readingFrames)

			# if ignoring stop codons, break into all coding fragments of each
			if ( ! clipAtStop) {
				pepFrags <- unlist( strsplit( pepsIn, split=STOP_CODON_PATTERN, fixed=FALSE), use.names=F)
			} else {
				pepFrags <- pepsIn
			}

			if ( tieBreakMode == "sample") {
				# longest wins
				nc <- nchar( pepFrags)
				pepFrags <- pepFrags[ nc == max(nc)]
				if ( length(pepFrags) > 1) pepFrags <- sample( pepFrags, size=1)
				return( pepFrags[1])
			}

			# other modes will consider the 'longer' subset
			nc <- nchar( pepFrags)
			pepFrags <- pepFrags[ nc > (max(nc)*0.9)]
			if ( length(pepFrags) < 2) return( pepFrags[1])

			if ( tieBreakMode == "reference") {
				paScores <- pairwiseAlignment( pepFrags, refAA, type="local", substitutionMatrix=BLOSUM62, scoreOnly=T)
				return( pepFrags[ which.max( paScores)])
			}

			if ( tieBreakMode == "evalue") {
				return( pepFrags[ chooseReadingFrame( pepFrags, allFrames=pepsIn)])
			}


		}, USE.NAMES=FALSE)

	return( out)
}


`DNAtoCleanChromatogramDNA` <- function( dna, referenceDNA, referenceAA=NULL, verbose=T) {

	require( Biostrings)

	# DNA from chromatograms can have repeats or gaps that destroy reading frame.
	# try to find and fix
	dnaNow <- dna[1]

	# accumulate a set of locations already visited
	maskBases <- vector()

	# and try to build a table of details if the caller wants to know exactly what changed
	fixLoc <- fixWas <- fixNew <- vector()
	nFix <- 0

	repeat {
		pa <- pairwiseAlignment( dnaNow, referenceDNA, type="local", scoreOnly=F)
		dnaStart <- start( pattern( pa))
		dnaStr <- as.character( pattern( pa))
		refStart <- start( subject( pa))
		refStr <- as.character( subject( pa))
		scoreNow <- score(pa)
		# if we don't have a good score, don't even try
		if ( scoreNow < nchar(dnaNow)*0.1) break

		# do we see any context that looks like a 1bp gap?
		# use a fixed size flank of 9bp
		N_FLANK <- 9
		dnaGapAns <- as.numeric( gregexpr( "[ACGTN]{9}\\-[ACGTN]{9}", dnaStr)[[1]])
		refGapAns <- as.numeric( gregexpr( "[ACGTN]{9}\\-[ACGTN]{9}", refStr)[[1]])

		# don't revisit any sites already done
		if ( length( maskBases)) {
			dnaGapAns <- setdiff( dnaGapAns, maskBases)
			if( ! length( dnaGapAns)) dnaGapAns <- -1
			refGapAns <- setdiff( refGapAns, maskBases)
			if( ! length( refGapAns)) refGapAns <- -1
		}

		# in the case there is more than one, pick the one nearest the middle
		if ( length(dnaGapAns) > 1) {
			midPt <- nchar( dnaStr) / 2
			dx <- abs( dnaGapAns + N_FLANK - midPt)
			best <- which.min( dx)
			dnaGapAns <- dnaGapAns[ best]
		}
		if ( length(refGapAns) > 1) {
			midPt <- nchar( refStr) / 2
			dx <- abs( refGapAns + N_FLANK - midPt)
			best <- which.min( dx)
			refGapAns <- refGapAns[ best]
		}
		if ( dnaGapAns < 0 && refGapAns < 0) break

		# which one? and does removing it make things better?
		if ( refGapAns > 0) {
			# gap in reference means extra base in chromatogram
			gapLocation <- refGapAns + N_FLANK
			nOtherGaps <- sum( gregexpr( "-", substr( dnaStr, 1, refGapAns), fixed=T)[[1]] > 0)
			dnaLocation <- dnaStart + gapLocation - 1 - nOtherGaps
			newLeft <- substr( dnaNow, 1, dnaLocation-1)
			newRight <- substr( dnaNow, dnaLocation+1, nchar(dnaNow))
			dnaNew <- paste( newLeft, newRight, sep="")

			# is this better?
			if ( is.null( referenceAA)) {
				scoreNew <- pairwiseAlignment( dnaNew, referenceDNA, type="local", scoreOnly=T)
			} else {
				aaNow <- DNAtoBestPeptide(dnaNow, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				aaNew <- DNAtoBestPeptide(dnaNew, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				scoreNow <- pairwiseAlignment( aaNow, referenceAA, type="local", scoreOnly=T)
				scoreNew <- pairwiseAlignment( aaNew, referenceAA, type="local", scoreOnly=T)
			}
			extraBase <- substr( dnaNow, dnaLocation, dnaLocation)
			oldContext <- substr( dnaNow, dnaLocation-1, dnaLocation+1)
			newContext <- paste( substr(oldContext,1,1), substr(oldContext,3,3), sep="")
			pctBetter <- round( (scoreNew-scoreNow) * 100 / scoreNow, digits=2)
			if (verbose) {
				cat( "\n--------------------------")
				cat( "\nCleaning Chromatogram DNA:  found extra base in Chromatogram.")
				cat( "\nDNA Context:", substr(dnaStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nRef Context:", substr(refStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nExtra Base: ", extraBase,"   Was: ", oldContext, "   Now: ", newContext)
			}
			if ( scoreNew > scoreNow) {
				if (verbose) cat( "\nScore:    Was: ", scoreNow, "   Now: ", scoreNew, "   PctImprove: ", pctBetter, "%")
				dnaNow <- dnaNew
				nFix <- nFix + 1
				fixLoc[nFix] <- dnaLocation
				fixWas[nFix] <- extraBase
				fixNew[nFix] <- ""
			} else {
				if (verbose) cat( "\nBut Score no better!    Was: ", scoreNow, "   Now: ", scoreNew)
			}
			if (verbose) cat( "\n--------------------------\n")
			maskBases <- sort( unique( c( maskBases, (refGapAns - N_FLANK %/% 2):(refGapAns + N_FLANK %/% 2))))
			next
		}
		if ( dnaGapAns > 0) {
			# gap in chromatogram means extra base in reference 
			gapLocation <- dnaGapAns + N_FLANK
			nOtherGaps <- sum( gregexpr( "-", substr( dnaStr, 1, dnaGapAns), fixed=T)[[1]] > 0)
			extraRefBase <- substr( refStr, gapLocation, gapLocation)
			dnaLocation <- dnaStart + gapLocation - 1 - nOtherGaps
			newLeft <- substr( dnaNow, 1, dnaLocation-1)
			# since the gap is in my string, don't step past it
			#newRight <- substr( dnaNow, dnaLocation+1, nchar(dnaNow))
			newRight <- substr( dnaNow, dnaLocation, nchar(dnaNow))
			dnaNew <- paste( newLeft, extraRefBase, newRight, sep="")

			# is this better?
			if ( is.null( referenceAA)) {
				scoreNew <- pairwiseAlignment( dnaNew, referenceDNA, type="local", scoreOnly=T)
			} else {
				aaNow <- DNAtoBestPeptide(dnaNow, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				aaNew <- DNAtoBestPeptide(dnaNew, readingFrames=1:3, tieBreakMode="reference", reference=referenceAA)
				scoreNow <- pairwiseAlignment( aaNow, referenceAA, type="local", scoreOnly=T)
				scoreNew <- pairwiseAlignment( aaNew, referenceAA, type="local", scoreOnly=T)
			}
			oldContext <- substr( dnaNow, dnaLocation-1, dnaLocation)
			newContext <- paste( substr(oldContext,1,1), extraRefBase, substr(oldContext,2,2), sep="")
			pctBetter <- round( (scoreNew-scoreNow) * 100 / scoreNow, digits=2)
			if (verbose) {
				cat( "\n--------------------------")
				cat( "\nCleaning Chromatogram DNA:  found missing base in Chromatogram.")
				cat( "\nDNA Context:", substr(dnaStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nRef Context:", substr(refStr, gapLocation-N_FLANK,gapLocation+N_FLANK))
				cat( "\nMissing Base: ", extraRefBase,"   Was: ", oldContext, "   Now: ", newContext)
			}
			if ( scoreNew > scoreNow) {
				if (verbose) cat( "\nScore:    Was: ", scoreNow, "   Now: ", scoreNew, "   PctImprove: ", pctBetter, "%")
				dnaNow <- dnaNew
				nFix <- nFix + 1
				# let's call it as being appended to previous instead of prepended to the current base
				fixLoc[nFix] <- dnaLocation - 1
				fixWas[nFix] <- substr(oldContext,1,1)
				fixNew[nFix] <- substr(newContext,1,2)
			} else {
				if (verbose) cat( "\nBut Score no better!    Was: ", scoreNow, "   Now: ", scoreNew)
			}
			if (verbose) cat( "\n--------------------------\n")
			maskBases <- sort( unique( c( maskBases, (dnaGapAns - N_FLANK %/% 2):(dnaGapAns + N_FLANK %/% 2))))
			next
		}
	}
	
	dnaOut <- dnaNow
	detailsOut <- NULL
	if (nFix) {
		detailsOut <- data.frame( "Location"=fixLoc, "PreviousBase"=fixWas, "CleanedBase"=fixNew,
					stringsAsFactors=F)
		# note that every cleaning modification we did throws off the locations w.r.t. the original DNA
		# we want the 'Location' data to be in caller's original units
		if ( nFix > 1) {
			ord <- order( detailsOut$Location)
			detailsOut <- detailsOut[ ord, ]
			rownames(detailsOut) <- 1:nFix
			curDelta <- nchar(detailsOut$PreviousBase[1]) - nchar(detailsOut$CleanedBase[1])
			for ( i in 2:nFix) {
				detailsOut$Location[i] <- detailsOut$Location[i] + curDelta
				curDelta <- curDelta + nchar(detailsOut$PreviousBase[i]) - nchar(detailsOut$CleanedBase[i])
			}
		}
	}

	return( list( "DNA"=dnaOut, "CleaningDetails"=detailsOut))
}


`peptide.Evalue` <- function( pepSet) {

	aaScores <- c("A"=4, "B"=4, "C"=9, "D"=6, "E"=5, "F"=6, "G"=6, "H"=8, "I"=4, "J"=3, "K"=5, "L"=4,
			"M"=5, "N"=6, "O"=0, "P"=7, "Q"=5, "R"=5, "S"=4, "T"=5, "U"=0, "V"=4, "W"=11,
			"X"= -1, "Y"=7, "Z"=4, "*"=1, "?"= -1)

	bases <- strsplit( pepSet, split="")
	scores <- sapply( bases, function(x) {
				wh <- match( x, names(aaScores), nomatch=0)
				return( sum( aaScores[ wh]))
			})
	
	# don't let empty peptides have a good score
	scores[ nchar(pepSet) < 1] <- 1

	return( 1/scores)
}


`chooseReadingFrame` <- function( aaSet, allFrames=aaSet) {

	aaLevels <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
			"M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W",
			"X", "Y", "Z", "*", "?")

	# get the overall distribution of aa, knowing that 5 of the 6 reading frames are not the right one
	aaChars <- strsplit( aaSet, split="")
	totalDist <- table( factor( unlist( strsplit(allFrames,split="")), levels=aaLevels))

	# see how similar each is to the overall distribution
	pvals <- sapply( aaChars, FUN=function(x) {
			if ( ! length(x)) return(0)
			thisDist <- table( factor( x, levels=aaLevels))
			return( cor.test( totalDist, thisDist)$p.value)
		})

	return( which.max( as.numeric(pvals)))
}


`DNAtoFrameShiftingPeptides` <- function( dnaSet, min.aa.length=10) {

	dnaIn <- as.character( dnaSet)
	nIn <- length( dnaIn)
	outList <- vector( mode="list", length=nIn)

	for ( i in 1:nIn) {
		dna <- dnaIn[ i]
		outFrame <- outPeptide <- outLength <- outDnaStart <- outDnaStop <- vector()
		N <- 0
		# repeatedly find the best peptide remaining in the DNA string
		repeat {
			pepStrings <- DNAtoAA( dna, clipAtStop=FALSE, readingFrames=1:3)
			pepTerms <- strsplit( pepStrings, split=STOP_CODON_PATTERN, fixed=FALSE)
	
			# which reading frame has the best?
			bigLength <- sapply( pepTerms, function(x) return( max( nchar(x))))
			bigPepTerm <- sapply( pepTerms, function(x) return( x[ which.max( nchar(x))]))
			bestFrame <- which.max( bigLength)
			bestPep <- bigPepTerm[ bestFrame]
			bestLen <- bigLength[ bestFrame]
			if ( bestLen < min.aa.length) break
	
			N <- N + 1
			outFrame[N] <- bestFrame
			outPeptide[N] <- bestPep
			outLength[N] <- bestLen
	
			# where was this peptide in the DNA, 'regexpr' has a limit on size of pattern
			if ( nchar( bestPep) < 999) {
				where <- regexpr( bestPep, pepStrings[bestFrame])[1]
			} else {
				where <- start( subject( pairwiseAlignment( bestPep, pepStrings[bestFrame], 
						type="global-local")))
			}
			dnaStart <- (where-1) * 3 + bestFrame
			dnaStop <- dnaStart + (bestLen * 3) - 1
			outDnaStart[N] <- dnaStart
			outDnaStop[N] <- dnaStop
	
			# mask it out and go around again
			polyN <- paste( rep.int( "N", (dnaStop-dnaStart+1)), collapse="")
			substr( dna, dnaStart, dnaStop) <- polyN
		}
	
		out <- data.frame( "Peptide"=outPeptide, "Frame"=outFrame, "AA_Length"=outLength, 
				"DNA_Start"=outDnaStart, "DNA_Stop"=outDnaStop, stringsAsFactors=F)
		if ( nrow(out)) {
			ord <- order( out$DNA_Start, decreasing=F)
			out <- out[ ord, ]
			rownames(out) <- 1:nrow(out)
		}
		outList[[i]] <- out
	}
	if ( nIn > 1) return( outList)
	return( outList[[1]])
}


`peptide2BestProtein` <- function( peptides, proteinsFastaFile, nBest=1, tieBreakMode=c("sample", "all", "topN"),
				substitutionMatrix=NULL, details=FALSE, verbose=T) {

	require( Biostrings)
	if( is.null( substitutionMatrix)) {
		data( PAM70MS, envir=environment())
		substitutionMatrix <- PAM70MS
	}
	tieBreakMode <- match.arg( tieBreakMode)

	# load the protein targets
	# allow it to be a pre-loaded FASTA object
	if ( is.character( proteinsFastaFile)) {
		if ( ! file.exists( proteinsFastaFile)) {
			cat( "\nProteins FASTA file not found: ", proteinsFastaFile)
			return( NULL)
		}
		fa <- loadFasta( proteinsFastaFile, verbose=verbose, short.desc=F)
	} else if ( is.list( proteinsFastaFile)) {
		if ( ! all( c("desc", "seq") %in% names( proteinsFastaFile))) {
			cat( "\nProteins FASTA object missing 'desc' and/or 'seq' fields.")
			return( NULL)
		}
		fa <- proteinsFastaFile
	} else {
		cat( "\nInvalid 'proteinsFastaFile' object.")
		return( NULL)
	}

	# use AAStrings for speedup
	protstrings <- AAStringSet( fa$seq)
	nProteins <- length( protstrings)
	protNames <- fa$desc
	pepstrings <- AAStringSet( peptides)
	nPeptides <- length( peptides)

	# a second method of speedup using the 'PDict' functions, to pre-find perfect matches
	# one giant string is the 'dictionary', and we make a pointer list to extract who is where
	proteinLengths <- nchar( fa$seq)
	spacerString <- "XXXXXXXXXX"
	proteinStarts <- c( 1, (cumsum( proteinLengths + nchar(spacerString)) + 1)[ 1:(nProteins-1)])
	giantPDictString <- AAString( paste( fa$seq, collapse=spacerString))
	if (verbose) cat( "\nPre-scanning for perfect matches..")
	pdAns <- matchPDict( pepstrings, giantPDictString)
	pdStarts <- startIndex( pdAns)
	pdLengths <- sapply( pdStarts, length)
	if (verbose) cat( "  ", sum( pdLengths > 0), "out of", nPeptides, "\n")


	# local function to test one peptide against all sequences, returning the N best local alignments
	bestHit <- function( ipep, nBest=1, tieBreakMode=c("all", "sample", "topN"), details=FALSE, verbose=T) { 
	
		tieBreakMode <- match.arg( tieBreakMode)

		peptide <- pepstrings[ipep]
		thisPep <- as.character(peptide)
		thisLen <- nchar( thisPep)

		# use the pre-screen PDict answer to perhaps use only a local subset of the proteins
		if ( pdLengths[ipep] > 0) {
			pdHits <- sort( unique( fastFindInterval( pdStarts[[ipep]], proteinStarts)))
			protstrings <- protstrings[ pdHits]
			protNames <- protNames[ pdHits]
		} else {
			pdHits <- 0
		}

		# find the local scores for all
		scores <- pairwiseAlignment( protstrings, peptide, type="local", 
					substitutionMatrix=substitutionMatrix, 
					gapOpening=-5, gapExtension=-3, scoreOnly=T); 
		best <- which.max(scores); 
		topScore <- scores[ best]
		topName <- protNames[best]
		topScorePerAA <- topScore / thisLen
		# if more than K all are 'best', choose which (subset?) to keep
		whoBest <- which( scores == topScore)
		bestPtrs <- whoBest
		nBestReported <- nBestFound <- length( bestPtrs)
		# default is same as 'all'
		if ( length( whoBest) > nBest) {
			if ( tieBreakMode == "sample") {
				whoBest <- sample( whoBest, size=nBest)
				bestPtrs <- whoBest
				nBestReported <- length( bestPtrs)
			}
		} else {
			if ( tieBreakMode == "topN") {
				ord <- order( scores, decreasing=T)
				whoBest <- ord[ 1:nBest]
				bestPtrs <- whoBest
				nBestReported <- nBestFound <- length( bestPtrs)
			}
		}
		bestScores <- scores[ bestPtrs]
		bestNames <- protNames[ bestPtrs]
		bestScoresPerAA <- bestScores / thisLen
		#if (verbose) cat( "\r", ipep, as.character(peptide), bestScores[1])
		if (verbose) {
			if ( length( pdHits) > 5) pdHits <- pdHits[1:5]
			cat( "\r", ipep, "PDict:", pdHits, "  Score:", bestScores[1])
		}
		
		if ( ! details) {
			out <- list( "Score"=topScore, "ScorePerAA"=topScorePerAA, "Protein"=bestNames)
			return( out)
		}

		# details mode returns more details after complete alignment
		bestPattFrag <- bestSubjFrag <- proteinRegion <- rep.int( "", nBestReported)
		protStart <- protStop <- editDist <- rep.int( NA, nBestReported)
		for (j in 1:nBestReported) {
			ans <- pairwiseAlignment( peptide, protstrings[ bestPtrs[j]], type="local", 
						substitutionMatrix=substitutionMatrix, 
						gapOpening=-5, gapExtension=-3, scoreOnly=F)
			bestPattFrag[j] <- as.character( pattern( ans))
			bestSubjFrag[j] <- as.character( subject( ans))
			from <- protStart[j] <- start( subject(ans)) - start( pattern(ans)) + 1
			to <- protStop[j] <- protStart[j] + width( peptide) - 1
			thatPep <- proteinRegion[j] <- substr( as.character( protstrings[ bestPtrs[j]]), from, to)
			editDist[j] <- stringDist( c( thisPep, thatPep))[1]
		}
	
		out <- data.frame("PeptideID"=rep.int(ipep,nBestReported), "Peptide"=as.character(peptide), 
				"ProtName"=bestNames, "Score"=bestScores, 
				"ScorePerAA"=bestScoresPerAA, "Weight"=1/nBestFound, 
				"PepFrag"=bestPattFrag, "ProtFrag"=bestSubjFrag, 
				"ProtStart"=protStart, "ProtStop"=protStop, 
				"ProtSeq"=proteinRegion, "EditDistance"=editDist, 
				row.names=1:nBestReported, stringsAsFactors=F)
		return(out)
	}


	bigP2ProAns <- multicore.lapply( 1:nPeptides, bestHit, nBest=nBest, tieBreakMode=tieBreakMode,
				details=details, preschedule=TRUE, verbose=verbose)

	return( bigP2ProAns)
}


`lowComplexityPeptides` <- function( peps, min.aa=4, max.top1.pct=51, max.top2.pct=76, 
				max.top3.pct=95, dropPattern="KKKKK|FFFFF|IYIYIY") {

	# given a set of peptide, and some rules for complexity, flag which are too low complexity
	isLOW <- vector()

	# part 1:  any test that work on strings
	if ( ! is.null( dropPattern)) {
		hits <- grep( dropPattern, peps)
		if ( length( hits)) {
			isLOW <- c( isLOW, hits)
		}
	}

	TABLE <- base::table
	WHICH <- base::which
	SAPPLY <- base::sapply
	UNION <- base::union

	# part 2:  any tests that need vectors of single characters
	lens <- nchar( peps)
	pepVs <- strsplit( peps, split="")
	pepAAtables <- lapply( pepVs, function(x) sort.int( TABLE(x), decreasing=T))
	AAtableLens <- SAPPLY( pepAAtables, length)

	tooFewAA <- WHICH( AAtableLens < min.aa)
	isLOW <- UNION( isLOW, tooFewAA)

	top1count <- SAPPLY( pepAAtables, function(x) x[1])
	top2count <- SAPPLY( pepAAtables, function(x) sum( x[1:2], na.rm=T))
	top3count <- SAPPLY( pepAAtables, function(x) sum( x[1:3], na.rm=T))

	tooTop1 <- WHICH( top1count*100/lens > max.top1.pct)
	isLOW <- UNION( isLOW, tooTop1)
	tooTop2 <- WHICH( top2count*100/lens > max.top2.pct)
	isLOW <- UNION( isLOW, tooTop2)
	tooTop3 <- WHICH( top3count*100/lens > max.top3.pct)
	isLOW <- UNION( isLOW, tooTop3)

	isLOW <- sort.int( isLOW)
	return( isLOW)
}


peptideChopper <- function( proteinSet, len=15, overlap=9, remove.duplicates=TRUE, verbose=TRUE) {

	# find all the peptides of length "len" with overlap "overlap" in the
	# given protein sequences...

	if (verbose) cat( "\n\n-----------------------\nPeptide Length: ", len, 
				"\tOverlap: ", overlap, "\n")

	# visit each protein one at a time
	nProtein <- length( proteinSet)
	out <- vector( mode="list", length=nProtein)

	for (i in 1:length( proteinSet)) {
		pseq <- proteinSet[i]
		nAA <- nchar( pseq)

		# turn length and overlap into our step size and where to end
		stepSize <- len - overlap
		lastStart <- nAA - len + 1

		# make all the begin:end pairs
		pepStart <- seq( 1, nAA, by=stepSize)
		pepEnd <- pepStart + len - 1
		npeps <- length( pepStart)

		# after the end, make one last peptide that has greater overlap, just to
		# make sure we end on the last a.a. of the protein
		perfectEnding <- ( pepEnd[ npeps] == nAA)
		if ( ! perfectEnding) {
			npeps <- npeps + 1
			pepStart[npeps] <- (nAA-len+1)
			pepEnd[npeps] <- nAA
		}

		# extract all those peptides
		pepSet <- substr( rep.int( pseq, npeps), start=pepStart, stop=pepEnd)

		# eliminate duplicates
		nBefore <- length( pepSet)
		if ( remove.duplicates) {
			duplicates <- which( duplicated( pepSet))
			if ( length( duplicates) > 0) {
				pepSet <- pepSet[ -duplicates]
				pepStart <- pepStart[ -duplicates]
				pepEnd <- pepEnd[ -duplicates]
			}
		}
		nAfter <- length( pepSet)
		if (verbose) cat( "\rGene: ", names(proteinSet)[i], "\tN_Peptides: ", nAfter)
		if (verbose && remove.duplicates) cat( "\tN_Duplicates_Removed:", nBefore-nAfter)
		
		# make a small data frame that has "where info"
		myDF <- data.frame( "Peptide"=pepSet, "Position"=pepStart, "End"=pepEnd, stringsAsFactors=FALSE)
		out[[i]] <- myDF
	}
	names( out) <- names(proteinSet)
	return( out)
}


peptideTable <- function( gSet, len=15, overlap=9, verbose=TRUE) {

	out <- peptideChopper( gSet, len, overlap, verbose=verbose)
	# just the list of data frames
	ans <- out$peptideSets

	# build a data.frame...
	nrows <- 0
	ncols <- length( ans)
	totalPeps <- 0
	for ( i in 1:ncols) {
		# "chopper" returns DFs now, not just vectors...
		myDF <- ans[[i]]
		nr <- length( myDF$Peptide)
		nrows <- max( c(nrows, nr))
		totalPeps <- totalPeps + nr
	}
	peps <- matrix( "", nrow=nrows, ncol=ncols)
	allPepNames <- vector()
	for ( i in 1:ncols) {
		myDF <- ans[[i]]
		pset <- myDF$Peptide
		np <- length( pset)
		peps[ 1:np, i] <- pset
		allPepNames <- union( allPepNames, pset)
	}
	colnames(peps) <- names( ans)
	write.table( peps, file="peptideTable.tsv", sep="\t", quote=FALSE)
	cat("\n\nSummary:  Length=", len, "\tOverlap=", overlap)
	cat("\n\tTotal Peptides=", totalPeps, "\tTotal Unique Peptides=", length(allPepNames),"\n")
}


peptideTableSet <- function( gSet, len=15, overlap=c(7:9)) {

	# call the peptide table and chopper iteratively to make a data frame of counts
	cnts <- matrix( nrow=length( gSet), ncol=length(overlap))

	for( i in 1:length(overlap)) {
		thislap <- overlap[i]
		out <- peptideChopper( gSet, len, thislap, remove.duplicates=FALSE)
		# we want the counts only
		cnts[ ,i] <- out$peptideCounts
	}

	myDF <- data.frame( gSet, cnts, stringsAsFactors=FALSE)
	colnames( myDF) <- c( "Gene", paste( "Overlap_", overlap, sep=""))
	rownames( myDF) <- 1:nrow( myDF)
	write.table( myDF, file="peptideOverview.txt", sep="\t", quote=FALSE)
}


peptideViewer <- function( gene, AAoffset=0, len=15, overlap=9, overlayDF=NULL, asPNG=FALSE) {

	# draw where the peptides land...
	# protein <- extractProteinSequenceFromWeb( gene, verbose=FALSE)
	protein <- gene2Protein( gene)
	out <- peptideChopper( gene, len=len, overlap=lap, remove.duplicates=TRUE)
	ans <- out$peptideSets

	if (asPNG) {
		jpeg( file=paste( gene, ".jpg", sep=""), width=1000, height=650, pointsize=12)
	}

	# plot setup
	myCex <- 1.0
	par( family="mono")
	nc <- nchar( protein[1])
	ncDraw <- par("pin") / par("cin") * 2

	# try to overlay other peptides to see where they land...  Christian's for now
	overlay <- FALSE
	if ( ! is.null( overlayDF)) {
		otherPeps <- subset( overlayDF, GeneID==gene)
		if ( nrow( otherPeps) > 0) {
			overlay <- TRUE
			otherStarts <- vector()
			for( i in 1:nrow( otherPeps)) otherStarts[i] <- regexpr( otherPeps$Sequence[i], protein, fixed=TRUE)
			missing <- which( otherStarts < 0)
			if (length( missing) > 0) cat( "\nN_overlay peptides not matching protein", length(missing))
		}
	}

	# total number of plots based on how many aa we can fit across X axis
	nPlots <- ceiling( nc / ncDraw)[1]
	nPlotsPerScreen <- 3
	par(mfcol=c(nPlotsPerScreen,1))
	for ( iplot in 1:nPlots) {
		thisXlim <- c( ((iplot-1)*ncDraw[1]), iplot*ncDraw[1]) + AAoffset
		plot( 0,0, xlim=thisXlim, ylim=c(0,5),type="n", xlab=paste( gene, "  amino acid position"), ylab="")
	
		# blindly draw the whole protein at the bottom
		text( (1:nc), rep(0,times=nc), labels=unlist(strsplit(protein,"")), adj=0, cex=myCex, col=1)
		xwid <- 1.0
		ygap <- 0.7

		# step through the peptides and draw them
		maxDeep <- floor( len / (len-lap)) + 1
		myDF <- ans[[1]]
		for ( i in 1:nrow(myDF)) {
			pep <- myDF$Peptide[i]
			ifrom <- myDF$Position[i]
			ito <- myDF$End[i]
			ordinal <- myDF$Ordinal[i]
			curDeep <- ( (ordinal-1) %% maxDeep) + 1
			text( (((ifrom-1)*xwid)+1), (curDeep*ygap), labels=pep, adj=0, cex=myCex, col=4)
			if ( ordinal %% 10 == 0) {
				text( (((ifrom-1)*xwid)+2), ((curDeep+0.5)*ygap), labels=ordinal, adj=0, cex=myCex*0.95, col=2)
			}
		}

		# now try to overlay the 'others'
		if ( overlay) {
			overSet <- which( (otherStarts > 0) & (otherStarts >= (thisXlim[1]-len)) & (otherStarts <= thisXlim[2]))
			if ( length( overSet) > 0) {
				text( (((otherStarts[overSet]-1)*xwid)+1), ((maxDeep+1)*ygap), 
					labels=otherPeps$Sequence[ overSet], adj=0, cex=myCex, col=2)
				text( (((otherStarts[overSet]-1)*xwid)+2), ((maxDeep+2)*ygap), 
					labels=otherPeps$PeptideID[overSet], adj=0, cex=myCex, col=2)
			}
		}

		# do we need to pause for screen capture?
		if ( asPNG) {

		} else {
			if ( iplot %% nPlotsPerScreen == 0 && iplot < nPlots) {
				cat( "\npausing for screeen capture.  Click on plot to continue.")
				locator(1)
			}
		}
	}
	if ( asPNG) dev.off()
}


peptideRedundancySniffer <- function() {

	# this will see if we can save many peptides from being ordered, by checking
	# for big partial match from one end of peptides
	cat( "\nReading in peptide table:     peptideTable.txt")
	pepTable <- read.delim( "peptideTable.txt", as.is=TRUE)

	# gather up all peptides as one big set
	pepSet <- vector()
	for ( i in 1:ncol(pepTable)) {
		oneSet <- pepTable[ ,i]
		# because the table is ragged, some are empty slots
		good <- which( oneSet != "")
		cat( "\nReading Gene: ", colnames(pepTable)[i], "\tN_peptides in: ", length( good))
		pepSet <- append( pepSet, oneSet[ good])
	}
	cat( "\nTotal Peptides read in: ", length( pepSet))
	cat( "\nTotal Unique Peptides (full length): ", length( unique( pepSet)))

	fullLen <- nchar( pepSet[1])
	for ( iLess in 1:3) {
		iend <- fullLen - iLess
		smlSet <- substr( pepSet, 1, iend)
		cat( "\n\nUnique over N-terminal ",iend, "-mer:  ", length( unique( smlSet)), sep="")
	}
	for ( iLess in 1:3) {
		iend <- fullLen
		ibeg <- 1 + iLess
		smlSet <- substr( pepSet, ibeg, iend)
		cat( "\n\nUnique over C-terminal ",(fullLen-iLess), "-mer:  ", length( unique( smlSet)), sep="")
	}
}


peptide.mergeOverlaps <- function( peptides, max.tail=5, min.count=2) {

	if ( is.data.frame( peptides)) {
		if ( ! all( colnames( peptides) == c( "Peptide", "Count"))) 
				stop( "Expected data frame with {Peptide, Count} columns")
		ord <- order( peptides$Count, decreasing=TRUE)
		allPeps <- peptides$Peptide[ord]
		allCnts <- peptides$Count[ord]
	} else {
		allPeps <- peptides
		allCnts <- rep.int( 1, length(peptides))
	}

	cat( "\nGiven Peptides:  ", N <- length( allPeps))

	tapply( 1:N, factor( allPeps), FUN=function(x) {
			if ( length(x) < 2) return()
			newcnt <- sum( allCnts[x])
			allCnts[ x[1]] <<- newcnt
			allCnts[ x[2:length(x)]] <<- 0
			return()
		}, simplify=FALSE)

	peps <- allPeps[ allCnts > 0]
	cnts <- allCnts[ allCnts > 0]
	cat( "\nUnique Peptides: ", length( peps))

	# only keep the ones with at least K reads to support it
	keep <- which( cnts >= min.count)
	peps <- peps[ keep]
	cnts <- cnts[ keep]
	len <- nchar( peps)
	N <- length( peps)
	cat( "\nCounts >= ",min.count,":  ", N)


	# build a set of A-Z 'first AA' sets
	first <- substr( peps, 1,1)
	letterSets <- lapply( LETTERS, function(x) which( first == x))
	names( letterSets) <- LETTERS

	cat( "\n")
	for ( iter in 1:max.tail) {
		nmerge <- 0
		for ( i in 1:N) {

			if ( len[i] < iter) next
			# extract the suffix from one peptide, and see all it hits
			mykey <- substr( peps[i], 1, 1)
			prefix <- substr( peps[i], 1, iter)
			suffix <- substr( peps[i], iter+1, len[i])
			theirkey <- substr( suffix, 1,1)
			totest <- letterSets[[ theirkey]]
			hits <- grep( suffix, peps[totest], fixed=T)
			#targets <- substr( peps[totest], 1, nchar(suffix))
			#hits <- which( targets == suffix)
			if ( length(hits) < 1) next
			hits <- totest[ hits]
			hits <- setdiff( hits, i)
			if ( length(hits) < 1) next

			# for each other peptide that starts with my suffix, prepend my prefix
			for ( j in hits) {
				loc <- regexpr( suffix, peps[j], fixed=TRUE)
				if ( loc != 1) next
				newpep <- paste( prefix, peps[j], sep="")
				letterSets[[ theirkey]] <- setdiff( letterSets[[ theirkey]], j)
				letterSets[[ mykey]] <- c( letterSets[[ mykey]], j)
				peps[j] <- newpep
				len[j] <- nchar(newpep)
				cnts[j] <- cnts[j] + cnts[i]
				nmerge <- nmerge + 1
			}

			# remove me fron the set
			letterSets[[mykey]] <- setdiff( letterSets[[mykey]], i)
			peps[i] <- ""
			len[i] <- 0
			cnts[i] <- 0
			if ( i %% 100 == 0) cat( "\r", i, nmerge)
		}
		cat( "\nIteration: ", iter, "\nTable of peptide lengths:\n")
		print( ans <- table( nchar( peps)))
		cat( "\nTotal:  ", sum( ans), "\n")
	}

	keep <- which( len > 0 & cnts > 0)
	out <- data.frame( "Peptide"=peps[keep], "Count"=cnts[keep], stringsAsFactors=FALSE)
	ord <- order( out$Peptide, -out$Count)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}


`peptides2Proteome` <- function( tbl, geneColumn="GeneID", peptideColumn="Peptide", countColumn="Count",
				minAAperSample=100000) {

	if ( ! ( peptideColumn %in% colnames(tbl))) stop( "Peptide Column Name not found")
	peptides <- tbl[[ peptideColumn]]
	if ( ! ( geneColumn %in% colnames(tbl))) stop( "Gene Column Name not found")
	genes <- tbl[[ geneColumn]]

	if ( is.null(countColumn)) {
		counts <- rep.int( 1, length(peptides))
	} else {
		if ( ! ( countColumn %in% colnames(tbl))) stop( "Count Column Name not found")
		counts <- tbl[[ countColumn]]
	}

	gmap <- subset.data.frame( getCurrentGeneMap(), REAL_G == T)
	mySpecies <- getCurrentSpecies()

	# we will use the length of each peptide
	peplen <- nchar( peptides)
	aaPerSample <- sum(  counts * peplen)
		
	geneFac <- factor( genes)
	genesOut <- levels(geneFac)

	genePtr <- match( genesOut, gmap$GENE_ID, nomatch=0)
	if ( all( genePtr == 0)) warning( "No Peptide GeneIDs match current SpeciesID...")

	totSpecOut <- uniSpecOut <- aapkmOut <- pctOut <- vector()

	if ( aaPerSample < minAAperSample) {
		cat( "\nNot enough peptides in the proteome..")
		cat( "\nN_AA_Found: ", aaPerSample)
		cat( "\nN_AA_Minimu: ", minAAperSample)
		aaPerSample <- minAAperSample
	}

	for ( i in 1:length(genesOut)) {

		thisGene <- genesOut[i]
		who <- which( genes == thisGene)
		uniSpecOut[i] <- length( unique( peptides[who]))
		totSpecOut[i] <- sum( counts[who])
		pctOut[i] <- totSpecOut[i] * 100 / sum(counts)

		myAAcount <- sum( counts[who] * peplen[who])
		if ( genePtr[i] == 0) {
			myProteinLength <- 1000
		} else {
			myProteinLength <- round( gmap$N_EXON_BASES[ genePtr[i]] / 3)
		}
		aapkmOut[i] <- myAAcount * (1000 / myProteinLength) * (1000000 / aaPerSample)
	}

	if ( mySpecies != "Pf3D7") {
		pfGenes <- ortholog( genesOut, from=mySpecies, to="Pf3D7")
		setCurrentSpecies( "Pf3D7")
		origGenes <- gene2OrigID( pfGenes)
		setCurrentSpecies( mySpecies)
		out <- data.frame( "GENE_ID"=genesOut, "PF_ID"=pfGenes, "ORIG_ID"=origGenes, 
				"PRODUCT"=gene2ProductAllSpecies(genesOut), "AAPKM"=aapkmOut, "PCT_SPECTRA"=pctOut,
				"TOTAL_SPECTRA"=totSpecOut, "UNIQUE_SPECTRA"=uniSpecOut, 
				stringsAsFactors=FALSE)
	} else {
		origGenes <- gene2OrigID( genesOut)
		out <- data.frame( "GENE_ID"=genesOut, "ORIG_ID"=origGenes, "PRODUCT"=gene2ProductAllSpecies(genesOut), 
				"AAPKM"=aapkmOut,"PCT_SPECTRA"=pctOut, "TOTAL_SPECTRA"=totSpecOut, 
				"UNIQUE_SPECTRA"=uniSpecOut, stringsAsFactors=FALSE)
	}
	ord <- order( out$AAPKM, decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	out$AAPKM <- round( out$AAPKM, digits=3)
	out$PCT_SPECTRA <- round( out$PCT_SPECTRA, digits=4)

	return( out)
}


`proteomeDiffAbundance` <- function( file1, file2, gene1Column="GENE_ID", gene2Column=gene1Column,
					value1Column="AAPKM", value2Column=value1Column,
					count1Column="TOTAL_SPECTRA", count2Column=count1Column,
					minimumValue=NULL, sep="\t", dropLowCountCutoff=NULL) {

	tbl1 <- read.delim( file1, as.is=T, sep=sep)
	genes1 <- tbl1[[ gene1Column]]
	if ( is.null( genes1)) stop( paste( "GeneID column not found.  Tried: ", gene1Column, "\nFound: ", colnames(tbl1)))
	value1 <- tbl1[[ value1Column]]
	if ( is.null( value1)) stop( paste( "Abundance column not found.  Tried: ", value1Column, "\nFound: ", colnames(tbl1)))
	count1 <- tbl1[[ count1Column]]
	if ( is.null( count1)) stop( paste( "Count column not found.  Tried: ", count1Column, "\nFound: ", colnames(tbl1)))

	tbl2 <- read.delim( file2, as.is=T, sep=sep)
	genes2 <- tbl2[[ gene2Column]]
	if ( is.null( genes2)) stop( paste( "GeneID column not found.  Tried: ", gene2Column, "\nFound: ", colnames(tbl2)))
	value2 <- tbl2[[ value2Column]]
	if ( is.null( value2)) stop( paste( "Abundance column not found.  Tried: ", value2Column, "\nFound: ", colnames(tbl2)))
	count2 <- tbl2[[ count2Column]]
	if ( is.null( count2)) stop( paste( "Count column not found.  Tried: ", count2Column, "\nFound: ", colnames(tbl2)))

	# we will keep the union of both gene sets
	allGenes <- base::sort( unique.default( c( genes1, genes2)))
	wh1 <- base::match( allGenes, genes1, nomatch=0)
	wh2 <- base::match( allGenes, genes2, nomatch=0)
	NG <- length( allGenes)
	allProd <- gene2ProductAllSpecies( allGenes)

	# get the values, and use a minimum to prevent divide by zero
	v1 <- as.numeric( value1)
	v2 <- as.numeric( value2)
	if ( is.null( minimumValue)) {
		minV <- max( min(v1[v1 > 0]), min(v2[v2 > 0]))
		minimumValue <- minV * 2
	}
	cat( "\nAdding small linear offset to prevent divide by zero:  ", minimumValue)
	bigv1 <- bigv2 <- rep( 0, times=NG)
	bigv1[ wh1 > 0] <- v1[wh1]
	bigv2[ wh2 > 0] <- v2[wh2]
	usev1 <- bigv1 + minimumValue
	usev2 <- bigv2 + minimumValue

	# grab the peptide count to show with results
	n1 <- as.numeric( count1)
	n2 <- as.numeric( count2)
	bign1 <- bign2 <- rep( 0, times=NG)
	bign1[ wh1 > 0] <- n1[wh1]
	bign2[ wh2 > 0] <- n2[wh2]
	if ( ! is.null( dropLowCountCutoff)) {
		cat( "\nRemoving proteins with less than", dropLowCountCutoff, "peptide/spectra counts: ")
		bignMax <- pmax( bign1, bign2)
		dropTooLow <- which( bignMax < dropLowCountCutoff)
		cat( " ", length( dropTooLow))
		if ( ! length( dropTooLow)) rm( dropTooLow)
	}

	log2fold <- log2( usev1 / usev2)

	out <- data.frame( "GENE_ID"=allGenes, "PRODUCT"=allProd, "LOG2FOLD"=log2fold, 
				"AAPKM_1"=bigv1, "AAPKM_2"=bigv2,
				"TOTAL_SPECTRA_1"=bign1, "TOTAL_SPECTRA_2"=bign2, 
				stringsAsFactors=T)

	if ( exists( "dropTooLow")) out <- out[ -dropTooLow, ]

	ord <- order( out$LOG2FOLD, decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	# let's give the output the column names we used as input
	colnames(out)[4:5] <- paste( c(value1Column,value2Column), 1:2, sep="_")
	colnames(out)[6:7] <- paste( c(count1Column,count2Column), 1:2, sep="_")

	return( out)
}


molecularWeight <- function( peptides) {

	# wt per AA
	AA_WTS <- c( "A"=71.09, "B"=NA, "C"=103.15, "D"=115.09, "E"=129.12, "F"=147.18, "G"=57.05, "H"=137.14, "I"=113.16, "J"=NA,
		"K"=128.17, "L"=113.16, "M"=131.19, "N"=114.11, "O"=NA, "P"=97.12, "Q"=128.14, "R"=156.19, "S"=87.08, "T"=101.11,
		"U"=168.07, "V"=99.14, "W"=186.21, "X"=NA, "Y"=163.18, "Z"=NA)

	N <- length( peptides)
	terms <- strsplit( peptides, split="")
	ans <- sapply( terms, FUN=function(x) {
			# get how many of each AA, and sum up that many of each weight
			aaTbl <- table( x)
			aaIDs <- names( aaTbl)
			where <- match( aaIDs, names(AA_WTS), nomatch=0)
			good <- which( where > 0)
			myWts <- as.vector( aaTbl[good]) * AA_WTS[where[good]]
			return( sum( myWts, na.rm=T))
		})

	names(ans) <- names(peptides)
	return( round( ans))
}

