# geneModelTools.R -- functions to automate exon splice boundary detection for ortholog proteins

require(Biostrings)
require(pwalign)
data(BLOSUM62)

# error flag for any invalid function return value
ERROR <- "ERROR"


# one top level wrapper...?
# optimizeGeneExonSpliceBoundaries <- function( starts, stops, gdna, idealAA) {
#}


# optimizer functions for each class of issue

bestPhaseShift <- function( starts, stops, strand, gdna, idealAA) {

	# jitter each exon back and forth, to see if phase shift of any one exon helps
	
	# see what no phase shift at all score is
	bestAA <- joinExonsToAA( starts, stops, strand, gdna)
	bestScore <- aaSimilarityScore( bestAA, idealAA, type="global")
	bestDesc <- "No change"
	bestStarts <- starts
	bestStops <- stops
	nex <- length(starts)
	
	for ( i in 1:nex) {
		startsNow <- starts
		stopsNow <- stops
		for (jit in c(-1, -2, 1, 2)) {
			startsNow[i] <- starts[i] + jit
			stopsNow[i] <- stops[i] + jit
			cdnaNow <- joinExonsToCDNA( startsNow, stopsNow, gdna)
			if (cdnaNow == ERROR) next
			aaNow <- translateToAA( cdnaNow, strand)
			scoreNow <- aaSimilarityScore( aaNow, idealAA, type="global")
			if ( scoreNow > bestScore) {
				bestDesc <- paste( "Exon", i, "phase shift of", formatC( jit, format="d", flag="+"))
				bestScore <- scoreNow
				bestStarts <- startsNow
				bestStops <- stopsNow
				bestAA <- aaNow
			}
		}
	}
	# all phase test shifts done.  Which was the best?
	out <- list( "starts"=bestStarts, "stops"=bestStops, "description"=bestDesc, "score"=bestScore, "aa"=bestAA)
	return( out)
}


bestSpliceSiteShift <- function( starts, stops, strand, gdna, idealAA) {

	# jitter each exon pair back and forth, to see if an in-phase shift of the splice site helps
	
	# see what no phase shift at all score is
	bestAA <- joinExonsToAA( starts, stops, strand, gdna)
	bestScore <- aaSimilarityScore( bestAA, idealAA, type="global")
	bestDesc <- "No change"
	bestStarts <- starts
	bestStops <- stops
	nsplice <- length(starts) - 1
	if (nsplice < 1) return( list( "starts"=bestStarts, "stops"=bestStops, "description"=bestDesc, "score"=bestScore, "aa"=bestAA))
	
	for ( i in 1:nsplice) {
		startsNow <- starts
		stopsNow <- stops
		for (jit in c(-1, -2, 1, 2)) {
			stopsNow[i] <- stops[i] + jit
			startsNow[i+1] <- starts[i+1] + jit
			cdnaNow <- joinExonsToCDNA( startsNow, stopsNow, gdna)
			if (cdnaNow == ERROR) next
			aaNow <- translateToAA( cdnaNow, strand)
			scoreNow <- aaSimilarityScore( aaNow, idealAA, type="global")
			if ( scoreNow > bestScore) {
				bestDesc <- paste( "Splice", i, "shift of", formatC( jit, format="d", flag="+"))
				bestScore <- scoreNow
				bestStarts <- startsNow
				bestStops <- stopsNow
				bestAA <- aaNow
			}
		}
	}
	# all phase test shifts done.  Which was the best?
	out <- list( "starts"=bestStarts, "stops"=bestStops, "description"=bestDesc, "score"=bestScore, "aa"=bestAA)
	return( out)
}


bestExonLength <- function( starts, stops, strand, gdna, idealAA) {

	# grow/shrink each exon, to see if elongation or truncation helps
	
	# see what no phase shift at all score is
	bestAA <- joinExonsToAA( starts, stops, strand, gdna)
	bestScoreG <- aaSimilarityScore( bestAA, idealAA, type="global")
	bestDesc <- "No change"
	bestStarts <- starts
	bestStops <- stops
	nex <- length(starts)
	# also track scores as per AA, to make sure any change is truly better
	bestScoreL <- aaSimilarityScore( bestAA, idealAA, type="local")
	bestScoreGperAA <- bestScoreG / nchar(bestAA)
	bestScoreLperAA <- bestScoreL / nchar(bestAA)
	
	# use a fibonacci series to try various lengths, in AA units
	fibo <- c( 1, 2, 3, 5, 8, 13, 21)
	for ( i in 1:nex) {
		startsNow <- starts
		stopsNow <- stops
		for (delta in c( fibo, -fibo)) {
			startsNow[i] <- starts[i] + delta*3
			cdnaNow <- joinExonsToCDNA( startsNow, stopsNow, gdna)
			if (cdnaNow == ERROR) next
			aaNow <- translateToAA( cdnaNow, strand)
			scoreNowG <- aaSimilarityScore( aaNow, idealAA, type="global")
			if ( scoreNowG > bestScoreG) {
				# do other tests to make sure we keep a truly better choice
				scoreNowGperAA <- scoreNowG / nchar( aaNow)
				scoreNowL <- aaSimilarityScore( aaNow, idealAA, type="local")
				scoreNowLperAA <- scoreNowL / nchar( aaNow)
				if ( scoreNowL > bestScoreL && scoreNowGperAA > bestScoreGperAA && scoreNowLperAA > bestScoreLperAA ) {
					bestDesc <- paste( "Exon", i, "shift start by", formatC( delta, format="d", flag="+"), "aa")
					bestScoreG <- scoreNowG
					bestScoreL <- scoreNowL
					bestStarts <- startsNow
					bestStops <- stopsNow
					bestAA <- aaNow
					bestScoreGperAA <- bestScoreG / nchar(bestAA)
					bestScoreLperAA <- bestScoreL / nchar(bestAA)
				}
			}
			startsNow[i] <- starts[i]
			stopsNow[i] <- stops[i] + delta*3
			cdnaNow <- joinExonsToCDNA( startsNow, stopsNow, gdna)
			if (cdnaNow == ERROR) next
			aaNow <- translateToAA( cdnaNow, strand)
			scoreNowG <- aaSimilarityScore( aaNow, idealAA, type="global")
			if ( scoreNowG > bestScoreG) {
				# do other tests to make sure we keep a truly better choice
				scoreNowGperAA <- scoreNowG / nchar( aaNow)
				scoreNowL <- aaSimilarityScore( aaNow, idealAA, type="local")
				scoreNowLperAA <- scoreNowL / nchar( aaNow)
				if ( scoreNowL > bestScoreL && scoreNowGperAA > bestScoreGperAA && scoreNowLperAA > bestScoreLperAA ) {
					bestDesc <- paste( "Exon", i, "shift stop by", formatC( delta, format="d", flag="+"), "aa")
					bestScoreG <- scoreNowG
					bestScoreL <- scoreNowL
					bestStarts <- startsNow
					bestStops <- stopsNow
					bestAA <- aaNow
					bestScoreGperAA <- bestScoreG / nchar(bestAA)
					bestScoreLperAA <- bestScoreL / nchar(bestAA)
				}
			}
		}
	}
	# all phase test shifts done.  Which was the best?
	out <- list( "starts"=bestStarts, "stops"=bestStops, "description"=bestDesc, "score"=bestScoreG, "aa"=bestAA)
	return( out)
}


# visualizers, etc

aaAlignment <- function( aa, idealAA, line.width=100, type=c("global","local")) {
	
	type <- match.arg( type)
	out <- pairwiseAlignment( aa, idealAA, type=type, scoreOnly=FALSE, substitutionMatrix=BLOSUM62,
								gapOpening = -8, gapExtension = -2)
	return( writePairwiseAlignments( out, block.width=line.width, Matrix="BLOSUM62"))	
}

asGeneModel <- function( starts, stops, strand, gdna, line.width=100) {

	# write out the best model we found, with some hints about exon/CDS junctions
	nex <- length(starts)
	cdna <- substr( rep.int(gdna,nex), starts, stops)
	if ( strand[1] == "+") {
		cdnaWithJunct <- paste( cdna, collapse="|")
		cdnaNoJunct <- paste( cdna, collapse="")
	} else {
		cdna <- sapply( cdna, myReverseComplement)
		cdna <- rev( cdna)
		cdnaWithJunct <- paste( cdna, collapse="|")
		cdnaNoJunct <- paste( cdna, collapse="")
	}
	aaWithJunct <- aaNoJunct <- DNAtoAA( cdnaNoJunct, clipAtStop=FALSE, readingFrames=1)
	# add exon junction marker hints to the AA too
	if ( nex > 1) {
		whereJunct <- gregexpr( "|", cdnaWithJunct, fixed=T)[[1]]
		aaGapPts <- floor( as.numeric( whereJunct) / 3)
		aaWithJunct <- aaNoJunct
		offset <- 0
		for ( k in aaGapPts) {
			aaWithJunct <- paste( substr( aaWithJunct, 1, k+offset), substring( aaWithJunct, k+offset+1), sep="|")
			offset <- offset + 1
		}
	}
	# now write them out 
	outTxt1 <- vector()
	nOut <- 0
	ncOut <- 0
	nbp <- nchar( cdnaWithJunct)
	while (ncOut < nbp) {
		now <- min( nbp, ncOut+line.width)
		nOut <- nOut + 1
		outTxt1[nOut] <- substr( cdnaWithJunct, (ncOut+1), now)
		ncOut <- now
	}
	outTxt2 <- vector()
	nOut <- 0
	ncOut <- 0
	naa <- nchar( aaWithJunct)
	while (ncOut < naa) {
		now <- min( naa, ncOut+line.width)
		nOut <- nOut + 1
		outTxt2[nOut] <- substr( aaWithJunct, (ncOut+1), now)
		ncOut <- now
	}
	return( list( "cdna"=outTxt1, "aa"=outTxt2))
}


# core low level extractor/translation/scoring functions
joinExonsToCDNA <- function( starts, stops, gdna) {

	nex <- length(starts)
	ngdna <- nchar( gdna)[1]
	if ( length(stops) != nex) return(ERROR)
	if ( any( starts < 1)) return(ERROR)
	if ( any( starts > stops)) return(ERROR)
	if ( any( stops > ngdna)) return(ERROR)
	if ( any( diff( starts) < 1)) return(ERROR)
	if ( any( diff( stops) < 1)) return(ERROR)
	cdna <- substr( rep.int(gdna,nex), starts, stops)
	out <- paste( cdna, collapse="")
	return( out)
}

joinExonsToAA <- function( starts, stops, strand, gdna) {

	cdna <- joinExonsToCDNA( starts, stops, gdna)
	if (cdna == ERROR) return(ERROR)
	aa <- translateToAA( cdna, strand)
	return( aa)
}

translateToAA <- function( cdna, strand) {
	
	if ( cdna == ERROR) return(ERROR)
	out <- DNAtoAA( cdna, clipAtStop=FALSE, readingFrames=if (strand[1] == "+") 1 else 4)
	return(out)
}

aaSimilarityScore <- function( aa, idealAA, type=c("global","local")) {

	if ( aa == ERROR) return(0)
	type <- match.arg( type)
	out <- pairwiseAlignment( aa, idealAA, type=type, scoreOnly=TRUE, substitutionMatrix=BLOSUM62,
								gapOpening = -8, gapExtension = -2)
	return(out)	
}
