# codonTools.R

APPEND <- base::append
GREP <- base::grep
GSUB <- base::gsub
MAPPLY <- base::mapply
MATCH <- base::match
ORDER <- base::order
PASTE <- base::paste
REV <- base::rev
SAPPLY <- base::sapply
STRSPLIT <- base::strsplit
SUB <- base::sub
SUBSTR <- base::substr
TAPPLY <- base::tapply
TOUPPER <- base::toupper
WHICH <- base::which


`codon.defaults` <- function() {
	# the data object is called 'codonMap'
	data( CodonMap, envir=environment())
	CodonEnv[[ "CodonMap"]] <- codonMap
	if ( nrow( codonMap) < 64) warning("CodonMap:  bad amino acid lookup table")
}


`getCodonMap` <- function() return( CodonEnv[[ "CodonMap"]])


codonUsageFrequency <- function( dna, allAA=FALSE) {

	dnaV <- STRSPLIT( dna, split="")[[1]]
	Ndna <- length( dnaV)
	dnaV <- TOUPPER( dnaV)

	# we only want the true codons, none of the extended
	codonMap <- getCodonMap()[ 1:64, ]

	beg <- seq.int( 1, Ndna, 3)
	Naa <- length(beg)
	end <- beg + 2
	if (end[Naa] > Ndna) {
		Naa <- Naa - 1
		length(beg) <- length(end) <- Naa
	}
	aa <- rep.int( "?", Naa)

	triplets <- MAPPLY( beg, end, FUN=function(b,e) { return( PASTE( dnaV[b:e], collapse=""))},
			USE.NAMES=FALSE)
	where <- MATCH( triplets, codonMap$DNA, nomatch=0)
	aa[ where > 0] <- codonMap$AA[ where]
	out <- PASTE( aa, triplets, sep="_")

	# we may be asked to make sure every AA is seen at least once
	if (allAA) {
		missingA <- base::setdiff( codonMap$AA, aa)
		if ( length(missingA)) {
			smlMap <- subset.data.frame( codonMap, AA %in% missingA)
			out2 <- PASTE( smlMap$AA, smlMap$DNA, sep="_")
			out <- c( out, out2)
		}
	}

	return( table(out))
}


# turn amino acids to DNA
`AAtoCodonOptimizedDNA` <- function( aa, dnaRef) {

	codonFreq <- codonUsageFrequency( dnaRef, allAA=TRUE)
	out <- rep.int( "", length(aa))
	
	aa <- TOUPPER( aa)
	aaV <- STRSPLIT( aa, split="")

	# test to provide a clean error
	allDefinedAA <- unique.default( SUBSTR( names(codonFreq), 1,1))
	allGivenAA <- unique.default( unlist( aaV))
	missing <- base::setdiff( allGivenAA, allDefinedAA)
	if ( length( missing)) {
		cat( "\n\nError:  given some undefined AA characters: ", missing)
		stop( "Unable to convert AA to Codon Optimized DNA.")
	}

	for ( j in 1:length(aa)) {
		thisAAvec <- aaV[[j]]
		aaPattern <- PASTE( "^",thisAAvec, sep="")
		aaPattern[ thisAAvec == "*"] <- "^\\*"
		dnaV <- SAPPLY( 1:length(thisAAvec), function(i) {
				who <- GREP( aaPattern[i], names(codonFreq))
				# use the codon frequencies as the sampling probabilities
				prob <- codonFreq[who]
				dnaSet <- SUB( "^[A-Z\\*]_", "", names(codonFreq)[who])
				return( sample( dnaSet, size=1, prob=prob))
			})
		out[j] <- PASTE( dnaV, collapse="")
	}
	return( out)
}


# turn DNA to AA
DNAtoAA <- function( dna, clipAtStop=TRUE, readingFrames=1:6) {

	dnaV <- STRSPLIT( dna, split="")[[1]]
	nc <- length( dnaV)
	readingFrames <- base::intersect( 1:6, readingFrames)
	if ( any( readingFrames %in% 4:6)) {
		dnaRV <- myReverseComplement( dna, as.vector=TRUE)
	}
	nFrames <- length(readingFrames)
	out <- rep( "", times=nFrames)
	names( out) <- readingFrames
	if ( nc < 3) return( out)


	for ( i in 1:length(readingFrames)) {
		k <- readingFrames[i]
		if ( k <= 3) {
			aaV <- DNAVtoAAV( dnaV[ k:nc])
		} else {
			aaV <- DNAVtoAAV( dnaRV[ (k-3):nc])
		}

		if ( clipAtStop) {
			where <- MATCH( STOP_CODON, aaV, nomatch=0)
			if ( where > 0) length(aaV) <- where
		}
		out[i] <- PASTE( aaV, collapse="")
	}
	out
}


DNAVtoAAV <- function( dnaVec) {


	# faster simple case:  'dna' is a vector, in frame, forward strand.
	Ndna <- length(dnaVec)
	if ( Ndna < 3) return("")

	dna <- TOUPPER( dnaVec)
	codonMap <- getCodonMap()

	beg <- seq.int( 1, Ndna, 3)
	Naa <- length(beg)
	end <- beg + 2
	if (end[Naa] > Ndna) {
		Naa <- Naa - 1
		length(beg) <- length(end) <- Naa
	}
	out <- rep.int( "?", Naa)

	triplets <- MAPPLY( beg, end, FUN=function(b,e) PASTE( dna[b:e], collapse=""),
			USE.NAMES=FALSE)
	where <- MATCH( triplets, codonMap$DNA, nomatch=0)
	out[ where > 0] <- codonMap$AA[ where]
	out
}


myReverseComplement <- function( dna, as.vector=FALSE) {

	if ( length(dna) > 1) {
		warning( "'reverseComplement' requires a single charater string...dropping some.")
		dna <- dna[1]
	}
	nc <- base::nchar(dna)
	dnaV <- base::unlist( STRSPLIT( dna, split=""))
	newV <- REV( dnaV)
	back <- nc:1
	newV[ dnaV[ back] == "A"] <- "T"
	newV[ dnaV[ back] == "C"] <- "G"
	newV[ dnaV[ back] == "G"] <- "C"
	newV[ dnaV[ back] == "T"] <- "A"
	newV[ dnaV[ back] == "N"] <- "N"
	newV[ dnaV[ back] == "a"] <- "t"
	newV[ dnaV[ back] == "c"] <- "g"
	newV[ dnaV[ back] == "g"] <- "c"
	newV[ dnaV[ back] == "t"] <- "a"
	newV[ dnaV[ back] == "n"] <- "n"
	if (as.vector) return( newV)
	return( PASTE( newV, collapse=""))
}


myReverse <- function( dna, as.vector=FALSE) {

	if ( length(dna) > 1) {
		warning( "'reverse' requires a single charater string...dropping some.")
		dna <- dna[1]
	}
	dnaV <- base::unlist( STRSPLIT( dna, split=""))
	newV <- REV( dnaV)
	if (as.vector) return( newV)
	return( PASTE( newV, collapse=""))
}


fastAAreadingFrame <- function( peptideTriple, protein, AAprotein, max.mismatch=3) {

	# can we find the best reading frame, given a 'target' protein?
	require( Biostrings)

	# try a fast way first, see if one exact match
	hits <- vector( length=3)
	hits[1] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[2] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[3] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	if ( sum( hits > 0) == 1) {
		best <- which.max( hits)
		bestPep <- peptideTriple[best]
		bestFirstAA <- hits[best]
		return( list( "Position"=bestFirstAA, "Length"=base::nchar(bestPep),
			"Mismatches"=0))
	}

	# see how long the AA chain is...
	firstStopCodon <- regexpr( STOP_CODON, peptideTriple, fixed=TRUE)
	lenAA <- ifelse( firstStopCodon > 0, firstStopCodon - 1, base::nchar( peptideTriple))

	# count up the size of the matches, with various mismatch threshold...
	nMismatch <- 0
	repeat {
		scores <- firstAAinProtein <- lenMatch <- nMisMatch <- rep( 0, times=3)
		for (j in 1:3) {
			pep <- peptideTriple[j]
			len <- lenAA[j]
			if ( len > 0) {
				mp <- matchPattern( pep, AAprotein, max.mismatch=nMismatch)
				if (length(mp) > 0) {
					lens <- nchar( mp)
					bestView <- which.max( lens)
					bestLen <- lens[bestView]
					scores[j] <- bestLen / width(mp)[bestView]
					firstAAinProtein[j] <- start(mp)[bestView]
					lenMatch[j] <- bestLen
					nMisMatch[j] <- nMismatch
				}
			}
		}
		# perfect result is one peptide exactly 1.0, and two at 0.0
		who1 <- (scores == 1.0)
		if ( sum(who1) == 1) {
			curBest <- which.max( scores)
			break
		}
		# not quite perfect,
		curBest <- which.max( scores)
		nMismatch <- nMismatch + 1
		if (nMismatch > max.mismatch) break
	}
	bestPep <- peptideTriple[curBest]
	bestFirstAA <- firstAAinProtein[ curBest]
	return( list( "Position"=bestFirstAA, "Length"=lenMatch[ curBest], 
			"Mismatches"=nMisMatch[curBest]))
}


bestAAreadingFrame <- function( peptideTriple, protein, max.mismatch=3) {

	# can we find the best reading frame, given a 'target' protein?
	require( Biostrings)

	# try a fast way first, see if one exact match
	hits <- vector( length=3)
	hits[1] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[2] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	hits[3] <- regexpr( peptideTriple[1], protein, fixed=TRUE)
	if ( sum( hits > 0) == 1) {
		best <- which.max( hits)
		bestPep <- peptideTriple[best]
		bestFirstAA <- hits[best]
		return( list( "Peptide"=bestPep, "Position"=bestFirstAA, "Length"=base::nchar(bestPep),
			"Mismatches"=0))
	}

	# see how long the AA chain is...
	firstStopCodon <- regexpr( STOP_CODON, peptideTriple, fixed=TRUE)
	lenAA <- ifelse( firstStopCodon > 0, firstStopCodon - 1, base::nchar( peptideTriple))

	# count up the size of the matches, with various mismatch threshold...
	nMismatch <- 0
	repeat {
		scores <- rep( 0, times=3)
		firstAAinProtein <- rep( 0, times=3)
		lenMatch <- rep( 0, times=3)
		nMisMatch <- rep( 0, times=3)
		for (j in 1:3) {
			pep <- peptideTriple[j]
			len <- lenAA[j]
			if ( len > 0) {
				mp <<- matchPattern( pep, protein, max.mismatch=nMismatch)
				if (length(mp) > 0) {
					bestView <- which.max( nchar( mp))
					scores[j] <- nchar(mp)[bestView] / width(mp)[bestView]
					firstAAinProtein[j] <- start(mp)[bestView]
					lenMatch[j] <- nchar(mp)[bestView]
					nMisMatch[j] <- nMismatch
				}
			}
		}
		# perfect result is one peptide exactly 1.0, and two at 0.0
		who1 <- (scores == 1.0)
		if ( sum(who1) == 1) {
			curBest <- which.max( scores)
			break
		}
		# not quite perfect,
		curBest <- which.max( scores)
		nMismatch <- nMismatch + 1
		if (nMismatch > max.mismatch) break
	}
	bestPep <- peptideTriple[curBest]
	bestFirstAA <- firstAAinProtein[ curBest]
	#cat( "\nScores: ", scores, "\tBest ReadingFrame: \t", bestPep, "\tStart location:\t", bestFirstAA)
	return( list( "Peptide"=bestPep, "Position"=bestFirstAA, "Length"=lenMatch[ curBest], 
			"Mismatches"=nMisMatch[curBest]))
}


# try to convert DNA to AA given the current annotation info
`convertGenomicBasesToCodingAA` <- function( seqID, position, end, strand="+", dnaQuery, genomeDNA,
			geneMap=NULL, cdsMap=NULL, geneID=NULL) {

	# allow strings as well as expected vectors of bases
	if ( length(dnaQuery) == 1 && nchar(dnaQuery[1]) > 1) {
		dnaQuery <- STRSPLIT( dnaQuery, split="")[[1]]
	}
	if ( length(genomeDNA) == 1 && nchar(genomeDNA[1]) > 1) {
		genomeDNA <- STRSPLIT( genomeDNA, split="")[[1]]
	}

	# given a vector of DNA bases and the its genomic location  < dnaQuery, position, end >,
	# convert that to the amino acid sequence for the coding strand that covers...
	nBases <- length( dnaQuery)
	outG <- outQ <- rep.int( "", nBases)
	aa.ordinal <- rep.int( NA, nBases)
	out <- list( "genomic"=outG, "query"=outQ, "aa.ordinal"=aa.ordinal)

	if ( nBases != (end - position + 1)) {
		warning( "bad DNA query string and/or size")
		return(out)
	}

	# get the gene(s) covered
	if ( is.null( geneMap)) geneMap <- getCurrentGeneMap()
	gmap <- subset.data.frame( geneMap, SEQ_ID == seqID & POSITION < end & END > position & 
				REAL_G == TRUE & STRAND == strand)
	if ( ! is.null(geneID)) gmap <- subset.data.frame( gmap, GENE_ID == geneID)
	if ( nrow(gmap) < 1) return( out)

	# using the CDS instead!!
	if ( is.null( cdsMap)) cdsMap <- getCurrentCdsMap()
	cdsMap <- subset.data.frame( cdsMap, GENE_ID %in% gmap$GENE_ID)

	for( ig in 1:nrow( gmap)) {
		thisG <- gmap$GENE_ID[ ig]
		thisStrand <- gmap$STRAND[ ig]
		cdsmap <- subset.data.frame( cdsMap, GENE_ID == thisG)
		if ( nrow(cdsmap) < 1) next

		# build the string of coding nucleotides from the forward strand
		forwardDNA <- vector()
		for ( ie in 1:nrow(cdsmap)) {
			thisCDS <- genomeDNA[ cdsmap$POSITION[ie] : cdsmap$END[ie] ]
			names( thisCDS) <- cdsmap$POSITION[ie] : cdsmap$END[ie] 
			forwardDNA <- APPEND( forwardDNA, thisCDS)
		}
		codingDNA <- forwardDNA
		if ( thisStrand == "-") {
			tmp <- PASTE( codingDNA, collapse="")
			tmp <- myReverseComplement(tmp)
			tmp <- STRSPLIT( tmp, split="")[[1]]
			codingDNA <- tmp
			names( codingDNA) <- REV( names( forwardDNA))
		}

		# force a trim to multiple of 3 bases
		bigN <- floor( length(codingDNA)/3) * 3
		nCodingBases <- length( codingDNA) <- bigN
		if ( nCodingBases < 3) next

		# convert to AA, ( we know we are in Frame, so just keep the first reading frame)
		thisAA <- DNAVtoAAV( codingDNA)

		# build this back into a vector of single letters, with the AA at the center of its 3 bases
		newstr <- rep.int( "", nCodingBases)
		newAAordinal <- rep.int( NA, nCodingBases)
		newstr[ seq( 2, nCodingBases, by=3)] <- thisAA
		names( newstr) <- names( codingDNA)
		newAAordinal[ seq( 2, nCodingBases, by=3)] <- 1:length(thisAA)
		names( newAAordinal) <- names( codingDNA)

		# lastly, we can put these AA back into the correct location of the output genomic string of text
		for ( ic in 1:nCodingBases) {
			if ( newstr[ic] == "") next
			genLocation <- as.integer( names( newstr)[ic])
			outLocation <- genLocation - position + 1
			if ( outLocation < 1 || outLocation > nBases) next
			outG[ outLocation] <- newstr[ic]
			aa.ordinal[ outLocation] <- newAAordinal[ic]
		}

		# now put the query bases in, and do it again
		outQ <- outG
		codingDNA <- forwardDNA
		codingLocs <- as.integer( names( codingDNA))
		queryLocs <- position : end
		where <- MATCH( codingLocs, queryLocs, nomatch=0)
		codingHits <- WHICH( where > 0)
		if ( length( codingHits) > 0) {
			newstr2 <- newstr
			codingDNA[ codingHits] <- dnaQuery[ where]
			# with the possibility of Indels, these query bases may not be singletons any more
			isDiff <- WHICH( codingDNA != forwardDNA)
			if ( length(isDiff)) {
				baseLen <- nchar( codingDNA)
				# try to estimate how many 'post-indel' AA for each 'pre-indel' codon, pad the ends for edge cases
				baseCumSum <- cumsum( c( 1, baseLen, 1))
				snpDNA <- PASTE( codingDNA, collapse="")
				if ( thisStrand == "-") {
					snpDNA <- myReverseComplement(snpDNA)
					baseCumSum <- cumsum( c( 1, REV(baseLen), 1))
				}
				snpDNA <- STRSPLIT( snpDNA, split="")[[1]]
				if ( thisStrand == "-") {
					names( snpDNA) <- rep( REV( names( forwardDNA)), times=baseLen)
				} else {
					names( snpDNA) <- rep( names( forwardDNA), times=baseLen)
				}
				snpAA <- DNAVtoAAV( snpDNA)
				# this length may be different due to indels, step along by hand...
				newNbases <- min( nCodingBases, length(snpDNA))
				newNaa <- floor( newNbases/3)
				nAAused <- 0
				for ( j in seq( 2, newNbases, by=3)) {
					# account for the padding around the cum sum vector
					ndna <- base::diff( baseCumSum[ c(j-1,j+2)])
					naa <- round( ndna/3)
					if (naa == 1) {
						newstr2[ j] <- snpAA[nAAused+1]
						nAAused <- nAAused + 1
					} else if ( naa < 1) {
						newstr2[ j] <- ""
					} else {
						newstr2[ j] <- PASTE( snpAA[ (nAAused+1):(nAAused+naa)], collapse="")
						if ( thisStrand == "-") newstr2[j] <- myReverse( newstr2[j])
						nAAused <- nAAused + naa
					}
				}
				whodiff <- WHICH( newstr2 != newstr)
				for ( ic in whodiff) {
					if ( newstr2[ic] == "") next
					genLocation <- as.integer( names( newstr2)[ic])
					outLocation <- genLocation - position + 1
					if ( outLocation < 1 || outLocation > nBases) next
					outQ[ outLocation] <- newstr2[ic]
				}
			} else {
				# no Indels, so do it fast and like before...
				if ( thisStrand == "-") {
					tmp <- PASTE( codingDNA, collapse="")
					tmp <- myReverseComplement(tmp)
					tmp <- STRSPLIT( tmp, split="")[[1]]
					codingDNA <- tmp
					names( codingDNA) <- REV( names( forwardDNA))
				}
				thisAA <- DNAVtoAAV( codingDNA)
				newstr2[ seq( 2, nCodingBases, by=3)] <- thisAA
				whodiff <- WHICH( newstr2 != newstr)
				for ( ic in whodiff) {
					if ( newstr2[ic] == "") next
					genLocation <- as.integer( names( newstr2)[ic])
					outLocation <- genLocation - position + 1
					if ( outLocation < 1 || outLocation > nBases) next
					outQ[ outLocation] <- newstr2[ic]
				}
			}
		}
	}

	out <- list( "genomic"=outG, "query"=outQ, "aa.ordinal"=aa.ordinal)
	return( out)
}


# try to convert AA to DNA given the current annotation info
`convertAApositionToGenomicDNAposition` <- function( geneID, AAposition, AAlength=1) {

	outSID <- outPos <- outEnd <- NA
	out <- list( "SEQ_ID"=outSID, "SEQ_POSITION"=outPos, "SEQ_END"=outEnd)

	if ( any( c( length(geneID), length(AAposition), length( AAlength)) > 1)) {
		warning( "Only using the first geneID/position elements")
		geneID <- geneID[1]
		AAposition <- AAposition[1]
		AAlength <- AAlength[1]
	}

	cdsmap <- subset.data.frame( getCurrentCdsMap(), GENE_ID == geneID)
	if ( nrow( cdsmap) < 1) return( out)
	outSID <- cdsmap$SEQ_ID[1]

	# make a little band of relative DNA positions
	vNow <- vector()
	for ( j in 1:nrow(cdsmap)) {
		thisBeg <- cdsmap$POSITION[j]
		thisEnd <- cdsmap$END[j]
		sml <- thisBeg:thisEnd
		vNow <- APPEND( vNow, sml)
	}
	names( vNow) <- 1:length(vNow)
	if ( cdsmap$STRAND[1] == "-") names( vNow) <- REV( 1:length(vNow))

	firstAAbase <- (AAposition-1) * 3 + 1
	lastAAbase <- (AAposition+AAlength-1) * 3 

	where <- MATCH( firstAAbase, names(vNow), nomatch=0)
	if ( where > 0) outPos <- vNow[ where]
	where <- MATCH( lastAAbase, names(vNow), nomatch=0)
	if ( where > 0) outEnd <- vNow[ where]
	if ( outPos > outEnd) { tmp <- outPos; outPos <- outEnd; outEnd <- tmp }
	names(outSID) <- names(outPos) <- names(outEnd) <- ""

	out <- list( "SEQ_ID"=outSID, "SEQ_POSITION"=outPos, "SEQ_END"=outEnd)
	return( out)
}


# try to convert DNA to AA given the current annotation info
`convertGenomicDNApositionToAAposition` <- function( seqID, DNAposition) {

	outPos <- outGene <- NA
	out <- list( "GENE_ID"=outGene, "AA_POSITION"=outPos)

	if ( any( c( length(seqID), length(DNAposition)) > 1)) {
		warning( "Only using the first seqID/position elements")
		seqID <- seqID[1]
		DNAposition <- DNAposition[1]
	}

	gmap <- subset.data.frame( getCurrentGeneMap(), SEQ_ID == seqID)
	if ( nrow( gmap) < 1) return( out)
	where <- findInterval( DNAposition, gmap$POSITION)
	if ( where < 1 || where > nrow(gmap)) return(out)
	outGene <- geneID <- gmap$GENE_ID[ where[1]]

	cdsmap <- subset.data.frame( getCurrentCdsMap(), GENE_ID == geneID)
	if ( nrow( cdsmap) < 1) return( out)

	# make a little band of relative DNA positions
	vNow <- vector()
	for ( j in 1:nrow(cdsmap)) {
		thisBeg <- cdsmap$POSITION[j]
		thisEnd <- cdsmap$END[j]
		sml <- thisBeg:thisEnd
		vNow <- APPEND( vNow, sml)
	}
	names( vNow) <- 1:length(vNow)
	if ( cdsmap$STRAND[1] == "-") names( vNow) <- REV( 1:length(vNow))

	where <- base::match( DNAposition, vNow, nomatch=0)
	if ( where > 0) {
		localDNApos <- as.numeric( names( vNow)[where])
		outPos <- floor( (localDNApos-1) / 3) + 1
	}

	out <- list( "GENE_ID"=outGene, "AA_POSITION"=outPos)
	return( out)
}


# try to convert genomic DNA to coding DNA given the current annotation info
`convertGenomicDNAtoCodingDNA` <- function( geneID, genomicDNA=NULL) {

	outDNA <- ""

	cdsmap <- subset.data.frame( getCurrentCdsMap(), GENE_ID == geneID)
	# CDS annotation is primary...  but if not there, see if gene is in exon map as a backup
	if ( nrow( cdsmap) < 1) {
		cdsmap <- subset.data.frame( getCurrentExonMap(), GENE_ID == geneID)
		if ( nrow( cdsmap) < 1) return( outDNA)
	}
	if ( is.null( genomicDNA)) {
		cat( "\nRequired 'genomicDNA' argument is missing...")
		return( outDNA)
	}
	seqID <- cdsmap$SEQ_ID[1]

	# make a little band of absolute DNA positions
	vNow <- vector()
	for ( j in 1:nrow(cdsmap)) {
		thisBeg <- cdsmap$POSITION[j]
		thisEnd <- cdsmap$END[j]
		sml <- thisBeg:thisEnd
		vNow <- APPEND( vNow, sml)
	}

	# get that chunk of genomic DNA
	beg <- min( vNow)
	end <- max( vNow)
	myDNA <- STRSPLIT( SUBSTR( genomicDNA, beg, end), split="")[[1]]

	# convert to relative positions
	vNow <- vNow - beg + 1
	myDNA <- myDNA[ vNow]
	outDNA <- PASTE( myDNA, collapse="")

	# flip if reverse strand
	if ( cdsmap$STRAND[1] == "-") outDNA <- myReverseComplement( outDNA)

	return( outDNA)
}


`buildCodonFreqMap` <- function( AAfasta, speciesID="Hs_grc", fraction=1.0, genomicFastaFile=NULL) {

	codonTable <- measureCodonFreq( AAfasta=AAfasta, speciesID=speciesID, fraction=fraction, 
					genomicFastaFile=genomicFastaFile)
	codonMap <- calculateCodonScores( codonTable)

	return( list( "table"=codonTable, "map"=codonMap))
}


`measureCodonFreq` <- function( AAfasta, speciesID="Hs_grc", fraction=1.0, 
				GImap=NULL, genomicFastaFile=NULL) {

	setCurrentSpecies( speciesID)
	geneMap <- getCurrentGeneMap()
	geneMap <- subset.data.frame( geneMap, REAL_G == TRUE)

	# given a pair of matching fasta files, one amino acids and one cDNA,
	# build the table of most frequent codons for each AA

	aaList <- loadFasta( AAfasta)
	aaNames <- aaList$desc
	aaNames <- SUB( "(gi.+ref\\|)(NP_.+)(\\.[1-9]\\|.+$)", "\\2", aaNames)
	aaSeqs <- aaList$seq
	cat( "\n  N_Proteins:  ", length( aaNames), head( aaNames))

	# convert to GeneIDs now...
	if ( is.null( GImap)) {
		aaGeneIDs <- aaNames
	} else {
		where <- MATCH( aaNames, GImap$protAcc, nomatch=0)
		aaGeneIDs <- rep( "", times=length(where))
		thisGeneName <- GImap$name[where]
		thisGInumber <- GImap$locusLinkId[where]
		allGeneID <- PASTE( thisGeneName, thisGInumber, sep=":GI")
		aaGeneIDs[ where > 0] <- allGeneID
	}
	aaGeneIDs <- intersect( aaGeneIDs, geneMap$GENE_ID)
	cat( "\n  N_Matching_Genes:   ", length( aaGeneIDs), head( aaGeneIDs))

	cat( "\nLooking up genes in geneMap...\n")
	aaGenePtrs <- SAPPLY( aaGeneIDs, function(x) {
				if ( x == "") return( 0)
				where <- GREP( x, geneMap$GENE_ID, fixed=T)[1]
				if ( is.na(where)) return( 0)
				return( where)
			})
	cat( "\nN_found:  ", sum( aaGenePtrs > 0), "  N_missing:  ", sum( aaGenePtrs == 0))


	codonMap <- getCodonMap()
	codonTable <- vector()

	visitOrder <- ORDER( aaGenePtrs)

	cat( "\n\nVisiting all named proteins")
	curSeqID <- ""
	genomicDNA <- ""

	for (i in visitOrder) {

		if ( i == 0) next
		thisname <- aaNames[i]
		thisGeneID <- aaGeneIDs[i]
		thisAA <- aaSeqs[i]
		thisGenePtr <- aaGenePtrs[i]
		if ( thisGenePtr < 1) next
		geneID <- geneMap$GENE_ID[ thisGenePtr]
		seqID <- geneMap$SEQ_ID[ thisGenePtr]
		if ( seqID != curSeqID) {
			genomicDNA <- getFastaSeqFromFilePath( filePath=genomicFastaFile, seqID=seqID)
			curSeqID <- seqID
		}
		thisDNA <- convertGenomicDNAtoCodingDNA( geneID, genomicDNA)

		# for this AA and DNA, see what codon is for each AA
		allAA <- STRSPLIT( thisAA, split="")[[1]]
		allDNA <- STRSPLIT( thisDNA, split="")[[1]]

		allCODONs <- SAPPLY( seq( 1, length(allDNA), by=3), function(x) {
					PASTE( allDNA[ x:(x+2)], collapse="")
				})
		nUse <- round( fraction * min( length( allAA), length(allCODONs)))

		smallTable <- base::table( PASTE( allAA[1:nUse], allCODONs[1:nUse], sep=":"))
		if ( length( codonTable) > 0) {
			codonTable <- mergeTables( codonTable, smallTable)
		} else {
			codonTable <- smallTable
		}

		if ( i %% 200 == 0) {
			cat( "\n", i, thisname, "\n")
			print( head( base::sort( codonTable, decreasing=T)))
		}
	}
	cat( "\n")

	return( codonTable)
}


`calculateCodonScores` <- function( codonTable) {

	hasN <- GREP( ":.*N", names( codonTable))
	if ( length(hasN) > 0) codonTable <- codonTable[ -hasN]

	# now get the top codons for each AA
	AAout <- SUB( ":.+", "", names( codonTable))

	allAA <- allCodons <- allCounts <- allPercents <- allPct1 <- allPct2 <- allPct3 <- vector()
	nCodons <- 0
	
	TAPPLY( 1:length(AAout), INDEX=factor( AAout), FUN=function(x) {

			# given the set of rows in 'codonTable' that all share one AA
			# turn the relative abundance of each base into a phred score
			m <- matrix( 0, nrow=4, ncol=3)
			rownames(m) <- c( "A","C","G","T")
			sumCnts <- sum( codonTable[ x])
			myCodons <- myX <- myCnts <- vector()
			myBases <- vector( mode="list")
			nNow <- 0
			for ( i in x) {
				thisCnt <- codonTable[i]
				# keep all that are more than 5% of the time
				if ( thisCnt < sumCnts * 0.05) next
				thisCodon <- SUB( ".:", "", names( codonTable)[i])
				theseBases <- STRSPLIT( thisCodon, split="")[[1]]
				if ( ! all( theseBases %in% rownames(m))) next
				m[ theseBases[1], 1] <-  m[ theseBases[1], 1] + thisCnt
				m[ theseBases[2], 2] <-  m[ theseBases[2], 2] + thisCnt
				m[ theseBases[3], 3] <-  m[ theseBases[3], 3] + thisCnt
				nNow <- nNow + 1
				myX[nNow] <- i
				myCodons[nNow] <- thisCodon
				myBases[[nNow]] <- theseBases
				myCnts[nNow] <- thisCnt
			}
			totalCnts <- apply( m, MARGIN=2, FUN=sum)
			visitOrd <- ORDER( myCnts, decreasing=T)
			for ( j in 1:length(myX)) {
				i <- visitOrd[j]
				thisCodon <- myCodons[i]
				thisCnt <- myCnts[i]
				thisBases <- myBases[[i]]
				scores <- c( 0, 0, 0)
				for ( j in 1:3) {
					if ( ! thisBases[j] %in% rownames(m)) next
					scores[j] <- round( m[ thisBases[j], j] * 100 / totalCnts[j])
				}
				nCodons <<- nCodons + 1
				allAA[nCodons] <<- AAout[myX[1]]
				allCodons[nCodons] <<- thisCodon
				allCounts[nCodons] <<- thisCnt
				allPercents[nCodons] <<- thisCnt * 100 / sumCnts
				allPct1[nCodons] <<- scores[1]
				allPct2[nCodons] <<- scores[2]
				allPct3[nCodons] <<- scores[3]
			}
			return()
		})

	# turn the list into a Nx4 matrix
	basePcts <- data.frame( "Pct1"=allPct1, "Pct2"=allPct2, "Pct3"=allPct3, stringsAsFactors=F)
	rownames(basePcts) <- 1:nrow(basePcts)
	phredIntScore <- apply( basePcts, MARGIN=1, function(x) {
					phred <- round( x * 40 / 100)
					return( PASTE( phred, collapse=" "))
				})

	phredAsciiScore <- SAPPLY( phredIntScore, solexaToPhred, scoreType="Phred33")

	out <- data.frame( "AA"=allAA, "Codon"=allCodons, "Count"=allCounts, "Percent"=allPercents,
				basePcts, "PhredIntegers"=phredIntScore, "PhredString"=phredAsciiScore,
				stringsAsFactors=F)
	return( out)
}


`buildExtendedCodonMap` <- function() {

	# there are 64 standard 3-letter codons, but it can be extended...

	# add the resolvable extended IUPAC nucleotides
	bases <- c( "A", "C", "G", "T")
	ambiguityLetters <- c( "W", "S", "M", "K", "R", "Y")
	ambiguityBases <- list( "W"=c("A","T"), "S"=c("C","G"), "M"=c("A","C"), "K"=c("G","T"), 
			"R"=c("A","G"), "Y"=c("C","T"))

	codonMap <- getCodonMap()
	if ( nrow(codonMap) > 64) codonMap <- codonMap[ 1:64, ]
	altMap <- data.frame()

	# loop over all possible triplets
	allbases <- c( bases, ambiguityLetters)
	NB <- length(allbases)

	cat( "\nLooking for unique ambiguous IUPAC triples..\n")
	for( i in 1:NB) {
		dna1 <- allbases[i]
		alts1 <- if ( is.na( wherei <- MATCH( dna1, ambiguityLetters))) dna1 else ambiguityBases[[wherei]]
		len1 <- length(alts1)
		for( j in 1:NB) {
			dna2 <- allbases[j]
			alts2 <- if ( is.na( wherej <- MATCH( dna2, ambiguityLetters))) dna2 else ambiguityBases[[wherej]]
			len2 <- length(alts2)
			for( k in 1:NB) {
				dna3 <- allbases[k]
				alts3 <- if ( is.na( wherek <- MATCH( dna3, ambiguityLetters))) dna3 else ambiguityBases[[wherek]]
				len3 <- length(alts3)

				# nothing to do for the standard triples
				if ( all( c(len1, len2, len3) == 1)) next

				thisAlternateTriple <- PASTE( dna1, dna2, dna3, sep="")

				# the test is "do all alternatives give the same one AA?"
				aaHits <- vector()
				for ( c1 in alts1)
				for ( c2 in alts2)
				for ( c3 in alts3) {
					thisStandardTriple <- PASTE( c1, c2, c3, sep="")
					who <- MATCH( thisStandardTriple, codonMap$DNA, nomatch=0)
					aaHits <- c( aaHits, codonMap$AA[who])
				}
				aaHits <- unique.default( aaHits)
				cat( "\r", thisAlternateTriple, " Hits: ", length(aaHits), aaHits)

				if ( length(aaHits) == 1) {
					bestRow <- MATCH( aaHits[1], codonMap$AA)
					smlDF <- codonMap[ bestRow, ]
					smlDF$DNA <- thisAlternateTriple
					smlDF$RNA <- GSUB( "T","U", thisAlternateTriple, fixed=T)
					altMap <- base::rbind( altMap, smlDF)
					cat( "\nNew: ", nrow( altMap), smlDF$DNA, smlDF$AA, smlDF$AA_Long, "\n")
					next
				}

				# there are a few ambiguous AA that are also allowed
				if ( length( aaHits) == 2) {
					altAA <- ""
					if ( all( aaHits %in% c("N","D"))) {
						altAA <- "B"
						altAAlong <- "Asx"
						rowKey <- "N"
					}
					if ( all( aaHits %in% c("Q","E"))) {
						altAA <- "Z"
						altAAlong <- "Glx"
						rowKey <- "Q"
					}
					if ( all( aaHits %in% c("I","L"))) {
						altAA <- "J"
						altAAlong <- "Xle"
						rowKey <- "I"
					}
					if ( altAA != "") {
						bestRow <- MATCH( rowKey, codonMap$AA)
						smlDF <- codonMap[ bestRow, ]
						smlDF$DNA <- thisAlternateTriple
						smlDF$RNA <- GSUB( "T","U", thisAlternateTriple, fixed=T)
						smlDF$AA <- altAA
						smlDF$AA_Long <- altAAlong
						altMap <- base::rbind( altMap, smlDF)
						cat( "\nNew + Ambig AA: ", nrow( altMap), smlDF$DNA, smlDF$AA, smlDF$AA_Long, "\n")
						next
					}
				}
			}
		}
	}
	cat( "\nDone.  N_Unique: ", nrow(altMap), "\n")

	if ( nrow(altMap) > 0) rownames(altMap) <- 1:nrow(altMap)

	return( altMap)
}


`test.codonTools` <- function() {
	checkEquals( STRSPLIT(DNAtoAA( "ATGTAG"),split="")[[1]], 
			DNAVtoAAV( STRSPLIT("ATGTAG", split="")[[1]]))
	checkEquals( DNAtoAA( "ATGTAG"), c("F1"="M*", "F2"="C", "F3"="V"))
}
