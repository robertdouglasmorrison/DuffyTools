# fingerprintTools.R -- utilities for DNA fingerprinting of parasite isolates using CLAG2 and VAR2CSA DBL4
#			Intended for use with PCR/Sanger sequencing




`getFingerprintSubstitutionMatrix` <- function( type=c("DNA","AA")) {

	type <- match.arg( type)

	if ( type == "DNA") {
		if ( ! exists( "DNA_MATRIX")) {
			DNA_MATRIX <<- makeFingerprintSubstitutionMatrix( "DNA")
		}
		return( DNA_MATRIX)
	} else {
		if ( ! exists( "AA_MATRIX")) {
			AA_MATRIX <<- makeFingerprintSubstitutionMatrix( "AA")
		}
		return( AA_MATRIX)
	}
}


`makeFingerprintSubstitutionMatrix` <- function( type=c("DNA","AA"), blanks.OK=TRUE, gaps.OK=TRUE) {

	type <- match.arg( type)

	require( Biostrings)
	require( pwalign)


	if ( type == "DNA") {
		m <- nucleotideSubstitutionMatrix( match=0, mismatch= -1, baseOnly=F, type="DNA")
		N <- ncol(m)
		# invert it into a Levenstein type distances
		#for( i in 1:N) {
			#m[ , i] <- 1
			#m[ i, i] <- 0
		#}
		# let's not penalize for missing data
		isN <- which( colnames(m) == "N")
		m[ isN, ] <- 0
		m[ , isN] <- 0
	}
	if ( type == "AA") {
		if ( ! exists( "BLOSUM62")) data( BLOSUM62)
		m <- BLOSUM62
		N <- ncol(m)
		# invert it into a Levenstein type distances
		for( i in 1:N) {
			m[ , i] <- -1
			m[ i, i] <- 0
		}
		# let's not penalize for missing data
		isX <- which( colnames(m) == "X")
		m[ isX, ] <- 0
		m[ , isX] <- 0
	}

	# allow blanks on the flanks to not cause any penalty
	if ( blanks.OK) {
		isBLANK <- which( colnames(m) == " ")
		N <- ncol(m)
		if ( ! length(isBLANK)) {
			# add a column/row for blanks
			mnew <- matrix( 0, nrow=N+1, ncol=N+1)
			for ( j in 1:N) mnew[ 1:N, j] <- m[ ,j]
			colnames(mnew) <- rownames(mnew) <- c( colnames(m), " ")
			m <- mnew
			isBLANK <- N + 1
		}
		m[ isBLANK, ] <- 0
		m[ , isBLANK] <- 0
	}

	# allow indel gaps in the sequences
	if ( gaps.OK) {
		isGAP <- which( colnames(m) == "-")
		N <- ncol(m)
		if ( ! length(isGAP)) {
			# add a column/row for indel gaps
			mnew <- matrix( 0, nrow=N+1, ncol=N+1)
			for ( j in 1:N) mnew[ 1:N, j] <- m[ ,j]
			colnames(mnew) <- rownames(mnew) <- c( colnames(m), "-")
			m <- mnew
			isGAP <- N + 1
		}
		m[ isGAP, ] <- -1
		m[ , isGAP] <- -1
		m[ isGAP, isGAP] <- 0
	}

	return( m)
}


`fingerprintDistance` <- function( strs, type=c("DNA","AA")) {

	# given a set of strings, return the distance matrix
	type <- match.arg( type)

	subM <- getFingerprintSubstitutionMatrix( type)

	ans <- stringDist( toupper(strs), method="substitutionMatrix", substitutionMatrix=subM,
				type="global")
	# the score tool finds maximum, and we used negative scores, so invert
	# return( -ans)
	# as of Oct 2019, string distance seems to now return positive distance....
	return( ans)
}


bestMatchingFingerprint <- function( tbl, seq, finger=c("CLAG2","DBL4","COMBO"), type=c("DNA","AA")) {

	finger <- match.arg( finger)
	type <- match.arg( type)

	wantColumn <- paste( finger, type, sep="_")
	seqSet <- tbl[[ wantColumn]]
	names(seqSet) <- tbl$ISOLATE_ID
	N <- length( seqSet)
	seqDist <- rep.int( NA, N)

	for ( i in 1:N) seqDist[i] <- fingerprintDistance( c( seq, seqSet[i]), type=type)

	best <- which.min( seqDist)
	bestSeq <- seqSet[best]
	bestDist <- seqDist[best]

	paAns <- pairwiseAlignment( seq, bestSeq, type="global", scoreOnly=F)

	# there may be 2+ all equally best
	best <- which( seqDist == bestDist)

	out <- data.frame( "N_Best"=length(best), "BestMatch"=paste(tbl$ISOLATE_ID[best], collapse="; "), 
				"EditDist"=bestDist, "Score"=score(paAns), "BestSeq"=as.character(alignedSubject(paAns)),
				"GivenSeq"=as.character( alignedPattern(paAns)), stringsAsFactors=F)

	out
}


