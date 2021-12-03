# pexelTools.R

# do a search of protein sequences for PEXELs

PEXEL.search <- function( x, max.aa.start=150, strict=TRUE, verbose=FALSE) {

	# we can be given any of several formats
	fa <- seqs <- desc <- NULL
	if ( is.character(x) && length(x) == 1) {
		# a FASTA filename?
		fa <- loadFasta(x, short=F)
		seqs <- fa$seq
		desc <- fa$desc
	}
	if ( is.list( x) && all( c("desc", "seq") %in% names(x))) {
		fa <- x
		seqs <- fa$seq
		desc <- fa$desc
	}
	if ( is.character(x) && length(x) > 1) {
		seqs <- x
		desc <- names(x)
	}

	# standarize
	N <- length( seqs)
	seqs <- toupper( seqs)

	# set up to find the various 'strick-ness' levels of PEXEL
	if (verbose) cat( "\nPexel search of the first ", max.aa.start, " amino acids.\n")

	nPexel <- 0
	startAA <- motif <- gName <- typeOut <- vector( length=N)
	pexelType <- rep( 0, times=N)

	pexelTerm1 <- "R[A-Z][IL][A-Z][EQD]"
	pexelTerm2 <- "[RK][NQSTYAVLIPFMWGC][LIVMAF][NQSTYAVLIPFMWGC][EQD]"
	pexelLikeTerm1 <- "R[A-Z][LMF][ASTCHYV][EQDSAVR]"
	pexelLikeTerm2 <- "[RK][A-Z][LIAVMF][ASTCHYV][EQDSAVR]"
	typeNames <- c( "PEXEL-Strict", "PEXEL-Lax", "PexelLike-VTS-Strict", "PexelLike-VTS-Lax")

	for ( i in 1:N) {
		bestType <- bestHit <- 0
		hit1 <- regexpr( pexelTerm1, seqs[i])
		if ( hit1 > 0 && hit1 <= max.aa.start) {
			bestHit <- hit1
			bestType <- 1
			bestWho <- typeNames[1]
		}
		if ( ! strict && bestHit == 0) {
			hit2 <- regexpr( pexelTerm2, seqs[i])
			if ( hit2 > 0 && hit2 <= max.aa.start) {
				bestHit <- hit2
				bestType <- 2
				bestWho <- typeNames[2]
			}
		}
		if ( ! strict && bestHit == 0) {
			hit3 <- regexpr( pexelLikeTerm1, seqs[i])
			if ( hit3 > 0 && hit3 <= max.aa.start) {
				bestHit <- hit3
				bestType <- 3
				bestWho <- typeNames[3]
			}
		}
		if ( ! strict && bestHit == 0) {
			hit4 <- regexpr( pexelLikeTerm2, seqs[i])
			if ( hit4 > 0 && hit4 <= max.aa.start) {
				bestHit <- hit4
				bestType <- 4
				bestWho <- typeNames[4]
			}
		}
		if ( bestHit > 0) {
			pexelType[ i] <- bestType
			motifAA <- substr( seqs[i], bestHit, (bestHit+4))
			nPexel <- nPexel + 1
			gName[ nPexel] <- desc[i]
			typeOut[ nPexel] <- bestWho
			startAA[ nPexel] <- bestHit
			motif[ nPexel] <- motifAA
			if (verbose) cat( "\r", i, nPexel, desc[i], motifAA)
		}
	}
	if (verbose) {
		cat( "\n\nTotal 'Traditional PEXEL' (Marti, et al, Fig1)(strictest): ", sum( pexelType == 1), "\n")
		if ( ! strict) {
		cat( "\nTotal 'Traditional PEXEL' (Marti, et al, Text):              ", sum( pexelType == 2), "\n")
		cat( "\nTotal 'Pexel-Like (VTS)' (Hiller, et al, strict):            ", sum( pexelType == 3), "\n")
		cat( "\nTotal 'Pexel-Like (VTS)' (Hiller, et al, lax):               ", sum( pexelType == 4), "\n")
		}
	}

	length( gName) <- length( typeOut) <- length( startAA) <- length( motif) <- nPexel
	out <- data.frame( gName, typeOut, motif, startAA, stringsAsFactors=FALSE)
	colnames( out) <- c( "GENE_ID", "TYPE", "MOTIF", "POSITION")
	out$PRODUCT <- gene2Product( out$GENE_ID)

	if (strict) {
		out <- subset( out, TYPE == typeNames[1])
	}
	rownames( out) <- 1:nrow( out)
	return( out)
}

