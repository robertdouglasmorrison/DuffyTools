# snpTools.R -- routines for dealing with SNPs


`loadKnownSNPtable` <- function( seqID, verbose=TRUE) {

	# lets try to be smarter / faster
	if ( exists( "SNP_curSNPtable") && nrow(SNP_curSNPtable) > 0) {
		if ( SNP_curSNPtable$SEQ_ID[1] == seqID) return()
	}

	# look for the table for this species/seqID
	dataFiles <- data( package="DuffyTools")$results[ , "Item"]
	dataFiles <- grep( "SNP", dataFiles, value=T)

	# first look for a SEQ_ID hit
	SNPtable <- NULL
	snpSeqFile <- paste( seqID, "SNP_Table", sep=".")
	if ( length( grep( snpSeqFile, dataFiles, fixed=T)) > 0) {
		data( list=snpSeqFile, envir=environment())
		if ( is.null( SNPtable)) {
			if (verbose) cat( "\nFailed to load DuffyTools SNP object: ", snpSeqFile)
			SNP_curSNPtable <<- data.frame()
			return()
		}
		# got it
		SNP_curSNPtable <<- checkSNPtable( SNPtable)
		return()
	}

	# next look for a Prefix hit
	prefix <- getCurrentSpeciesFilePrefix()
	snpPrefixFile <- paste( prefix, "SNP_Table", sep=".")
	if ( length( grep( snpPrefixFile, dataFiles, fixed=T)) > 0) {
		data( list=snpPrefixFile, envir=environment())
		if ( is.null( SNPtable)) {
			if (verbose) cat( "\nFailed to load DuffyTools SNP object: ", snpPrefixFile)
			SNP_curSNPtable <<- data.frame()
			return()
		}
		# got it
		SNPtable <- subset( SNPtable, SEQ_ID == seqID)
		if (nrow(SNPtable) < 1) {
			SNP_curSNPtable <<- data.frame()
			return()
		}
		rownames( SNPtable) <- 1:nrow( SNPtable)
		SNP_curSNPtable <<- checkSNPtable( SNPtable)
		return()
	}

	# last, look for an explicit folder of SNP tables
	if ( exists( "SNP_curSNPpath")) {
		snpSeqFile <- file.path( SNP_curSNPpath, paste(snpSeqFile,"rda", sep="."))
		if ( file.exists( snpSeqFile)) {
			if (verbose) cat( "  loading SNP file: ", basename(snpSeqFile))
			load( snpSeqFile)
			SNP_curSNPtable <<- checkSNPtable( SNPtable)
			return()
		}
		snpPrefixFile <- file.path( SNP_curSNPpath, paste(snpPrefixFile,"rda", sep="."))
		if ( file.exists( snpPrefixFile)) {
			if (verbose) cat( "  loading SNP file: ", basename(snpSeqFile))
			load( snpPrefixFile)
			SNPtable <- subset( SNPtable, SEQ_ID == seqID)
			rownames( SNPtable) <- 1:nrow( SNPtable)
			SNP_curSNPtable <<- checkSNPtable( SNPtable)
			return()
		}
	}
	if (verbose) cat( "\nFailed to find DuffyTools SNP object: ", basename(snpPrefixFile))
	SNP_curSNPtable <<- data.frame()
	return()
}


`checkSNPtable` <- function( SNPtable) {

	# make sure we have the columns we expect
	if ( all( c("MAJOR_ALLELE", "MINOR_ALLELE") %in% colnames( SNPtable))) return( SNPtable)

	out <- SNPtable
	# the MAJOR is the reference...
	if ( !("MAJOR_ALLELE" %in% colnames( SNPtable))) {
		out$MAJOR_ALLELE <- SNPtable$REF_ALLELE
	}
	if ( !("MINOR_ALLELE" %in% colnames( SNPtable))) {
		out$MINOR_ALLELE <- SNPtable$ALT_ALLELE
	}
	out
}

