# lymphocyteTools.R -- various tools for Bcell and Tcell genes, etc.


`getLymphocyteRegions` <- function( type=c( "IG", "TCR")) {

	type <- match.arg( type)

	if ( getCurrentSpecies() == "Hs_grc") {
		# regions for IG genes  (IGK, IGL, IGH)
		IGseqids <- c( "Hs_grc_02", "Hs_grc_22", "Hs_grc_14")
		IGstarts <- c( 88850000, 22020000, 105560000)
		IGstops <-  c( 90250000, 22930000, 106890000)

		# regions for TR genes  (TRA, TRB, TRg)
		# TRD@ is fully inside TRA@, so no explicitly given
		TCRseqids <- c( "Hs_grc_14", "Hs_grc_07", "Hs_grc_07")
		TCRstarts <- c( 21620000, 142290000, 38230000)
		TCRstops <-  c( 22560000, 142820000, 38370000)
	}

	if ( ! exists( "IGseqids")) {
		cat( "\nCurrent species has no defined IG/TCR regions.")
		cat( "\nCurrently set to:  ", getCurrentSpecies())
		return( NULL)
	}

	if ( type == "IG") return( list( "SEQ_ID"=IGseqids, "START"=IGstarts, "STOP"=IGstops))
	if ( type == "TCR") return( list( "SEQ_ID"=TCRseqids, "START"=TCRstarts, "STOP"=TCRstops))
	return( NULL)
}


`getLymphocyteGeneMap` <- function() {

	# grab all the defined genes in the B cell and T cell loci

	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
	out <- data.frame()

	IG <- getLymphocyteRegions( "IG")
	N <- length( IG$SEQ_ID)
	for ( i in 1:N) {
		sml <- subset( gmap, SEQ_ID == IG$SEQ_ID[i] & POSITION >= IG$START[i] & END <= IG$STOP[i])
		if ( nrow(sml)) out <- rbind( out, sml)
	}

	TCR <- getLymphocyteRegions( "TCR")
	N <- length( TCR$SEQ_ID)
	for ( i in 1:N) {
		sml <- subset( gmap, SEQ_ID == TCR$SEQ_ID[i] & POSITION >= TCR$START[i] & END <= TCR$STOP[i])
		if ( nrow(sml)) out <- rbind( out, sml)
	}

	# we may want to hand clean this a bit...
	# 1:  no locus genes that are the whole thing
	isLocus <- grep( "@", out$GENE_ID, fixed=T)
	# 2:  regular generic human gene names "LOCxxxxx"
	isGeneric <- grep( "^LOC[0-9]+", out$GENE_ID)

	drops <- sort( unique( c( isLocus, isGeneric)))
	if ( length(drops)) out <- out[ -drops, ]

	# specific inclusion steps
	if ( getCurrentSpecies() == "Hs_grc") {
		keepIG <- grep( "^IG[HKL]", out$GENE_ID)
		keepTCR <- grep( "^TR[ABDG]", out$GENE_ID)
		keep <- sort( unique( c( keepIG, keepTCR)))
		out <- out[ keep, ]
	}

	rownames(out) <- 1:nrow(out)
	return( out)
}

