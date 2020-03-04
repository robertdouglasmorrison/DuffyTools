# imgtTools.R -- various tools for manipulating IMGT B cell and T cell data (IG & TR)



`getIMGTdata` <- function() {

	file <- paste( "Hs", "IMGT.VDJC.rda", sep=".")
	imgt <- NULL
	data( list=file, envir=environment())
	if ( is.null( imgt)) {
		# try local folder
		if ( file.exists( file)) load( file)
		if ( is.null( imgt)) {
			cat( "\nFailed to load IMGT data object: ", file)
			return( NULL)
		}
	}
	if ( ! all( colnames(imgt) == c("IMGT_ID", "DNA", "AA"))) {
		cat( "\nUnexpected column names in IMGT data object..")
	}
	return( imgt)
}


`getBestIMGTgeneData` <- function( genes) {

	imgt <- getIMGTdata()

	# best would be a perfect match
	N <- length( genes)
	where <- match( genes, imgt$IMGT_ID, nomatch=NA)
	out <- imgt[ where, ]
	out$GivenGene <- genes
	out$N_Match <- 0
	out$MatchNotes <- "No Match Found"
	rownames(out) <- 1:N

	# did some hit?
	out$N_Match[ where > 0] <- 1
	out$MatchNotes[ where > 0] <- "Unique Perfect Match"
	if ( ! any( is.na( where))) return( out)

	# general case to try and find the best match
	for ( k in which( is.na(where))) {
		out[ k, ] <- findBestIMGTgeneData( genes[k], imgt)
	}
	out <- out[ , c(4,5,1,6,2:3)]
	return( out)
}


`findBestIMGTgeneData` <- function( gene, imgt) {

	# some sources have leading zeros after the VDJ letter
	gene <- sub( "(^IG[HKL][VDJCAGM])(0)([1-9].+)", "\\1\\3", gene)
	gene <- sub( "(^TR[ABDG][VDJC])(0)([1-9].+)", "\\1\\3", gene)

	# try again for a perfect match
	where <- match( gene, imgt$IMGT_ID)
	if ( ! is.na( where)) {
		out <- imgt[ where, ]
		out$GivenGene <- gene
		out$N_Match <- 1
		out$MatchNotes <- "Match: removed leading zero"
		out <- out[ , c(4,5,1,6,2:3)]
		return(out)
	}
	# perhaps given is not specific enough
	patt <- paste( "^", gene, sep="")
	who <- grep( patt, imgt$IMGT_ID)
	if ( length(who)) {
		hit <- who[1]
		out <- imgt[ hit, ]
		out$GivenGene <- gene
		nHit <- length(who)
		out$N_Match <- nHit
		if ( nHit == 1) {
			out$MatchNotes <- "Unique Partial Match"
		} else {
			out$MatchNotes <- paste( "Multiple Partial Matches: ", 
						paste(imgt$IMGT_ID[who], collapse=";"))
		}
		out <- out[ , c(4,5,1,6,2:3)]
		return(out)
	}

	# some given names have an explicit '-01*' or '-1*' when IMGT just not
	# note the leading zero already taken out
	gTry <- sub( "-01*", "*", gene, fixed=T)
	if ( gTry != gene) {
		where <- match( gTry, imgt$IMGT_ID)
		if ( ! is.na( where)) {
			out <- imgt[ where, ]
			out$GivenGene <- gene
			out$N_Match <- 1
			out$MatchNotes <- "Match: removed expicit level '-01*'"
			out <- out[ , c(4,5,1,6,2:3)]
			return(out)
		}
		# perhaps given is not specific enough
		patt <- paste( "^", gTry, sep="")
		who <- grep( patt, imgt$IMGT_ID)
		if ( length(who)) {
			hit <- who[1]
			out <- imgt[ hit, ]
			out$GivenGene <- gene
			nHit <- length(who)
			out$N_Match <- nHit
			if ( nHit == 1) {
				out$MatchNotes <- "Unique Partial Match after remove level '-01*'"
			} else {
				out$MatchNotes <- paste( "Multiple Partial Matches: ", 
							paste(imgt$IMGT_ID[who], collapse=";"))
			}
			out <- out[ , c(4,5,1,6,2:3)]
			return(out)
		}
	}
	# getting desparate...
	# first make sure the root is valid
	geneRoot <- sub( "[\\-\\*].+", "", gene)
	nHits <- grep( paste( "^", geneRoot, sep=""), imgt$IMGT_ID)
	if ( ! nHits) {
		out <- imgt[ NA, ]
		out$GivenGene <- gene
		out$N_Match <- 0
		out$MatchReason <- "Unknown Gene Root. No Match"
		out <- out[ , c(4,5,1,6,2:3)]
		return(out)
	}
	# last gasp, repeatedly shorten until something hits
	gTry <- gene
	while (nchar( gTry)) {
		gTry <- substr( gTry, 1, nchar(gTry)-1)
		patt <- paste( "^", gTry, sep="")
		who <- grep( patt, imgt$IMGT_ID)
		if ( length(who)) {
			hit <- who[1]
			out <- imgt[ hit, ]
			out$GivenGene <- gene
			nHit <- length(who)
			out$N_Match <- nHit
			if ( nHit == 1) {
				out$MatchNotes <- "Unique Partial Match after removing suffix chars"
			} else {
				out$MatchNotes <- paste( "Multiple Partial Matches after removing suffix chars: ", 
							paste(imgt$IMGT_ID[who], collapse=";"))
			}
			out <- out[ , c(4,5,1,6,2:3)]
			return(out)
		}
	}
	out <- imgt[ NA, ]
	out$GivenGene <- gene
	out$N_Match <- 0
	out$MatchReason <- "Unknown Gene Root. No Match"
	out <- out[ , c(4,5,1,6,2:3)]
	return(out)
}

