# annotationTable.R

# a set of routines to query user-settable annotation details

`readAnnotationTable` <- function( file="Annotation.txt", neededColumns=NULL, sep="\t", verbose=FALSE) {

	tbl <- read.delim( file, as.is=TRUE, sep=sep, colClasses="character")
	if ( ncol(tbl) < 2) stop( paste( "Invalid annotation table:  expected at least 2 columns:  got", ncol(tbl)))

	if ( ! is.null( neededColumns)) {
		hasThem <- all( neededColumns %in% colnames( tbl))
		if ( ! hasThem) stop( paste( "Missing required annotation columns.  Expected: ", 
				paste( neededColumns, collapse="  ")))
	}

	# if there is an 'exclude' or 'include' column, trim by that
	trimColumns <- match( c( "include", "exclude"), tolower( colnames(tbl)), nomatch=0)
	if ( any( trimColumns > 0)) {
		incl <- trimColumns[1]
		if ( incl > 0) {
			keep <- which( tbl[[ incl]] == TRUE)
		}
		excl <- trimColumns[2]
		if ( excl > 0) {
			keep <- which( tbl[[ excl]] == FALSE)
		}
		tbl <- tbl[ keep, ]
	}
	
	return( tbl)
}


# get the value of one annotation entry as a charater string

`getAnnotationValue` <- function( tbl, key, columnArg, notfound="", verbose=TRUE) {

	if ( typeof( tbl) != "list") {
		try( tbl <- readAnnotationTable( tbl))
	}

	row <- base::match( key, tbl[ ,1], nomatch=0)
	if ( row == 0) stop( paste( "Not a valid Annotation Table key:   tried: ", key))

	# allow the arg to be almost anything, character string, object name, etc...
	name <- try( base::eval( columnArg), silent=TRUE)
	if ( class(name) == "try-error") {
		name <- as.character( substitute( columnArg));
	}
	if (verbose) cat( "\nAnnotation Value: \t", key, ":\t", name, sep="")

	col <- base::match( name, colnames(tbl), nomatch=0)
	if ( col == 0) {
		# when not found, see if it is 'almost' a match, e.g. a misspelling
		editDist <- adist( name, colnames(tbl))
		closeEnough <- max( 3, round(nchar(name)*0.25))
		if ( any( editDist <= closeEnough)) {
			cat( "\n\nWarning: wanted Annotation argument: ", name, 
				" is a close but inexact match to some Annotation columns.\n  Using 'default' value instead...\n")
		}
		out <- as.character(notfound)
		if (verbose) cat( ":\tDefault:\t", out)
	} else {
		out <- tbl[ row, col]
		out <- env.sub( out)
		if (verbose) cat( ":\tFound:\t", out)
	}

	return( out)
}


# get the value of one option as a boolean True/False

`getAnnotationTrue` <- function( tbl, key, columnArg, notfound=FALSE, verbose=TRUE) {

	if ( typeof( tbl) != "list") {
		try( tbl <- readAnnotationTable( tbl))
	}

	row <- base::match( key, tbl[ ,1], nomatch=0)
	if ( row == 0) stop( paste( "Not a valid Annotation Table key:   tried: ", key))

	# allow the arg to be almost anything, character string, object name, etc...
	name <- try( base::eval( columnArg), silent=TRUE)
	if ( class(name) == "try-error") {
		name <- as.character( substitute( columnArg));
	}
	if (verbose) cat( "\nAnnotation True: \t", key, ":\t", name, sep="")
	col <- base::match( name, colnames(tbl), nomatch=0)
	if ( col == 0) {
		# when not found, see if it is 'almost' a match, e.g. a misspelling
		editDist <- adist( name, colnames(tbl))
		closeEnough <- max( 3, round(nchar(name)*0.25))
		if ( any( editDist <= closeEnough)) {
			cat( "\n\nWarning: wanted Annotation argument: ", name, 
				" is a close but inexact match to some Annotation columns.\n  Using 'default' value instead...\n")
		}
		out <- as.logical(notfound)
		if (verbose) cat( ":\tDefault:\t", out)
	} else {
		out <- as.TRUEorFALSE( tbl[ row, col])
		if (verbose) cat( "\tFound:\t", out)
	}

	return( out)
}


# get the vector of row numbers that are True

`whichAnnotationTrue` <- function( tbl,  columnArg, negate=FALSE, verbose=TRUE) {

	if ( typeof( tbl) != "list") {
		try( tbl <- readAnnotationTable( tbl))
	}

	# allow the arg to be almost anything, character string, object name, etc...
	name <- try( base::eval( columnArg), silent=TRUE)
	if ( class(name) == "try-error") {
		name <- as.character( substitute( columnArg));
	}
	if (verbose) cat( "\nwhichAnnotationTrue:", "\tcolumn:\t", name)

	col <- base::match( name, colnames(tbl), nomatch=0)
	if ( col == 0) {
		warning( paste( "Annotation Table column not found:  looked for: ", name))
		out <- vector()
	} else {

		str <- tbl[ , col]
		out <- sapply( str, as.TRUEorFALSE)
	}

	# final step, they may want 'who is False'
	if ( negate) {
		out <- which( ! out)
		if (verbose) cat( "     N_False: ", length(out))
	} else {
		out <- which( out)
		if (verbose) cat( "     N_True: ", length(out))
	}

	return( out)
}
