# optionsTable.R

# a set of routines to query user-settable options

`readOptionsTable` <- function( file="Options.txt") {

	tbl <- read.delim( file, as.is=TRUE, colClasses="character")
	if ( ncol(tbl) < 2) stop( paste( "Invalid options table:  expected 2 columns:  got", ncol(tbl)))
	colnames(tbl)[1:2] <- c("OptionName", "Value") 

	# quick check for invalid non-tab separators
	blanksInNames <- grep( " +", tbl$OptionName)
	if ( length( blanksInNames)) {
		cat( "\nWarning:  Options Table may have non-tab separators.")
		cat( "\n  Check:  ", tbl$OptionName[ blanksInNames], "\n")
	}
	return( tbl)
}


# get the value of one option as a charater string

`getOptionValue` <- function( tbl, arg, speciesID=NULL, notfound=NA, verbose=TRUE) {

	if ( typeof( tbl) != "list") {
		try( tbl <- readOptionsTable( tbl))
	}

	if ( ! all( colnames(tbl)[1:2] == c("OptionName", "Value"))) stop( "Not a valid options table")

	# allow the arg to be almost anything, character string, object name, etc...
	bigWidth <- max( base::nchar( tbl$OptionName))
	bigWidth <- min( bigWidth, 24)

	name <- try( base::eval( arg), silent=TRUE)
	if ( class(name) == "try-error") {
		name <- as.character( substitute( arg));
	}

	# find that argument, with posibility of a speciesID suffix
	row <- base::match( name, tbl$OptionName, nomatch=0)
	if ( !is.null( speciesID)) {
		name2 <- paste( name, speciesID, sep=".")
		row2 <- base::match( name2, tbl$OptionName, nomatch=0)
		if ( row2 > 0) {
			row <- row2
			name <- name2
			bigWidth <- max( bigWidth, nchar(name2))
		}
	}
	if (verbose) cat( "\nOption Value:     \t", format( name, width=bigWidth), sep="")

	if ( is.na( notfound) && (row == 0)) {
		stop( paste( "Required 'OptionsTable' entry not found:  looked for: ", name))
	}

	if ( row > 0) {
		out <- tbl$Value[ row]
		# allow environment variables in the value field...
		out <- env.sub( out)
		if (verbose) cat( "\tFound:\t", out)
	} else {
		out <- notfound
		if (verbose) cat( "\tDefault:\t", out)
	}

	return( out)
}


# get the value of one option as a boolean True/False

`getOptionTrue` <- function( tbl, arg, speciesID=NULL, notfound=NA, verbose=TRUE) {

	if ( typeof( tbl) != "list") {
		try( tbl <- readOptionsTable( tbl))
	}

	if ( ! all( colnames(tbl)[1:2] == c("OptionName", "Value"))) stop( "Not a valid options table")

	# allow the arg to be almost anything
	bigWidth <- max( base::nchar( tbl$OptionName))
	bigWidth <- min( bigWidth, 24)

	name <- try( base::eval( arg), silent=TRUE)
	if ( class(name) == "try-error") {
		name <- as.character( substitute( arg));
	}

	# find that argument, with posibility of a speciesID suffix
	row <- base::match( name, tbl$OptionName, nomatch=0)
	if ( !is.null( speciesID)) {
		name2 <- paste( name, speciesID, sep=".")
		row2 <- base::match( name2, tbl$OptionName, nomatch=0)
		if ( row2 > 0) {
			row <- row2
			name <- name2
			bigWidth <- max( bigWidth, nchar(name2))
		}
	}

	if (verbose) cat( "\nOption True:      \t", format( name, width=bigWidth), sep="")
	if ( is.na( notfound) && (row == 0)) {
		stop( paste( "Required 'OptionsTable' entry not found:  looked for: ", name))
	}

	if ( ! is.na( notfound)) {
		if ( typeof( notfound) != "logical") stop( "Error: getOptionTrue expects logical 'notfound' arg.")
	}
	
	if ( row > 0) {
		out <- tbl$Value[ row]
		# allow environment variables in the value field...
		out <- env.sub( out)
		out <- as.TRUEorFALSE( out)
		if (verbose) cat( "\tFound: \t", out)
	} else {
		out <- notfound
		if (verbose) cat( "\tDefault:\t", out)
	}
	return( out)
}
