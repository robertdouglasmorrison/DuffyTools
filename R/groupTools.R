# groupTools.R -- utilities that manipulate "GroupID" terms

`checkGroupNames` <- function( grps) {

	grps <- as.character( grps)

	# they get used as R column names...
	R.grps <- make.names( grps)
	bad <- which( R.grps != grps)
	if ( length(bad)) {
		cat( "\n\nGroupID name error:  ", grps[bad], "\n")
		stop( "'Group' names must be valid R symbol names.  No blanks, math symbols, etc.")
	}

	# they get used in file & folder names
	file.grps <- file.cleanSpecialCharactersFromFileName( grps)
	bad <- which( file.grps != grps)
	if ( length(bad)) {
		cat( "\n\nGroupID name error:  ", grps[bad], "\n")
		stop( "'Group' names must be valid file names.  No blanks, math symbols, etc.")
	}

	return( grps)
}
