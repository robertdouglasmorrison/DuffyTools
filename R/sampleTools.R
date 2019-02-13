# sampleTools.R -- utilities that manipulate "SampleID" terms

`checkSampleNames` <- function( sids) {

	sids <- as.character( sids)

	# they get used as R column names...
	R.sids <- make.names( sids)
	bad <- which( R.sids != sids)
	if ( length(bad)) {
		cat( "\n\nSampleID name error:  ", sids[bad], "\n")
		stop( "'SampleID' names must be valid R symbol names.  No blanks, math symbols, leading digits, etc.")
	}

	# they get used in file & folder names
	file.sids <- file.cleanSpecialCharactersFromFileName( sids)
	bad <- which( file.sids != sids)
	if ( length(bad)) {
		cat( "\n\nSampleID name error:  ", sids[bad], "\n")
		stop( "'SampleID' names must be valid file names.  No blanks, math symbols, etc.")
	}

	return( sids)
}
