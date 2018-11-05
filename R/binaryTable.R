# binaryTable.R

# tools to speed up reading giant RNAseq tables

`makeBinaryTable` <- function( x, filename=NULL) {

	if ( is.null( filename)) stop( "makeBinaryTable:  required 'filename' argument is missing")

	fileToUse <- paste( filename, ".rda", sep="")
	save( x, file=fileToUse)
	cat( "\nWrote Binary Table object: \t", fileToUse,"\n")
	return()
}


`loadBinaryTable` <- function( filename) {

	status <- FALSE
	object <- NULL

	fileToTry <- paste( filename, ".rda", sep="")
	if ( file.exists( fileToTry)) {

		# yes, there is a binary, is it newer than its original?
		finfoFull <- file.info( filename)
		finfoBinary <- file.info( fileToTry)

		# if we have only the binary table, it is always good
		if ( all( is.na( finfoFull))) {
			isNewer <- TRUE
		} else {
			isNewer <- (finfoBinary$ctime > finfoFull$ctime)
			isNewer <- ( !is.na( isNewer) && isNewer)
		}

		# get that object
		if (isNewer) {
			who <- load( fileToTry, envir=as.environment( -1))
			object <- get( who[1])
			status <- TRUE
		}
	}

	return( list( "Status"=status, "Object"=object))
}
