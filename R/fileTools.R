# fileTools.R

# various file utilites....


# see if a file is readable...
`file.readable` <-
function(f) return( ( file.exists(f) && (file.access( f, mode=4) == 0)))


# delete a file if it already exists...
`file.delete` <-
function(files) {
	for( f in base::unlist(files)) {
		if ( (!is.null(f)) && file.exists(f) && ( file.info(f)$isdir == FALSE)) {
			file.remove(f)
		}
	}
}


`file.find` <- function( f, searchPath=c( ".", Sys.getenv('HOME')) ) {

	for ( path in searchPath) {
		fileTry <- file.path( path, f)
		if ( file.readable( fileTry)) return( fileTry)
	}
	cat( "\nFile not found or not read enabled: ", f)
	cat( "\nSearched in:  ", searchPath, "\n")
	return( NA)
}


# file combine, with Header Row checking
`file.combine` <- function( infiles, outfile=NULL, check.headers=FALSE, renameOK=TRUE,
			showWarnings=TRUE) {

	if ( is.null( outfile)) stop( "file.combine:  required 'outfile' argument is missing")
	if ( length( gsub( " ", "", outfile)) < 1) stop( "file.combine:  required 'outfile' argument is empty")

	fset <- infiles
	if ( typeof( fset) == "list") fset <- base::unlist( fset)
	if ( length( fset) < 1) {
		if (showWarnings) warning( "file.combine:  no input files specified")
		file.delete( outfile)
		return()
	}

	# first quietly drop non-existent files
	newFset <- vector()
	for ( f in fset) {
		if( file.exists( f)) newFset <- base::append( newFset, f)
	}
	if ( length( newFset) < 1) {
		if (showWarnings) warning( "file.combine:  none of input files exist")
		file.delete( outfile)
		return()
	}
	fset <- newFset

	# clear the output
	if ( ! any( fset == outfile)) file.delete( outfile)

	# names the same, already done...
	if ( length(fset) == 1 && fset == outfile) return()

	# if not doing header check and removal, use OS
	if ( ! check.headers) {

		if ( length( fset) == 1) {
		
			# names the same == do nothing
			if ( fset[1] == outfile) return
			if (renameOK) {
				system( paste( "mv  ", fset[1], " ", outfile))
			} else {
				system( paste( "cp  ", fset[1], " ", outfile))
			}
			return()
		}

		system( paste( "cat ", paste( fset, collapse=" "), "  > ", outfile))
		return()
	}

	# we are doing header test and removal
	header <- ""
	outCon <- file( outfile, open="w")

	for ( i in 1:length( fset)) {

		conIn <- file( fset[i], open="r")
		chunkSize <- 100000
		nChunks <- 0

		repeat { 
			txt <- readLines( con=conIn, n=chunkSize, warn=FALSE)
			if ( ! length(txt)) break
			nChunks <- nChunks + 1
			lines <- 1:length( txt)

			if ( i == 1 && nChunks == 1) {
				header <- txt[1]
			}
			if ( i > 1 && nChunks == 1) {
				if ( txt[1] != header) {
					if (showWarnings) warning( paste( "file.combine:  files do not have matching headers\n",
						fset[1], header,"\n",fset[i], txt[1]))
				}
				lines <- if ( length(txt) > 1) (2:length(txt)) else NULL
			}
			if ( ! is.null( lines)) writeLines( txt[ lines], con=outCon)

			if ( length( txt) < chunkSize) break
		}
		close( conIn)
	}
	close( con=outCon)
	return()
}


`file.cleanSpecialCharactersFromFileName` <- function( fnames, max.length=150) {

	# take out any characters that may break any OS file system....
	fnew <- fnames
	fnew <- gsub( ":", "-", fnew, fixed=TRUE)
	fnew <- gsub( "/", "-", fnew, fixed=TRUE)
	fnew <- gsub( "\n", "-", fnew, fixed=FALSE)
	fnew <- gsub( "\\", "-", fnew, fixed=TRUE)
	fnew <- gsub( "|", "-", fnew, fixed=TRUE)
	fnew <- gsub( "@", "-at-", fnew, fixed=TRUE)
	fnew <- gsub( ",", "-", fnew, fixed=TRUE)
	fnew <- gsub( "(", "-", fnew, fixed=TRUE)
	fnew <- gsub( ")", "-", fnew, fixed=TRUE)
	fnew <- gsub( ">", "-", fnew, fixed=TRUE)
	fnew <- gsub( "<", "-", fnew, fixed=TRUE)
	fnew <- gsub( ";", "-", fnew, fixed=TRUE)
	fnew <- gsub( "*", "-", fnew, fixed=TRUE)
	# lastly, trim out adjacent dashes
	fnew <- gsub( " ", "-", fnew, fixed=FALSE)
	fnew <- gsub( "\\-+", "-", fnew, fixed=FALSE)
	# and remove any "-." or ".-" connectors
	fnew <- gsub( "-.", ".", fnew, fixed=TRUE)
	fnew <- gsub( ".-", ".", fnew, fixed=TRUE)
	# don't end with a dash
	fnew <- sub( "-$", "", fnew, fixed=FALSE)

	# lastly prevent super long names
	fnew <- clipLongString( fnew, max.length=max.length)

	return( fnew)
}
			

# allow a compressed .GZ or .BZ2 file to be read implicitly, given a filename
# does a search for openable files, so only works for READING files...
`allowCompressedFileName` <- function( filename) {

	if ( file.readable( filename)) return( filename)

	if ( regexpr( "\\.gz$", filename) < 0) {
		tryFile <- paste( filename, "gz", sep=".")
		if ( file.readable( tryFile)) return( tryFile)
	}

	if ( regexpr( "\\.bz2$", filename) < 0) {
		tryFile <- paste( filename, "bz2", sep=".")
		if ( file.readable( tryFile)) return( tryFile)
	}

	return( filename)
}


`openCompressedFile` <- function( filename, open="r") {

	# watch for the standard compression suffixes
	if ( regexpr( "\\.gz$", filename) > 0) {
		return( gzfile( filename, open=open))
	}

	if ( regexpr( "\\.bz2$", filename) > 0) {
		return( bzfile( filename, open=open))
	}

	return( file( filename, open=open))
}


`tempFolder` <- function( foldername, tmpRoot="/tmp") {

	# first make sure the temp root exists
	# we may be given a comma separated list of folders
	tmpRootList <- strsplit( tmpRoot, split=", *")[[1]]
	tmpUse <- NULL
	for (tmpTry in tmpRootList) {
		if ( file.exists( tmpTry)) {
			tmpUse <- tmpTry
			break
		}
	}
	if ( is.null(tmpUse)) {
		cat( "\nTemporary storage location not found: ", tmpRootList, "\nDefaulting to '/tmp/'")
		tmpUse <- "/tmp"
	}
	folder <- file.path( tmpUse, foldername)

	# create the folder if we need to
	if ( ! file.exists( folder)) dir.create( folder, recursive=TRUE, showWarnings=FALSE)

	# now explicit test of existence
	if ( ! file.exists( folder)) stop( paste( "Failed to create temporary folder:  ", folder))

	return( folder)
}


`removeTempFolder` <- function( path) {

	cmdline <- paste( "rm -fr ", path)
	system( cmdline)
	return()
}
