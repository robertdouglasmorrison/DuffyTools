# getFileLineCount.R

# quick File Line counts for reading giant data files

fileLineCountFile <- "./.FileLineCounts.txt"


`getFileLineCount` <-
function( f, sampleID="", verbose=TRUE, mode="full", what=c("lineCount", "readCount")) {

	# only use the first file, no matter what
	f <- f[1]

	# allow the default filename
	if ( f == "") return( -1)

	f <- allowCompressedFileName(f)

	if ( ! file.exists(f)) {
		if (verbose) warning( paste( "getFileLineCount:  file not found:  ", f))
		return( 0)
	}
	if ( ! file.readable(f)) {
		if (verbose) warning( paste( "getFileLineCount:  not a read-able file:  ", f))
		return( 0)
	}
	
	what <- match.arg( what)

	# see if we've seen this file recently
	if ( (n <- quickFileLineCountLookup(f, sampleID, what=what)) > 0) return( n)

	# count the lines in chunks
	if (verbose) cat( "\nCounting lines in file: ", basename(f))
	chunk <- 500000
	nlines <- 0
	options( "warn"= -1)
	con <- openCompressedFile(f,open="r")
	repeat {
		nout <- length( readLines( con, n=chunk, warn=FALSE))
		nlines <- nlines + nout
		if ( mode != "full") break
		if ( nout < chunk) break
	}
	close( con)
	options( "warn"= 0)
	if (verbose && mode == "full") {
		cat( ".Done.  N=", formatC( round(nlines), big.mark=",", format="d"))
	}

	# save the count for later
	if ( mode == "full") quickFileLineCountRecord( f, sampleID=sampleID, 
				lineCount=nlines, readCount=nlines)
	
	nout <- if ( what == "lineCount") nlines else round( nlines/4)

	return( nout)
}


quickFileLineCountRecord <- function( f, sampleID="", lineCount, readCount=lineCount) {

	# only use the first file, no matter what
	# updating to handle paired end files, where both have indentical counts
	files <- f
	for (f in files) {

	# save this line count for later
	fInfo <- file.info( f)
	fInfo$lineCount <- lineCount
	fInfo$readCount <- readCount
	rownames( fInfo) <- basename( rownames(fInfo))
	# strip any compression suffix
	rownames( fInfo) <- sub( ".gz$", "", rownames(fInfo))
	
	useLineCountFile <- fileLineCountFile
	if ( base::nchar( sampleID) > 0) {
		useLineCountFile <- paste( "./.FileLineCounts", sampleID, "txt", sep=".")
	}
	file.lock( useLineCountFile, ID=f, sleeptime=1)
	on.exit( file.unlock( useLineCountFile))

	if ( file.exists( useLineCountFile)) {
		curInfo <- read.delim( useLineCountFile, as.is=TRUE)
		# make sure it has the columns we expect
		if ( ! ( "lineCount" %in% colnames(curInfo))) curInfo$lineCount <- NA
		if ( ! ( "readCount" %in% colnames(curInfo))) curInfo$readCount <- NA
	} else { 
		write.table( fInfo, file=useLineCountFile, sep="\t", quote=FALSE)
		next
	}

	row <- base::match( sub( ".gz$", "", basename(f)), rownames( curInfo), nomatch=0)
	if ( row == 0) {
		write.table( rbind( curInfo, fInfo), file=useLineCountFile, sep="\t", quote=FALSE)
		next
	}

	curInfo[ row, ] <- fInfo[ 1, ]
	write.table( curInfo, file=useLineCountFile, sep="\t", quote=FALSE)
	}  # for each file in 'f'
	return()
}


quickFileLineCountLookup <- function( f, sampleID="", what=c("lineCount", "readCount")) {

	# only use the first file, no matter what
	f <- f[1]

	what <- match.arg( what)
	useLineCountFile <- fileLineCountFile
	if ( base::nchar( sampleID) > 0) {
		useLineCountFile <- paste( "./.FileLineCounts", sampleID, "txt", sep=".")
	}

	if ( ! file.exists( useLineCountFile))  {
		return( 0)
	}

	file.lock( useLineCountFile, ID=f, sleeptime=1)
	on.exit( file.unlock( useLineCountFile))

	curInfo <- read.delim( useLineCountFile, as.is=TRUE)
	row <- base::match( basename(f), rownames( curInfo), nomatch=0)
	checkSize <- ( f == allowCompressedFileName(f))
	if ( row == 0) {
		ftry <- sub( ".gz$", "", f)
		row <- base::match( basename(ftry), rownames( curInfo), nomatch=0)
		checkSize <- FALSE
	}
	if ( row == 0) return( 0)
	#cat( "\nRow= ", row, "\n",f,"\n")

	# see if this file is the same
	if ( checkSize) {
		fInfo <- file.info( allowCompressedFileName(f))
		if ( fInfo$size == curInfo$size[row]) {
			return( curInfo[[ what]][ row])
		} else {
			return(0)
		}
	}

	return( curInfo[[ what]][ row])
}

