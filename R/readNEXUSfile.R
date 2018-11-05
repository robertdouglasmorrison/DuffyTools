# readNEXUS.R

readNEXUSfile <- function( file="3D7_DBL1domains.nxs") {

	txt <- readLines( file)
	start <- which( "matrix" == txt)[1]
	end <- which( ";" == txt)[1]
	gapLines <- which( txt == "")

	# we can build a giant flat table of the alignments
	gapLines <- gapLines[ gapLines > start]
	curLine <- start + 1
	gptr <- 1
	NG <- median( diff( gapLines)) - 1
	big <- rep( "", times=NG)

	while ( curLine < end) {
		thisLine <- txt[ curLine]
		if ( thisLine == "" || curLine %in% gapLines) {
			curLine <- curLine + 1
			gptr <- 1
			next
		}
		thisLine <- sub( "^.+ +","", thisLine)
		big[gptr] <- paste( big[gptr], thisLine, sep="")
		gptr <- gptr + 1
		curLine <- curLine + 1
	}

	cat( "\nN_Seq: ", NG)
	cat( "\nSeq Lengths: ", nchar( big))

	bigM <- matrix("", nrow=NG, ncol=nchar( big[1]))
	for ( i in 1:NG) bigM[ i, ] <- strsplit( big[i], split="")[[1]]

	colnames( bigM) <- 1:ncol(bigM)
	rnames <- rep("", times=NG)
	for ( i in 1:NG) {
		nameline <- start + i
		thisName <- sub( "(^.+)( +)(.+)","\\1", txt[nameline])
		rnames[i] <- thisName
	}
	cat( "\n", rnames)
	rownames(bigM) <- rnames

	return( bigM)
}

