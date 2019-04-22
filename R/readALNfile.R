# readALNfile.R -- read up the result consensus alignment filr from Clustalw2 


readALN <- function( file, verbose=TRUE) {

	txt <- readLines( file)
	start <- 2
	end <- length(txt)
	gapLines <- which( txt == "")
	gapLines <- gapLines[ gapLines > start]
	starLines <- grep( "^    ", txt)

	# we can build a giant flat table of the alignments, finding group IDs as we go
	curLine <- start
	gptr <- 1
	NG <- 1
	big <- rep( "", times=NG)
	gnames <- rep( "", times=NG)
	starText <- ""
	starHeadPtr <- vector()

	while ( curLine <= end) {
		thisLine <- txt[ curLine]
		if ( curLine %in% starLines) {
			if ( length( starHeadPtr) > 0) {
				starStart <- median( starHeadPtr)
				thisStars <- substr( thisLine, starStart, nchar(thisLine))
			} else {
				thisStars <- sub( "^ +", "", thisLine)
			}
			starText <- paste( starText, thisStars, sep="")
			curLine <- curLine + 1
			next
		}
		if ( thisLine == "" || substr( thisLine, 1, 4) == "    " || curLine %in% gapLines) {
			curLine <- curLine + 1
			gptr <- 1
			next
		}

		thisGene <- sub( "(^.+)( +)(.+$)","\\1", thisLine)
		thisSeq <- sub( "(^.+)( +)(.+$)","\\3", thisLine)
		thisGene <- gsub( " ", "", thisGene)
		thisSeq <- gsub( " ", "", thisSeq)
		if (gptr > NG) { 
			NG <- gptr
			big[NG] <- gnames[NG] <- ""
		}
		gnames[gptr] <- thisGene
		big[gptr] <- paste( big[gptr], thisSeq, sep="")
		gptr <- gptr + 1
		curLine <- curLine + 1
		starHeadPtr <- c( starHeadPtr, regexpr( thisSeq, thisLine, fixed=T))
	}
	allLens <- nchar(big)
	if (verbose) cat( "\n", NG, gnames, "\nLen: ", allLens)

	NB <- max( allLens)
	bigM <- matrix("", nrow=NG, ncol=NB)
	for ( i in 1:NG) {
		thisLen <- allLens[i]
		bigM[ i, 1:thisLen] <- strsplit( big[i], split="")[[1]][1:thisLen]
	}
	rownames(bigM) <- gnames
	colnames(bigM) <- 1:NB

	starText <- strsplit( starText, split="")[[1]]
	tbl <- table( starText)
	whichStar <- match( "*", names(tbl), nomatch=0)
	pctConserved <- 0
	if (whichStar > 0) pctConserved <- tbl[whichStar] / sum(tbl)
	if (verbose) cat( "\nPct_Conserved: ", as.percent( pctConserved))

	return( list( "alignment"=bigM, "consensus"=starText, "pct.conserved"=pctConserved))
}


`writeALN` <- function( x, outfile, line.width=60, 
			title="CLUSTAL format alignment from DuffyTools object",
			capitalize=NULL, blankMissingFlanks=FALSE, fastaToo=FALSE) {

	# make sure this object has the elements we need
	if ( ! all( c( "alignment", "consensus") %in% names(x))) {
		cat( "\nError.  Expected an object with 'alignment' and 'consensus' components")
		return(NULL)
	}

	align <- x$alignment
	consensus <- x$consensus
	ids <- rownames(align)
	nIDs <- nrow(align)
	nCHAR <- ncol(align)

	con <- file( outfile, open="wt")
	writeLines( title, con=con)
	writeLines( "\n", con=con)

	# possibly alter lower/upper case
	if ( ! is.null( capitalize)) {
		if ( is.logical(capitalize) && capitalize == TRUE) capitalize <- "upper"
		if ( is.logical(capitalize) && capitalize == FALSE) capitalize <- "lower"
		caseChoices <- c( "upper", "lower")
		capitalize <- caseChoices[ pmatch( capitalize, caseChoices, nomatch=0)]
		if (capitalize == "upper") align <- toupper( align)
		if (capitalize == "lower") align <- tolower( align)
	}

	# we can allow blanking out of the end flanks when no data was present in the MSA.
	# so missing is not the same as true insertions/deletions
	if (blankMissingFlanks) {
		N <- nrow(align)
		NC <- ncol(align)
		for ( i in 1:N) {
			v <- align[ i, ]
			for( j in 1:NC) {
				if ( v[j] != "-") break
				v[j] <- " "
			}
			for( j in NC:1) {
				if ( v[j] != "-") break
				v[j] <- " "
			}
			align[ i, ] <- v
		}

		# once we blank flanks, reassess consensus calls
		for ( j in 1:NC) {
			v <- align[ , j]
			if ( length( unique( setdiff( v, " "))) == 1) consensus[j] <- "*"
		}
	}

	# step along, writing out a chunk at a time
	chDone <- 0
	repeat {
		if ( chDone >= nCHAR) break
		istart <- chDone + 1
		istop <- min( nCHAR, (chDone + line.width))

		for ( i in 1:nIDs) {
			thisFrag <- paste( align[ i, istart:istop], collapse="")
			thisID <- "               "
			substr( thisID, 1, nchar(ids[i])) <- ids[i]
			outText <- paste( thisID, thisFrag, sep=" ")
			writeLines( outText, con=con)
		}
		thisFrag <- paste( consensus[ istart:istop], collapse="")
		thisID <- "               "
		outText <- paste( thisID, thisFrag, sep=" ")
		writeLines( outText, con=con)
		writeLines( "", con=con)
		chDone <- istop
	}
	close( con)

	# perhaps write as FASTA too?
	if (fastaToo) {
		outfile <- paste( outfile, "fasta", sep=".")
		strs <- apply( align, MARGIN=1, FUN=paste, collapse="")
		writeFasta( as.Fasta( ids, strs), outfile, line.width=line.width)
	}
}


`ALNtoFasta` <- function( f, outfile=sub( "aln$", "fasta", f), gap.character="-", line.width=100, reorder=NULL) {

	# convert a .ALN Clustal or MAFFT MSA alignment result back to FASTA with the spacing intact
	aln <- readALN( f, verbose=F)

	ch <- aln$alignment
	desc <- rownames(ch)

	strs <- apply( ch, MARGIN=1, paste, collapse="")

	if ( gap.character != "-") strs <- gsub( "-", gap.character, strs, fixed=T)

	if ( ! is.null( reorder)) {
		# put the ones with the least gaps at the top
		if ( reorder == "N_Gaps") {
			nGaps <- apply( ch, MARGIN=1, function(x) sum( x == "-"))
			ord <- order( nGaps)
			strs <- strs[ord]
			desc <- desc[ord]
		} else {
			cat( "\nUnsupported 'Reorder' option: ", reorder)
			cat( "\nChoices are:  'N_Gaps' ")
			stop()
		}
	}

	# ready to write out
	writeFasta( as.Fasta( desc, strs), outfile, line.width=line.width)
}
