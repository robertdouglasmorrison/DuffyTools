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


`writeALN` <- function( x, outfile, line.width=60, max.id.width=30, 
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
	observed.id.width <- max( nchar(ids), na.rm=T)

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
	max.id.width.use <- min( observed.id.width, max.id.width)
	blankID <- paste( rep.int(" ",max.id.width.use), collapse="")
	
	repeat {
		if ( chDone >= nCHAR) break
		istart <- chDone + 1
		istop <- min( nCHAR, (chDone + line.width))

		for ( i in 1:nIDs) {
			thisFrag <- paste( align[ i, istart:istop], collapse="")
			thisID <- blankID
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


`plotALN` <- function( aln, type=c("fill","dots"), cex.letter=1, cex.label=1, y.label.length=NULL, 
			codonMap=getCodonMap(), ref.row=1, number.from=NULL, range=NULL, ...) {

	# we may be given the top level ALN object or just the aligment matrix
	# or even just the filename
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	}

	# OK, we have one matrix of characters from a MSA
	# allow subsetting on a range of locations
	# if we were given an explicit numbering start, then the range is in those units
	range.offset <- 0
	referenceRowChars <- aln[ ref.row, ]
	refGaps <- which( referenceRowChars == "-")
	if ( ! is.null(range) && ! is.null(number.from)) {
		range <- range - as.integer(number.from) + 1
		# we need to account for gaps when we set the true local range
		nGapLeft <- sum( refGaps <= min(range))
		nGapRight <- sum( refGaps <= max(range))
		range[1] <- min(range) + nGapLeft
		range[length(range)] <- max(range) + nGapRight
		number.from <- number.from + range[1] - 1
	}
	if( ! is.null( range)) {
		low <- max( 1, min(range))
		high <- min( ncol(aln), max(range))
		aln <- aln[ , low:high]
		range.offset <- low - 1
	}
	aln <- toupper(aln)
	nch <- ncol(aln)
	niso <- nrow(aln)
	plot( 1,1, type="n", xlim=c(0,nch+1), ylim=c(0,niso+1), xlab="Amino Acid Location", 
			xaxs="i", xaxt="n", ylab=NA, yaxs="i", yaxt="n", ...)
	if ( is.null( number.from)) {
		axis( side=1, ...)
		axis( side=3, ...)
	} else {
		# when given a explicit starting number, we need to do a bunch of math to keep that numbering system consistent
		number.from <- as.integer( number.from)
		tickStep <- 10
		shown.ats <-  unique( floor( ((number.from+tickStep):(number.from+nch+range.offset)) / tickStep) * tickStep)
		xAts <- shown.ats - number.from + 1
		# we need to count and deal with gaps in the reference row
		prevGaps <- 0
		for ( k in 1:length(xAts)) {
			nGapNow <- sum( refGaps <= (xAts[k] + range.offset + prevGaps))
			xAts[k] <- xAts[k] + nGapNow
			prevGaps <- nGapNow
		}
		# make the axis that shows the given numbering with all it's  gaps
		axis( side=1, at=xAts, label=as.character(shown.ats), ...)
		axis( side=3, at=xAts, label=as.character(shown.ats), ...)
	}

	# allow a few types of plot images
	type <- match.arg( type)
	if ( type == "fill") {
		boxes <- TRUE
		showLetters <- "all"
		letterColor <- "black"
	} else if (type == "dots") {
		boxes <- FALSE
		showLetters <- "alternate"
		letterColor <- "black"
		if (ref.row < 1) ref.row <- 1
		if (ref.row > niso) ref.row <- niso
	}

	cexScaleX <- 20 / niso
	cexScaleY <- 100 / nch
	cex <- min( cexScaleX, cexScaleY)

	# draw the alignment letters
	for( i in 1:niso) {
		x <- 1:nch
		y <- rep.int( niso - i + 1, length(x))
		ch <- aln[ i, ]
		where <- match( substr(ch,1,1), codonMap$AA, nomatch=0)
		if ( boxes) { 
			use <- which( where > 0)
			fill <- codonMap$Color[where]
			rect( x[use]-0.5, y[use]-0.5, x[use]+0.5, y[use]+0.5, col=fill, border=fill)
		}
		if ( showLetters == "all") {
			text( x, y, ch, col=letterColor, font=2, cex=cex*cex.letter)
		} else if (showLetters == "alternate") {
			chToShow <- rep.int( ".", nch)
			if ( i == ref.row) chToShow <- ch
			toShow <- which( ch != aln[ ref.row, ])
			chToShow[ toShow] <- ch[ toShow]
			# ?? for now let's leave all changed letters as black -- easier to see
			letterColor <- rep.int( "black", nch)
			toColor <- intersect( toShow, which( where > 0))
			letterColor[ toColor] <- codonMap$Color[where[toColor]]

			# let's use bold for what's different, but not for the dots
			isDot <- which( chToShow == ".")
			notDot <- setdiff( 1:nch, isDot)
			if ( length(isDot)) {
				text( x[isDot], y[isDot], chToShow[isDot], col=letterColor[isDot], font=1, cex=cex*cex.letter, adj=c(0.5,0))
			}
			if ( length(notDot)) {
				text( x[notDot], y[notDot], chToShow[notDot], col=letterColor[notDot], font=2, cex=cex*cex.letter, adj=c(0.5,0))
			}
		}
	}

	# draw the names.  Recall that we draw fist at the top...
	ylabs <- rev( rownames(aln))
	if ( ! is.null( y.label.length)) ylabs <- substr( ylabs, 1, y.label.length)
	axis( side=2, at=1:niso, label=ylabs, cex.axis=cex*cex.label, las=2)
	dev.flush()
}


`renumberALN` <- function( aln, number.from, ref.row=1) {

	# we may be given the top level ALN object or just the aligment matrix
	out <- aln
	isObject <- FALSE
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
		isObject <- TRUE
	}

	# OK, we have one matrix of characters from a MSA
	# see where the gaps are
	referenceRowChars <- aln[ ref.row, ]
	refGaps <- which( referenceRowChars == "-")

	# we need to count and deal with gaps in the reference row as we hand out new numbers
	NC <- ncol(aln)
	prevGaps <- 0
	newColNames <- rep.int( "", NC)
	for ( k in 1:NC) {
		nGapNow <- sum( refGaps <= (k))
		isGapNow <- (k %in% refGaps)
		if (isGapNow) {
			newColNames[k] <- paste( "G", as.character(prevGaps+1),sep="")
		} else {
			newColNames[k] <- as.character( number.from + k - 1 - nGapNow)
		}
		prevGaps <- nGapNow
	}
	
	colnames(aln) <- newColNames
	if (isObject) out$alignment <- aln
	return(out)
}


`summarizeALN` <- function( aln, range=NULL) {

	# we may be given the top level ALN object or just the aligment matrix
	# or even just the filename
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	}

	# OK, we have one matrix of characters from a MSA
	# allow subsetting on a range of locations
	if( ! is.null( range)) {
		low <- max( 1, min(range))
		high <- min( ncol(aln), max(range))
		aln <- aln[ , low:high]
	}
	aln <- toupper(aln)
	nch <- ncol(aln)
	niso <- nrow(aln)

	# tally up the various metrics
	consensCh <- diversePcts <- character( nch)
	consensPct <- double( nch)
	lapply( 1:nch, function(x) {
		chTbl <- sort( table( aln[,x]), decreasing=T)
		chSet <- names(chTbl)
		chCnts <- as.numeric(chTbl)
		chPcts <- round( chCnts * 100 / niso)
		consensCh[x] <<- chSet[1]
		consensPct[x] <<- chPcts[1]
		diversePcts[x] <<- paste( chSet, ":", chPcts, "%", sep="", collapse="; ")
		return(NULL)
	})

	# put the column names on these
	names(consensCh) <- colnames(aln)
	names(consensPct) <- colnames(aln)
	names(diversePcts) <- colnames(aln)
	out <- data.frame( "Location"=colnames(aln), "Consensus.Call"=consensCh, "Pct.Conserved"=consensPct, 
			"Diversity.Calls"=diversePcts, stringsAsFactors=F)
	return( out)
}

