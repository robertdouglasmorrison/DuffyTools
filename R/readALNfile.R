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
	blankConsensus <- blankID
	
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
		thisID <- blankConsensus
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
			codonMap=getCodonMap(), ref.row=1, number.from=NULL, range=NULL, max.X=NULL, 
			number.shown.offset=NULL, bottom.axis=TRUE, top.axis=TRUE, 
			xLabel="Amino Acid Location", main="MSA Alignment", ...) {

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
	nchScaling <- nch <- ncol(aln)
	niso <- nrow(aln)
	bigX <- nch + 1
	if ( ! is.null( max.X)) {
		bigX <- max.X
		nchScaling <- max.X
	}

	plot( 1,1, type="n", xlim=c(0,bigX), ylim=c(0.4,niso+0.6), xlab=NA, main=NA, 
			xaxs="i", xaxt="n", ylab=NA, yaxs="i", yaxt="n", ...)
	# place the labels and main text more precisely
	if ( ! is.na( main)) title( main=main, line=1.25, ...)
	title( xlab=xLabel, line=2.05, ...)

	if ( is.null( number.from)) {
		if (bottom.axis) axis( side=1, ...)
		if (top.axis) axis( side=3, mgp=c(3,0.5,0), ...)
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
		# allow the caller to adjust the shown numbering, to put it all on a full protein location etc.
		if ( ! is.null( number.shown.offset)) shown.ats <- shown.ats + number.shown.offset
		if (bottom.axis) axis( side=1, at=xAts, label=as.character(shown.ats), ...)
		if (top.axis) axis( side=3, at=xAts, label=as.character(shown.ats), mgp=c(3,0.5,0), ...)
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

	cexScaleY <- 20 / niso
	cexScaleX <- 120 / nchScaling
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
			# try to lighten them a bit
			fill <- adjustColor( fill, adjust=0.35)
			rect( x[use]-0.5, y[use]-0.5, x[use]+0.5, y[use]+0.5, col=fill, border=fill)
		}
		if ( showLetters == "all") {
			# allow a small offset to better center the text font
			text( x-0.05, y, ch, col=letterColor, font=2, cex=cex*cex.letter)
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


`plotALN.Panels` <- function( aln, type=c("fill","dots"), n.per.panel=100, cex.letter=1, cex.label=1, y.label.length=NULL, 
				codonMap=getCodonMap(), ref.row=1, number.from=1, range=NULL, 
				bottom.axis=TRUE, top.axis=FALSE, xLabel="Amino Acid Location", 
				letter.col=NULL,  main="MSA Alignment", mai=c( 0.42,1,0.42,0.2), ...) {

	# version for longer ALN sequences, where we do N letters per panel
	# we may be given the top level ALN object or just the aligment matrix
	# or even just the filename
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	}
	nch <- ncol(aln)

	# show numbering based on the reference sequence
	referenceRowChars <- aln[ ref.row, ]
	refNumbering <- cumsum( referenceRowChars != "-") + number.from - 1
	colnames(aln) <- refNumbering
	
	# see how many panels we need
	nPanels <- ceiling( nch / n.per.panel)
	if (nPanels > (par("fin")[2]/1.25)) cat( "\nWarning: may need too many panels for current plot window height..")
	savMF <- par( "mfrow")
	on.exit( par( mfrow=savMF), add=TRUE)
	par( mfrow=c( nPanels, 1))
	if ( ! is.null(mai)) {
		savMAI <- par( "mai")
		on.exit( par( mai=savMAI), add=TRUE)
		par( mai=mai)
	}
		
	type <- match.arg( type)
	nDone <- 0
	nowNumberFrom <- number.from
	lastRow <- FALSE
	while (nDone < nch) {
		# get the bounds of the next panel to show
		nFrom <- nDone + 1
		nTo <- nDone + n.per.panel
		if ( nTo > nch) {
			lastRow <- TRUE
			max.X <- n.per.panel
			nTo <- nch
		} else {
			max.X <- NULL
		}
		smlALN <- aln[ , nFrom:nTo, drop=F]
		nowNumberFrom <- as.numeric( colnames(smlALN)[1])

		# do this chunk, only show the title on the top one
		mainText <- if (nDone < 1) main else NA
		plotALN( smlALN, type=type, codonMap=codonMap, number.from=nowNumberFrom, main=mainText,
						max.X=max.X, etter.col=letter.col, cex.letter=cex.letter, cex.label=cex.label,  
						ref.row=ref.row, bottom.axis=bottom.axis, top.axis=top.axis, xLabel=xLabel, 
						letter.col=letter.col, ...)
						
		# increment
		nDone <- nDone + n.per.panel
	}
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


`ALN2DistanceMatrix` <- function( aln, gapAdjust=FALSE) {

	# convert a sequence alignment into a Distance Matrix (i.e. how dissimilar)
	# by counting up the number of disagreements
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	} else {
		cat( "\nError: expected an ALN 'alignment' element.")
		return(NULL)
	}

	# rownames are the sequence discriptors
	nams <- rownames(aln)
	N <- nrow(aln)
	dm <- matrix( 0, nrow=N, ncol=N)
	rownames(dm) <- colnames(dm) <- nams
	
	# with all MSA work finished, just count up differences
	for ( i in 1:(N-1)) {
		ch1 <- aln[ i, ]
		if ( gapAdjust) isGap1 <- which( ch1 == "-")
		for ( j in (i+1):N) {
			ch2 <- aln[ j, ]
			ndif <- sum( ch1 != ch2)
			if (gapAdjust) {
				isGap2 <- which( ch2 == "-")
				# locations where one (but not both) were gaps, get only half the penalty
				only1 <- setdiff( isGap1, isGap2)
				only2 <- setdiff( isGap2, isGap1)
				nAdj <- (length(only1) + length(only2))/2
				ndif <- ndif - nAdj
			}
			dm[i,j] <- dm[j,i] <- ndif
		}
	}
	return(dm)
}


`plotALN.BitScore` <- function( aln, heightM=NULL, codonMap=getCodonMap(), ref.row=1, number.from=1, 
				max.X=NULL, max.Y=NULL, letter.col=NULL, min.bit.score=0.01, main="Sequence Logo",
				xLabel="Amino Acid Location (NF54)", yLabel="Bit Score", gap.x=0.08, gap.y=0.05, 
				col.axis=1, col.lab=1, ...) {

	# we may be given the top level ALN object or just the aligment matrix
	# or even just the filename
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	}

	# OK, we have one matrix of characters from a MSA
	
	# turn this into frequency matrix, Infomation Content, and letter heights
	icAns <- ALNtoInformationContent( aln)
	heightM <- icAns$height
	
	# show numbering based on the reference sequence
	referenceRowChars <- aln[ ref.row, ]
	refNumbering <- cumsum( referenceRowChars != "-") + number.from - 1
	colnames(aln) <- colnames(icAns$proportion) <- refNumbering
	
	# set up to plot letters scaled by info content, using a polygon database of a system font
	# data frame is called "fontPolygons", with columns: x, y, letter, order
	data( "Helvetica.Regular.Font.Polygons")
	nch <- ncol(aln)
	bigX <- number.from + nch - 1
	if ( ! is.null( max.X)) bigX <- max( bigX, max.X)
	yBig <- max( apply(heightM,2,sum,na.rm=T)) + (gap.y * 3) # allow for a bit extra between letters
	if ( ! is.null( max.Y)) yBig <- max.Y

	# place the labels and main text more precisely
	plot( 1,1, type="n", xlim=c(number.from, bigX), ylim=c(0,yBig), 
			xlab=NA, ylab=yLabel, frame.plot=F, main=NA, xaxt="n", yaxt="n", col.lab=col.lab, ...)
	if ( ! is.na( main)) title( main=main, line=1.25, col.main=col.lab)
	title( xlab=xLabel, line=2.05, col.lab=col.lab)
	axis( side=2, col=col.axis, col.ticks=col.axis, col.axis=col.axis, ...)
	
	# for the axis numbers, try to account for gaps, and try to make the range fill the actual sequence
	prettyXpts <- intersect( pretty(refNumbering), refNumbering)
	nAxisPts <- length(prettyXpts)
	if ( (prettyXpts[1] - refNumbering[1]) > 5) {
		prettyXpts <- c( refNumbering[1], prettyXpts)
		nAxisPts <- length(prettyXpts)
	}
	if ( ( refNumbering[nch] - prettyXpts[nAxisPts]) > 5) {
		prettyXpts <- c( prettyXpts, refNumbering[nch])
	}
	prettyXats <- match( prettyXpts, refNumbering) + refNumbering[1] - 1
	axis( side=1, at=prettyXats, label=prettyXpts, col=col.axis, col.ticks=col.axis, col.axis=col.axis, ...)

	# draw the alignment letters, with size proportional to abundance
	# set up the gapping details
	xScale <- (1 - gap.x)
	xScale2 <- xScale / 2
	yScale <- (1 - gap.y)
	for( i in 1:nch) {
		x <- i + number.from - 1
		# only draw the letters with non-zero height
		myHts <- heightM[ , i]
		show <- which( myHts > 0)
		# show them with the best at the top
		show <- show[ order( myHts[show])]
		# now visit each seen AA, and draw it proportional to it's abundance
		ynow <- 0
		for (k in show) {
			myCh <- rownames( heightM)[k]
			myHt <- myHts[k]
			# don't draw if too small to see
			if ( myHt < min.bit.score) next
			myCol <- codonMap$Color[ match( myCh, codonMap$AA)]
			if ( ! is.null( letter.col)) myCol <- letter.col[ match( myCh, names(letter.col))]
			if ( is.na(myCol)) myCol <- 'black'
			# get the polygon data for this letter
			smlFont <- subset( fontPolygons, letter == myCh)
			# the font is alway 0 to 1 in X and Y.  Shift to center X on the point, 
			# shrink X a bit to avoid touching, and scale the Y by Info Content
			xx <- (smlFont$x * xScale) + x - xScale2
			yy <- (smlFont$y * myHt)
			polygon( xx, (yy*yScale) + ynow, col=myCol, border=myCol)
			ynow <- ynow + max(yy) + gap.y
		}
	}
	dev.flush()
}


`plotALN.BitScore.Panels` <- function( aln, n.per.panel=100, codonMap=getCodonMap(), ref.row=1, number.from=1, 
					max.X=NULL, max.Y=NULL, letter.col=NULL, min.bit.score=0.01, main="Sequence Logo", 
					xLabel="Amino Acid Location (NF54)", gap.x=0.12, gap.y=0.1, 
					mai=c( 0.42,1,0.42,0.2), ...) {

	# version for longer ALN sequences, where we do N letters per panel
	# we may be given the top level ALN object or just the aligment matrix
	# or even just the filename
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	}
	nch <- ncol(aln)

	# show numbering based on the reference sequence
	referenceRowChars <- aln[ ref.row, ]
	refNumbering <- cumsum( referenceRowChars != "-") + number.from - 1
	colnames(aln) <- refNumbering
	
	# see how many panels we need
	nPanels <- ceiling( nch / n.per.panel)
	if (nPanels > (par("fin")[2]/1.25)) cat( "\nWarning: may need too many panels for current plot window height..")
	savMF <- par( "mfrow")
	on.exit( par( mfrow=savMF), add=TRUE)
	par( mfrow=c( nPanels, 1))
	if ( ! is.null(mai)) {
		savMAI <- par( "mai")
		on.exit( par( mai=savMAI), add=TRUE)
		par( mai=mai)
	}
	
	# calculate the Information Content on the full sequence
	icAns <- ALNtoInformationContent( aln)
	colnames(icAns$proportion) <- refNumbering
	heightM <- icAns$height
	bigY <- max( heightM)
	if ( ! is.null( max.Y)) bigY <- max( bigY, max.Y)
	
	nDone <- 0
	nowNumberFrom <- number.from
	lastRow <- FALSE
	while (nDone < nch) {
		# get the bounds of the next panel to show
		nFrom <- nDone + 1
		nTo <- nDone + n.per.panel
		if ( nTo > nch) {
			lastRow <- TRUE
			nTo <- nch
		} else {
			max.X <- NULL
		}
		smlALN <- aln[ , nFrom:nTo, drop=F]
		smlHT <- heightM[ , nFrom:nTo, drop=F]
		nowNumberFrom <- as.numeric( colnames(smlALN)[1])
		if ( lastRow) max.X <- nowNumberFrom + n.per.panel - 1

		# do this chunk, only show the title on the top one
		mainText <- if (nDone < 1) main else NA
		plotALN.BitScore( smlALN, heightM=smlHT, codonMap=codonMap, number.from=nowNumberFrom, main=mainText,
				max.X=max.X, max.Y=bigY, letter.col=letter.col, min.bit.score=min.bit.score, 
				gap.x=gap.x, gap.y=gap.y, ...)
						
		# increment
		nDone <- nDone + n.per.panel
	}
}
