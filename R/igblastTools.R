# igblastTools.R -- various code for working with IgBlast


`readIgBlastOutput` <- function( infile="igblastOut.txt", m=7, verbose=TRUE, min.blast.score=40) {

	# scan the output of Blast for the data columns we need
	# depending on the version of BLAST, the column mode number is different
	if ( m != 7) {
		warning( paste( "Unsupported IgBlast Output format...:  m=", m))
		return( NULL)
	}

	# the file has lots of comments, other details, etc.
	# need to parse by hand
	blastText <- readLines( infile)
	if ( length( blastText) < 6) return(NULL)

	keyLines <- grep( "^# IGBLAST", blastText)
	nKeys <- length( keyLines)
	keyEnds <- c( keyLines[2:nKeys] - 1, length(blastText))
	
	nout <- 0
	idOut <- vNamOut <- dNamOut <- jNamOut <- vector()
	vScoreOut <- dScoreOut <- jScoreOut <- vector()
	pctIdent <- vector()


	selectBestVDJ <- function( seg, nam, pct, eval, score) {

		vn <- dn <- jn <- ""
		vs <- ds <- js <- 0
		vp <- dp <- jp <- 0

		seg <- toupper( seg)
		isV <- which( seg == "V")
		if ( length(isV)) {
			best <- which.max( score[isV])
			vn <- nam[ isV[ best]]
			vs <- as.numeric( score[ isV[ best]])
			vp <- as.numeric( pct[ isV[ best]])
		}
		isD <- which( seg == "D")
		if ( length(isD)) {
			best <- which.max( score[isD])
			dn <- nam[ isD[ best]]
			ds <- as.numeric( score[ isD[ best]])
			dp <- as.numeric( pct[ isD[ best]])
		}
		isJ <- which( seg == "J")
		if ( length(isJ)) {
			best <- which.max( score[isJ])
			jn <- nam[ isJ[ best]]
			js <- as.numeric( score[ isJ[ best]])
			jp <- as.numeric( pct[ isJ[ best]])
		}

		bestScore <- max( vs, ds, js)
		pctIdent <- mean( c(vp, dp, jp)[ which( c(vs, ds, js) > 0)])

		return( list( "Vname"=vn, "Vscore"=vs, "VpctIdent"=vp, 
				"Dname"=dn, "Dscore"=ds, "DpctIdent"=dp, 
				"Jname"=jn, "Jscore"=js, "JpctIdent"=jp, 
				"BestScore"=bestScore, "PctIdentity"=pctIdent))
	}
	

	mapply( keyLines, keyEnds, FUN=function( start, stop) {

			myID <- sub( "# Query: ", "", blastText[ start+1])

			# the rearangement facts...
			vmatch <- jmatch <- ""
			rearange <- grep( "rearangement summary", blastText[ start:stop], fixed=T)
			#cat( "\nDebug: ", start, stop, myID, rearange)

			if ( ! is.na( who <- rearange[1])) {
				terms <- strsplit( blastText[ who+start], split="\t")[[1]]
				vmatch <- terms[1]
				jmatch <- terms[2]
			}

			# the blast hits for V,D,J
			hits <- grep( "Hit table", blastText[ start:stop], fixed=T)
			bestScore <- 0
			if ( ! is.na( who <- hits[1])) {
				n <- 0
				seg <- nam <- pct <- eval <- score <- vector()
				for (k in (who+start+2):stop) {
					terms <- strsplit( blastText[k], split="\t")[[1]]
					n <- n + 1
					seg[n] <- terms[1]
					nam[n] <- terms[3]
					pct[n] <- terms[4]
					eval[n] <- terms[13]
					score[n] <- terms[14]
				}
				ans <- selectBestVDJ( seg, nam, pct, eval, score)
				bestScore <- ans$BestScore
			} 

			if ( bestScore < min.blast.score) return()
			nout <<- nout + 1
			idOut[nout] <<- myID
			vNamOut[nout] <<- ans$Vname
			vScoreOut[nout] <<- ans$Vscore
			dNamOut[nout] <<- ans$Dname
			dScoreOut[nout] <<- ans$Dscore
			jNamOut[nout] <<- ans$Jname
			jScoreOut[nout] <<- ans$Jscore
			pctIdent[nout] <<- ans$PctIdentity
			#cat( "\t", unlist( ans))
			return()

		})

	blastDF <- data.frame( "CONTIG_ID"=idOut, "V_NAME"=vNamOut, "D_NAME"=dNamOut, "J_NAME"=jNamOut,
				 "V_SCORE"=vScoreOut, "D_SCORE"=dScoreOut, "J_SCORE"=jScoreOut,
				 "PctIdentity"=pctIdent, stringsAsFactors=FALSE)
	ttlScore <- blastDF$V_SCORE + blastDF$D_SCORE + blastDF$J_SCORE
	ord <- order( ttlScore, decreasing=T)
	blastDF <- blastDF[ ord, ]
	rownames(blastDF) <- 1:nrow(blastDF)

	# clean up and extract the details we need for later profiling
	blastDF$V_GROUP <- igblastGroupName( blastDF$V_NAME)
	blastDF$J_GROUP <- igblastGroupName( blastDF$J_NAME)
	blastDF$READ_DEPTH <- round( as.numeric( sub( "^NODE.+_cov_", "",  blastDF$CONTIG_ID)), digits=3)
	blastDF$PctIdentity <- round( as.numeric(  blastDF$PctIdentity), digits=3)

	return( blastDF)
}



`callIgBlast` <- function( fastafile, outfile="igblastOut.txt", program="igblastn", db=c("Both", "IG", "TR"),
			path="~/IgBlast", organism=c( "human", "mouse"), outfmt=7) {

	# validate the arguments
	useprogram <- file.path( path, program)
	if ( ! file.exists(useprogram)) stop( paste( "Can't find IgBlast executable: ", useprogram))
	cmdline <- useprogram 

	# push the directory of IgBlast databases to a system env variable
	path <- file.path( dirname(path), basename(path))
	Sys.setenv( "IGDATA"=path)

	organism <- match.arg( organism)
	db <- match.arg( db)
	if ( db == "Both") db <- "IGTR"
	dbArgs <- paste( " -germline_db_V ", path, "/IMGT_", organism, "_", db, "_Vdb", 
			" -germline_db_D ", path, "/IMGT_", organism, "_", db, "_Ddb", 
			" -germline_db_J ", path, "/IMGT_", organism, "_", db, "_Jdb", sep="")
	cmdline <- paste( cmdline, dbArgs)

	if ( db == "TR") cmdline <- paste( cmdline, " -ig_seqtype TCR ")

	cmdline <- paste( cmdline, " -organism ", organism, " -domain_system kabat ")
	cmdline <- paste( cmdline, " -auxiliary_data ", file.path( path, paste( organism, "_gl.aux", sep="")))
	cmdline <- paste( cmdline, " -num_threads ", 4)

	if ( ! file.exists( fastafile)) stop( "Can't find FASTA query file: ", fastafile)
	cmdline <- paste( cmdline, "  -query ", fastafile)
	cmdline <- paste( cmdline, "  -outfmt ", outfmt, "  > ", outfile)

	# call BLAST
	cat( "\nCalling IgBLAST: \nCommand line:  ", cmdline, "\n")
	system( cmdline)
	cat( "\nDone.\nWrote file: ", outfile, "\n")
	return()
}


`igblastGroupName` <- function( txt) {

	# the construct names have detailed suffixes that we want to lop off
	txt <- sub( "-.+", "", txt)
	txt <- sub( "\\*.+", "", txt)
	txt <- sub( "/.+", "", txt)
	#txt <- sub( "P$", "", txt)
	#txt <- sub( "P", "", txt)
	txt <- sub( "D$", "", txt)

	# some may be unusual...
	txt <- sub( "(I)", "1", txt, fixed=T)
	txt <- sub( "(II)", "2", txt, fixed=T)
	txt <- sub( "(III)", "3", txt, fixed=T)
	txt <- sub( "(IV)", "4", txt, fixed=T)
	txt <- sub( "(V)", "5", txt, fixed=T)
	txt <- sub( "(VI)", "6", txt, fixed=T)
	txt <- sub( "(VII)", "7", txt, fixed=T)

	# add some space between the Class and the variant
	txt <- paste( substr( txt, 1,3), substr( txt, 4, 8), sep="_")

	txt
}


`igblastVnames` <- function( mode=c( "IG", "TR")) {

	mode <- match.arg( mode)
	if ( mode == "IG") {
		return( c( paste( "IGH_V", 1:7, sep=""), 
				paste( "IGK_V", 1:7, sep=""),
				paste( "IGL_V", 1:9, sep="")))
	}
	if ( mode == "TR") {
		return( c( paste( "TRA_V", 1:41, sep=""), 
				paste( "TRB_V", 1:30, sep=""), "TRB_VA", "TRB_VB", "TRB_VC",
				paste( "TRD_V", 1:3, sep=""), 
				paste( "TRG_V", 1:11, sep="")))
	}
}


`igblastJnames` <- function( mode=c( "IG", "TR")) {

	mode <- match.arg( mode)
	if ( mode == "IG") {
		return( c( paste( "IGH_J", 1:6, sep=""), 
				paste( "IGK_J", 1:5, sep=""),
				paste( "IGL_J", 1:7, sep="")))
	}
	if ( mode == "TR") {
		return( c( paste( "TRA_J", 1:61, sep=""), 
				paste( "TRB_J", 1:2, sep=""),
				paste( "TRD_J", 1:4, sep=""), 
				paste( "TRG_J", 1:2, sep=""), "TRG_JP", "TRG_JP1", "TRG_JP2"))
	}
}


`profileIgBlastCoverage` <- function( tbl, min.blast.score=100) {

	# eat the IgBlast result data frame, and present it as read depth percentages and 
	# germline identity per lymphocyte construct
	tbl <- subset( tbl, V_SCORE >= min.blast.score)
	drops <- union( which( tbl$V_NAME == ""), which( tbl$J_NAME == ""))
	if ( length(drops)) tbl <- tbl[ -drops, ]
	if ( nrow(tbl) < 1) return( NULL)

	# extract the details we need
	vnam <- tbl$V_GROUP
	jnam <- tbl$J_GROUP
	readDepth <- tbl$READ_DEPTH
	pctIdent <- tbl$PctIdentity

	# make a table to accumulate the reads by group
	# use all possible names, not just those seen in this file
	mode <- substr( vnam[1], 1, 2)
	vlevels <- igblastVnames( mode)
	jlevels <- igblastJnames( mode)
	vfac <- factor( vnam, levels=vlevels)
	jfac <- factor( jnam, levels=jlevels)
	nv <- length( vlevels)
	nj <- length( jlevels)
	m <- n <- p <- matrix( 0, nrow=nj, ncol=nv)
	colnames(m) <- colnames(p) <- vlevels
	rownames(m) <- rownames(p) <- jlevels

	for ( k in 1:nrow(tbl)) {
		i <- jfac[k]
		j <- vfac[k]
		m[ i, j] <- m[i,j] + readDepth[k]
		n[ i, j] <- n[i,j] + 1
		p[ i, j] <- p[i,j] + pctIdent[k]
	}

	# average the percent Identity data
	for (i in 1:ncol(p)) {
		who <- which( n[ ,i] > 0)
		p[ who, i] <- p[ who, i] / n[ who, i]
	}
	totalReads <- sum( m)
	m <- m * 100 / totalReads

	return( list( "profile"=m, "identity"=p))
}


`plotIgBlastProfile` <- function( profile, sampleID="", label="", sublabel=NULL, min.read.pct=0.1, ...) {

	# make a 'Vlad-style' color plot of circles where size is percentage of reads and
	# color is deviation from germline...
	require( plotrix)

	pcts <- profile$profile
	ident <- profile$identity

	# lets drop rows/columns with 'not enough data'
	rowMax <- apply( pcts, 1, max)
	drops <- which( rowMax < min.read.pct)
	if ( length( drops)) {
		pcts <- pcts[ -drops, ]
		ident <- ident[ -drops, ]
	}
	colMax <- apply( pcts, 2, max)
	drops <- which( colMax < min.read.pct)
	if ( length( drops)) {
		pcts <- pcts[ , -drops]
		ident <- ident[ , -drops]
	}

	vIdent <- as.vector( ident)
	vIdent[ vIdent < 90] <- 90
	colorScale <- color.scale( vIdent, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1),
				xrange=c(90,100))

	NJ <- nrow(pcts)
	NV <- ncol(pcts)
	xlocs <- matrix( rep( 1:NV, each=NJ), NJ, NV)
	ylocs <- matrix( rep( 1:NJ, times=NV), NJ, NV)

	# put a minimum on the size for plotting
	NJ <- max( NJ, 7)
	NV <- max( NV, 9)

	# use the function that does it all at once
	# we may have more than one sample
	if ( length( sampleID) > 1) {
		textID <- paste( "AVG( ", paste( sampleID, collapse="; "), " )", sep="")
	} else {
		textID <- sampleID[1]
	}
	lymphKey <- names( sort( table( substr( c( colnames(pcts), rownames(pcts)), 1, 2)), decreasing=T))[1]
	if ( lymphKey == "TR") lymphKey <- "TCR"

	if ( is.null( sublabel)) {
		mainText <- paste( lymphKey, " Profile:    ", label, "\nSampleID:  ", textID)
	} else {
		mainText <- paste( lymphKey, " Profile:    ", label, "\n", sublabel)
	}

	my_size_n_color( x=xlocs, y=ylocs, size=pcts/100, col=colorScale,
			xaxlab=colnames(pcts), yaxlab=rownames(pcts), xcex=1.2, ycex=1.2,
			xgrid=T, ygrid=T, xlim=c(0.5,NV+0.5), ylim=c(-1.7,NJ+0.5), font.lab=2,
			xlas=3, main=mainText, mar=c( 6,6,4,2), ...)

	lines( c( -1,NV+2), c(0,0), lty=1, lwd=2, col=1)
	colorscale <- color.scale( 90:100, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1), xrange=c(90,100))
	colorXlo <- round( NV * 0.33)
	colorXhi <- NV - colorXlo + 1
	color.legend( colorXlo, -1.6, colorXhi, -0.9, seq(90, 100, by=1), rect.col=colorscale)
	text( colorXlo, -1.25, "GermLine Identity", pos=2, cex=1.2, font=2)

	sizeX <- round( NV * 0.90)
	draw.circle( sizeX, -1.3, sqrt(0.05), col="dodgerblue", border=1)
	text( sizeX, -1.01, "5% of Constructs", pos=3, font=2, cex=1.15)
}


`getIgBlastProfiles` <- function( profileFile) {

	# the profile is stored as a cbind of percentages then germline identity
	tbl <- read.delim( profileFile[1], as.is=T)
	N2 <- ncol(tbl)
	N <- as.integer( N2/2)
	profile <- as.matrix( tbl[ ,1:N, drop=FALSE])
	pctIdent <- as.matrix( tbl[ ,(N+1):N2, drop=FALSE])
	ans <- list( "profile"=profile, "identity"=pctIdent)

	# we may have be given a set of files, is which case we merge
	NP <- length(profileFile)
	if ( NP == 1) return( ans)

	profList <- identList <- vector( mode="list")
	profList[[1]] <- ans$profile
	identList[[1]] <- ans$identity
	for ( j in 2:NP) {
		tbl <- read.delim( profileFile[1], as.is=T)
		N2 <- ncol(tbl)
		N <- as.integer( N2/2)
		profile <- as.matrix( tbl[ ,1:N, drop=FALSE])
		pctIdent <- as.matrix( tbl[ ,(N+1):N2, drop=FALSE])
		profList[[j]] <- profile
		identList[[j]] <- pctIdent
	}

	# first decide the union for both D's
	newcols <- colnames( profList[[1]])
	newrows <- rownames( profList[[1]])
	for ( j in 2:NP) {
		colnam2 <- colnames( profList[[j]])
		newcols <- sort( union( newcols, colnam2))
		rownam2 <- rownames( profList[[j]])
		newrows <- sort( union( newrows, rownam2))
	}
	NC <- length(newcols)
	NR <- length(newrows)

	# make new bigger matrices to hold the average over all samples
	newpct <- newident <- matrix( 0, nrow=NR, ncol=NC)
	colnames(newpct) <- colnames(newident) <- newcols
	rownames(newpct) <- rownames(newident) <- newrows

	whRow2 <- match( newrows, rownam2, nomatch=0)

	# visit all cells, and get whoever has data for that cell
	for ( j in 1:NC) {
		myCol <- newcols[j]
		for ( i in 1:NR) {
			myRow <- newrows[i]
			vProf <- vIdent <- vector()
			for (k in 1:NP) {
				whereC <- match( myCol, colnames(profList[[k]]), nomatch=0)
				whereR <- match( myRow, rownames(profList[[k]]), nomatch=0)
				if ( any( c( whereC, whereR) == 0)) next
				vProf <- c( vProf, profList[[k]][ whereR, whereC])
				vIdent <- c( vIdent, identList[[k]][ whereR, whereC])
			}
			if ( length( vProf) == 0) next
			# store the average of what we found
			newpct[ i, j] <- mean( vProf, na.rm=T)
			newident[ i, j] <- mean( vIdent, na.rm=T)
		}
	}

	ans <- list( "profile"=newpct, "identity"=newident)


}


`diffIgBlastProfiles` <- function( profile1, profile2) {

	# given the IgBlast profiles from 2 time points, calculate the difference
	pct1 <- profile1$profile
	ident1 <- profile1$identity
	pct2 <- profile2$profile
	ident2 <- profile2$identity

	# we need to expand each to the union of both sets, before we can do the difference

	# first decide the union for both D's
	colnam1 <- colnames( pct1)
	colnam2 <- colnames( pct2)
	newcols <- sort( union( colnam1, colnam2))
	NC <- length(newcols)
	rownam1 <- rownames( pct1)
	rownam2 <- rownames( pct2)
	newrows <- sort( union( rownam1, rownam2))
	NR <- length(newrows)

	newpct1 <- newident1 <- newpct2 <- newident2 <- matrix( 0, nrow=NR, ncol=NC)
	colnames(newpct1) <- colnames(newident1) <- newcols
	colnames(newpct2) <- colnames(newident2) <- newcols
	rownames(newpct1) <- rownames(newident1) <- newrows
	rownames(newpct2) <- rownames(newident2) <- newrows

	# now decide where each new entry is in the original data
	whCol1 <- match( newcols, colnam1, nomatch=0)
	whCol2 <- match( newcols, colnam2, nomatch=0)
	whRow1 <- match( newrows, rownam1, nomatch=0)
	whRow2 <- match( newrows, rownam2, nomatch=0)

	# lastly, move the data that exists to its new location
	for ( i in 1:NC) {
		if ( whCol1[i] > 0) {
			newpct1[ whRow1 > 0, i] <- pct1[ , whCol1[i]]
			newident1[ whRow1 > 0, i] <- ident1[ , whCol1[i]]
		}
		if ( whCol2[i] > 0) {
			newpct2[ whRow2 > 0, i] <- pct2[ , whCol2[i]]
			newident2[ whRow2 > 0, i] <- ident2[ , whCol2[i]]
		}
	}

	# now the subtraction is trivial
	newpct <- newpct2 - newpct1
	newident <- newident2 - newident1

	return( list( "diffProfile"=newpct, "diffIdentity"=newident))
}


`plotIgBlastDiffProfile` <- function( diffProfile, sampleID1="sample1", sampleID2="sample2", 
				label="", sublabel=NULL, min.read.pct=0.1, ...) {

	# make a 'Vlad-style' color plot of circles where size is percentage of reads and
	# color is deviation from germline...
	require( plotrix)

	pcts <- diffProfile$diffProfile
	ident <- diffProfile$diffIdentity

	# lets drop rows/columns with 'no data'
	rowMax <- apply( abs(pcts), 1, max)
	drops <- which( rowMax < min.read.pct)
	if ( length( drops)) {
		pcts <- pcts[ -drops, ]
		ident <- ident[ -drops, ]
	}
	colMax <- apply( abs(pcts), 2, max)
	drops <- which( colMax < min.read.pct)
	if ( length( drops)) {
		pcts <- pcts[ , -drops]
		ident <- ident[ , -drops]
	}

	vIdent <- as.vector( ident)
	vIdent[ vIdent < -15] <- -15
	vIdent[ vIdent > 15] <- 15
	colorScale <- color.scale( vIdent, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1),
				xrange=c(-15,15))

	NJ <- nrow(pcts)
	NV <- ncol(pcts)
	xlocs <- matrix( rep( 1:NV, each=NJ), NJ, NV)
	ylocs <- matrix( rep( 1:NJ, times=NV), NJ, NV)

	# put a minimum on the size for plotting
	NJ <- max( NJ, 7)
	NV <- max( NV, 9)

	# for now, any group that decreased its profile, we won't show...
	#pcts[ pcts < 0] <- 0
	# try no fill for the ones that shrank...
	isNeg <- (pcts < 0)
	colorNeg <- colorScale
	colorNeg[ ! isNeg] <- NA
	colorScale[ isNeg] <- NA
	pcts[ isNeg] <- abs( pcts[ isNeg])

	# may be given multiple IDs per group...
	if ( length( sampleID1) > 1) {
		textID1 <- paste( "AVG( ", paste( sampleID1, collapse="; "), " )", sep="")
	} else {
		textID1 <- sampleID1[1]
	}
	if ( length( sampleID2) > 1) {
		textID2 <- paste( "AVG( ", paste( sampleID2, collapse="; "), " )", sep="")
	} else {
		textID2 <- sampleID2[1]
	}
	lymphKey <- names( sort( table( substr( c( colnames(pcts), rownames(pcts)), 1, 2)), decreasing=T))[1]
	if ( lymphKey == "TR") lymphKey <- "TCR"
	if ( is.null( sublabel)) {
		mainText <- paste( "Change in ", lymphKey, " Profile:     ", label, "\n", textID1, "  ->  ", textID2)
	} else {
		mainText <- paste( "Change in ", lymphKey, " Profile:     ", label, "\n", sublabel)
	}

	# use the function that does it all at once
	my_size_n_color( x=xlocs, y=ylocs, size=pcts/100, col=colorScale, border=colorNeg,
			xaxlab=colnames(pcts), yaxlab=rownames(pcts), xcex=1.2, ycex=1.2,
			xgrid=T, ygrid=T, xlim=c(0.5,NV+0.5), ylim=c(-1.7,NJ+0.5), font.lab=2,
			xlas=3, main=mainText, mar=c( 6,6,4,2), ...)

	lines( c( -1,NV+2), c(0,0), lty=1, lwd=2, col=1)
	colorscale <- color.scale( 90:100, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1), xrange=c(90,100))
	colorXlo <- round( NV * 0.33)
	colorXhi <- NV - colorXlo + 1
	color.legend( colorXlo, -1.6, colorXhi, -0.9, seq( -15, 15, by=3), rect.col=colorscale)
	text( (colorXlo+colorXhi)/2, -1.5, "Change in GermLine", pos=1, cex=1.05, font=2)

	sizeX <- round( NV * 0.10)
	draw.circle( sizeX, -1.3, sqrt(0.05), border="dodgerblue", col=NA)
	text( sizeX, -1.01, "Decrease 5%", pos=3, font=2, cex=1.1)
	sizeX <- round( NV * 0.90)
	draw.circle( sizeX, -1.3, sqrt(0.05), col="dodgerblue", border=1)
	text( sizeX, -1.01, "Increase 5%", pos=3, font=2, cex=1.1)
}


`plotIgBlastOverlayProfile` <- function( profile1, profile2, sampleID1="sample1", sampleID2="sample2", 
				label="", min.read.pct=0.1, ...) {

	# make a 'Vlad-style' color plot of circles where size is percentage of reads and
	# color is deviation from germline...
	require( plotrix)

	pcts1 <- profile1$profile
	ident1 <- profile1$identity
	pcts2 <- profile2$profile
	ident2 <- profile2$identity

	# lets drop rows/columns with 'no data'
	rowMax1 <- apply( abs(pcts1), 1, max)
	drops1 <- which( rowMax1 < min.read.pct)
	rowMax2 <- apply( abs(pcts2), 1, max)
	drops2 <- which( rowMax2 < min.read.pct)
	drops <- intersect( drops1, drops2)
	if ( length( drops)) {
		pcts1 <- pcts1[ -drops, ]
		ident1 <- ident1[ -drops, ]
		pcts2 <- pcts2[ -drops, ]
		ident2 <- ident2[ -drops, ]
	}
	colMax1 <- apply( abs(pcts1), 2, max)
	drops1 <- which( colMax1 < min.read.pct)
	colMax2 <- apply( abs(pcts2), 2, max)
	drops2 <- which( colMax2 < min.read.pct)
	drops <- intersect( drops1, drops2)
	if ( length( drops)) {
		pcts1 <- pcts1[ , -drops]
		ident1 <- ident1[ , -drops]
		pcts2 <- pcts2[ , -drops]
		ident2 <- ident2[ , -drops]
	}

	# we will draw #2 as normal, then overlay #1 to show change
	vIdent <- as.vector( ident2)
	vIdent[ vIdent < 90] <- 90
	colorScale <- color.scale( vIdent, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1), xrange=c(90,100))

	NJ <- nrow(pcts2)
	NV <- ncol(pcts2)
	xlocs <- matrix( rep( 1:NV, each=NJ), NJ, NV)
	ylocs <- matrix( rep( 1:NJ, times=NV), NJ, NV)

	# put a minimum on the size for plotting
	NJ <- max( NJ, 7)
	NV <- max( NV, 9)

	if ( length( sampleID1) > 1) {
		textID1 <- paste( "AVG( ", paste( sampleID1, collapse="; "), " )", sep="")
	} else {
		textID1 <- sampleID1[1]
	}
	if ( length( sampleID2) > 1) {
		textID2 <- paste( "AVG( ", paste( sampleID2, collapse="; "), " )", sep="")
	} else {
		textID2 <- sampleID2[1]
	}
	mainText <- paste( "Change in IG/TCR Profile:     ", label, "\n", textID1, "  ->  ", textID2)

	# use the function that does it all at once
	my_size_n_color( x=xlocs, y=ylocs, size=pcts2/100, col=colorScale, border=NA,
			xaxlab=colnames(pcts2), yaxlab=rownames(pcts2), xcex=1.2, ycex=1.2,
			xgrid=T, ygrid=T, xlim=c(0.5,NV+0.5), ylim=c(-1.7,NJ+0.5), font.lab=2,
			xlas=3, main=mainText, mar=c( 6,6,4,2), ...)

	vIdent <- as.vector( ident1)
	vIdent[ vIdent < 90] <- 90
	colorScale1 <- color.scale( vIdent, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1), xrange=c(90,100))
	my_size_n_color( x=xlocs, y=ylocs, size=pcts1/100, border='black', col=NA,
			xcex=1.2, ycex=1.2, xgrid=T, ygrid=T, add=TRUE, lwd=2, ...)

	lines( c( -1,NV+2), c(0,0), lty=1, lwd=2, col=1)
	colorscale <- color.scale( 90:100, cs1=c(1,0.01), cs2=c(0,0), cs3=c(0,1), xrange=c(90,100))
	colorXlo <- round( NV * 0.33)
	colorXhi <- NV - colorXlo + 1
	color.legend( colorXlo, -1.6, colorXhi, -0.9, seq( -15, 15, by=3), rect.col=colorscale)
	text( (colorXlo+colorXhi)/2, -1.5, "GermLine Identity", pos=1, cex=1.05, font=2)

	sizeX <- round( NV * 0.10)
	draw.circle( sizeX, -1.3, sqrt(0.02), col="dodgerblue", border=NA)
	draw.circle( sizeX, -1.3, sqrt(0.05), border="black", col=NA, lwd=2)
	text( sizeX, -1.01, "Down: 5 -> 2%", pos=3, font=2, cex=1.1)
	sizeX <- round( NV * 0.90)
	draw.circle( sizeX, -1.3, sqrt(0.05), col="dodgerblue", border=NA)
	draw.circle( sizeX, -1.3, sqrt(0.02), border="black", col=NA, lwd=2)
	text( sizeX, -1.01, "Up: 2 -> 5%", pos=3, font=2, cex=1.1)
}

