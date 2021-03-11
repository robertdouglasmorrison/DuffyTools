# pipe.MarkerGenes.R -- evaluate transcriptomes based on a set of marker genes


`pipe.ScoreMarkerGenes` <- function( sampleIDset, markerDF, optionsFile="Options.txt", results.path=NULL, 
				folder=NULL, mode=c("absolute", "relative"), 
				groupOrder=NULL, col=NULL, main="", legend.cex=1, nFDRsimulations=100,
				forceYmax=NULL) {
	
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound="./results", verbose=F)
	}

	NS <- length( sampleIDset)
	grpFac <- factor( markerDF$Group)
	NG <- nlevels( grpFac)

	m <- pvm <- matrix( NA, nrow=NG, ncol=NS)
	rownames(m) <- rownames(pvm) <- levels(grpFac)
	colnames(m) <- colnames(pvm) <- sampleIDset
	colnames(m) <- sampleIDset
	colnames(pvm) <- paste( "FDR", sampleIDset, sep="_")

	prefix <- getCurrentSpeciesFilePrefix()

	# use either transcriptome or DE results
	if ( is.null(folder)) {
		myFolder <- "transcript"
		mySuffix <- paste( prefix, "Transcript.txt", sep=".")
		intensityColumn <- "RPKM_M"
	} else {
		myFolder <- file.path( "MetaResults", paste( prefix, folder, sep="."))
		mySuffix <- paste( prefix, "Meta.UP.txt", sep=".")
		intensityColumn <- "LOG2FOLD"
	}

	for ( i in 1:NS) {
		s <- sampleIDset[i]
		f <- file.path( results.path, myFolder, paste( s, mySuffix, sep="."))
		tbl <- read.delim( f, as.is=T)
		ans <- scoreMarkerGenes( tbl, markerDF, intensityColumn=intensityColumn, nFDRsimulations=nFDRsimulations)
		if ( is.null( ans)) next
		ord <- order( ans$Group)
		m[ , i] <- ans$Score[ ord]
		pvm[ , i] <- ans$FDR[ ord]
		cat( "\n", i, s, "\tBest: ", ans$Group[1], ans$Score[1])
	}
	cat( "\n")

	if ( ! is.null( groupOrder)) {
		if ( length(groupOrder) != NG) stop( "'groupOrder' length does not match group count")
		m <- m[ groupOrder, , drop=F]
		pvm <- pvm[ groupOrder, , drop=F]
	} else {
		groupOrder <- 1:NG
	}
	
	if ( is.null(col)) {
		groupColor <- rainbow( NG, end=0.76)
	} else {
		groupColor <- rep( col, length.out=NG)
	}

	# spacing for the legend depends on counts and label length...
	gap <- 2.5
	NSsizing <- max( NS, 3)
	bigX <- NSsizing * (NG+gap) * (1 + ((max( nchar(levels(grpFac)))+6)/120))
	ylabel <- "Marker Gene Score   (Absolute Scale)"
	mode <- match.arg( mode)
	if ( mode == "relative") {
		for ( i in 1:nrow(m)) {
			mymed <- mean( m[ i, ])
			m[ i, ] <- m[ i, ] - mymed
		}
		ylabel <- "Marker Gene Score   (Relative to Group Means)"
	}
	yLimits <- range(m,0) * 1.15
	if ( ! is.null(forceYmax)) {
		yLimits[2] <- as.numeric( forceYmax)
		if ( yLimits[1] < -1.0) yLimits[1] <- -(as.numeric(forceYmax))
	}

	barAns <- barplot( m, beside=T, col=groupColor, main=paste( "Marker Gene Profile:    ", main),
			las=if(NS<5) 1 else 3, xlim=c(1,bigX), ylim=yLimits, space=c(0,gap), 
			ylab=ylabel, cex.lab=1.1, cex.axis=1.1, font.lab=2, font.axis=2)
	lines( c(-10,bigX*2), c(0,0), lwd=1, lty=1, col=1)
	if ( nFDRsimulations > 0) {
		text( as.vector(barAns), as.vector(m), paste( "P=",as.vector(pvm),sep=""), pos=ifelse( as.vector(m) > 0, 3, 1), cex=0.75)
	}
	legend( "topright", levels(grpFac)[groupOrder], fill=groupColor, bg='white', cex=legend.cex) 

	return( data.frame( "Group"=rownames(m), m, pvm, row.names=1:NG, stringsAsFactors=F))
}


`pipe.ShowMarkerGenes` <- function( sampleID, markerDF, optionsFile="Options.txt", results.path=NULL, 
				folder=NULL, groupOrder=NULL, col=NULL, main="", legend.cex=1,
				names.cex=1.2, nFDRsimulations=100) {
	
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound="./results", verbose=F)
	}

	grpFac <- factor( markerDF$Group)
	NG <- nlevels( grpFac)

	prefix <- getCurrentSpeciesFilePrefix()

	# use either transcriptome or DE results
	if ( is.null(folder)) {
		myFolder <- "transcript"
		mySuffix <- paste( prefix, "Transcript.txt", sep=".")
		intensityColumn <- "RPKM_M"
		ylabel <- "Gene Rank Percentile in Transcriptome"
	} else {
		myFolder <- file.path( "MetaResults", paste( prefix, folder, sep="."))
		mySuffix <- paste( prefix, "Meta.UP.txt", sep=".")
		intensityColumn <- "LOG2FOLD"
		ylabel <- "Gene Rank Percentile in Diff Expression"
	}

	f <- file.path( results.path, myFolder, paste( sampleID, mySuffix, sep="."))
	tbl <- read.delim( f, as.is=T)
	genes <- tbl$GENE_ID
	isNG <- grep( "(ng)", genes, fixed=T)
	if ( length( isNG)) {
		tbl <- tbl[ -isNG, ]
		genes <- tbl$GENE_ID
	}
	genes <- shortGeneName( genes, keep=1)

	# offset for log scale...
	v <- tbl[[ intensityColumn]]
	vRank <- as.rankPercentile(v)
	Ngenes <- length( genes)

	if ( ! is.null( groupOrder)) {
		if ( length(groupOrder) != NG) stop( "'groupOrder' length does not match group count")
	} else {
		groupOrder <- 1:NG
	}
	
	if ( is.null(col)) {
		groupColor <- rainbow( NG, end=0.76)
	} else {
		groupColor <- rep( col, length.out=NG)
	}

	xlabel <- if ( NG > 4) NA else "Marker Gene Group Name"

	plot( 1, 1,  main=paste( "Marker Gene Profile:      SampleID =", sampleID, main), type="n", 
			xlim=c(0.4,NG*1.2), ylim=c( -6,105), ylab=ylabel, xlab=xlabel, xaxt="n",
			cex.lab=1.1, cex.axis=1.1, font.lab=2, font.axis=2)

	# calc and show where the marker genes land
	outX <- outY <- outPCH <- outCol <- vector()
	for ( j in 1:NG) {
		myGrp <- levels(grpFac)[ groupOrder[ j]]
		myColor <- groupColor[j]
		who <- which( markerDF$Group == myGrp)
		myGenes <- markerDF$GENE_ID[ who]
		myPCH <- ifelse( markerDF$Direction[ who] == "UP", 24, 6)
		where <- match( myGenes, genes, nomatch=0)
		use <- which( where > 0)
		outX <- c( outX, rep.int( j, length(use)))
		outY <- c( outY, vRank[where[use]])
		outPCH <- c( outPCH, myPCH[use])
		outCol <- c( outCol, rep.int(myColor,length(use)))
	}

	# extra annotations at the left edge
	xLeft <- min( outX) - diff( range(outX))/(NG*1.25)
	lines( c(-10, NG*2), c(50,50), lty=2, lwd=2, col='gray50')
	text( c(xLeft,xLeft), c(80,28), c( "Positive Weighting", "Negative Weighting"), srt=90, 
		col='gray50', cex=1.03, font=2)

	points( jitter(outX, factor=1.5), outY, pch=outPCH, col=outCol, bg=outCol)
	axis( 1, at=1:NG, levels(grpFac)[groupOrder], font=2, cex.axis=names.cex, las=if (NG>4) 3 else 1)

	legend( "topright", levels(grpFac)[groupOrder], fill=groupColor, bg='white', cex=legend.cex) 
	legend( "bottomright", c( "Positve Marker", "Negative Marker"), pch=c(24,6), 
		col="brown", pt.bg="brown", pt.cex=2, bg='white', cex=1.1) 

	ans <- scoreMarkerGenes( tbl, markerDF, intensityColumn=intensityColumn, nFDRsimulations=nFDRsimulations)
	ord <- match( ans$Group, levels(grpFac)[groupOrder])

	text( c( xLeft, ord), -3.5, c( "SCORE:", as.character(round( ans$Score, digits=1))), font=2, cex=0.95,
		srt=if(NG>5) 90 else 0)
	if ( nFDRsimulations > 0 && NG <= 5) {
		text( c( xLeft, ord), -6, c( "P-value:", as.character(round( ans$FDR, digits=3))), font=2, cex=0.75)
	}

	return( ans)
}


`scoreMarkerGenes` <- function( transcript, markerDF, geneColumn="GENE_ID", intensityColumn="RPKM_M",
				nFDRsimulations=100) {

	genes <- transcript[[ geneColumn]]
	if ( is.null( genes)) {
		cat( "\nGeneID column not found: ", geneColumn)
		return( NULL)
	}
	
	isNG <- grep( "(ng)", genes, fixed=T)
	if ( length( isNG)) {
		transcript <- transcript[ -isNG, ]
		genes <- transcript[[ geneColumn]]
	}

	# we will only want the gene symbol...
	genes <- shortGeneName( genes, keep=1)

	inten <- transcript[[ intensityColumn]]
	if ( is.null( inten)) {
		cat( "\nIntensity column not found: ", intensityColumn)
		return( NULL)
	}

	# make sure the transcript is in order
	ord <- order( inten, decreasing=T)
	if ( any( ord != 1:nrow(transcript))) {
		cat( "\nNot in Intensity order..  Resorting..")
		genes <- genes[ ord]
		inten <- inten[ ord]
	}
	NG <- length(genes)
	generank <- as.rankPercentile( NG:1)
	
	neededColumns <- c( "GENE_ID", "Group", "Direction")
	if ( ! all( neededColumns %in% colnames(markerDF))) {
		cat( "Expected 'MarkGene' data frame columns not found.   Need: ", neededColumns)
		return(NULL)
	}

	# make sure most of the markers are seen
	markgenes <- markerDF$GENE_ID
	NMG <- nrow(markerDF)
	where <- match(markgenes, genes, nomatch=0)
	if ( sum( where != 0) < NMG*0.5) {
		cat( "\nDebug:  N_Missing: ", sum( where == 0), " \tN_Found: ", sum( where > 0),"\n")
		cat( "Too many 'Marker Genes' not found in transcriptome..  Check current species..")
		return( NULL)
	}
	markerDF$SCORE <- ifelse( toupper( markerDF$Direction) == "UP", 1, -1)
	markerDF$SCORE[ where == 0] <- 0
	markerFac <- factor( markerDF$Group)
	NGRPS <- nlevels( markerFac)

	RANK_MIDPT <- 50

	scores <- tapply( 1:NMG, markerFac, function(x) {

				mywhere <- where[x]
				myscores <- markerDF$SCORE[x]
				myrankPcts <- rep.int( RANK_MIDPT, length(x))
				myrankPcts[ mywhere > 0] <- generank[ mywhere]
				
				# turn these 1..N into centered about zero and scaled to -1,1
				myPtiles <- (myrankPcts - RANK_MIDPT) / RANK_MIDPT
	
				# multiply by their sign and scale to number of markers
				score <- sum( myPtiles * myscores) * 100 / sum( mywhere > 0)
				return( score)
			})

	# if we do FDR, permute and redo the scoring
	fdr <- rep.int( NA, NGRPS)
	if ( nFDRsimulations > 0) {
		randScores <- matrix( 0, nrow=NGRPS, ncol=nFDRsimulations)
		for ( j in 1:nFDRsimulations) {
			randGenes <- sample( genes)
			where <- match(markgenes, randGenes, nomatch=0)
			randScores[,j] <- tapply( 1:NMG, markerFac, function(x) {
				mywhere <- where[x]
				myscores <- markerDF$SCORE[x]
				myrankPcts <- rep.int( RANK_MIDPT, length(x))
				myrankPcts[ mywhere > 0] <- generank[ mywhere]
				myPtiles <- (myrankPcts - RANK_MIDPT) / RANK_MIDPT
				score <- sum( myPtiles * myscores) * 100 / sum( mywhere > 0)
				return( score)
			})
		}
		# now count up how often did we do better than real score
		for ( i in 1:NGRPS) {
			if ( scores[i] >= 0) {
				fdr[i] <- sum( randScores[ i, ] >= scores[i]) / nFDRsimulations
			} else {
				fdr[i] <- sum( randScores[ i, ] <= scores[i]) / nFDRsimulations
			}
		}
	}
	
	out <- data.frame( "Group"=levels(markerFac), "Score"=round(scores,digits=2), 
			"FDR"=round( fdr, digit=4), stringsAsFactors=F)
	ord <- order( scores, decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:NGRPS
	return( out)
}
