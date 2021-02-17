# markerGenes.R -- tools to evaluate transcriptomes based on a set of marker genes


`showMarkerGenes` <- function( x, markerDF, sampleID="", geneColumn="GENE_ID", intensityColumn="RPKM_M",
				mode=c("transcriptome","proteome"), groupOrder=NULL, col=NULL, label="", pt.cex=1.0, ...) {
	
	grpFac <- factor( markerDF$Group)
	NG <- nlevels( grpFac)
	if ( ! is.null( groupOrder)) {
		if ( length(groupOrder) != NG) {
			cat( "\n'groupOrder' length does not match marker group count: ", NG)
			return(NULL)
		}
	} else {
		groupOrder <- 1:NG
	}

	mode <- match.arg( mode)

	# arg 1 can be a data.frame or a file name
	if ( is.data.frame(x)) {
		tbl <- x
	} else {
		f <- x
		if ( length(f) > 1) warning( "'showMarkerGenes': Using only the first given filename")
		if ( ! file.exists( f[1])) {
			cat( "\n", mode, " file not found:  ", f[1])
			return(NULL)
		}
		tbl <- read.delim( f[1], as.is=T)
		if ( sampleID == "") sampleID <- basename(f)
	}

	if ( ! all( c( geneColumn, intensityColumn) %in% colnames(tbl))) {
		cat( "\nSome needed columns not seen:  ", geneColumn, intensityColumn)
		cat( "\nFound:  ", colnames(tbl))
		return( NULL)
	}
	genes <- tbl[[ geneColumn]]
	isNG <- grep( "(ng)", genes, fixed=T)
	if ( length( isNG)) {
		tbl <- tbl[ -isNG, ]
		genes <- genes[ -isNG]
	}
	genes <- shortGeneName( genes, keep=1)
	markerDF$GENE_ID <- shortGeneName( markerDF$GENE_ID, keep=1)

	v <- tbl[[ intensityColumn]]
	vRank <- as.rankPercentile(v)
	Ngenes <- length( genes)

	if ( is.null(col)) {
		groupColor <- rainbow( NG, end=0.76)
	} else {
		groupColor <- rep( col, length.out=NG)
	}

	xlabel <- if ( NG > 4) NA else "Marker Gene Group Name"
	ylabel <- "Gene Rank Percentile"

	if ( mode == "proteome") {
		ylabel <- "Proteome Rank Percentile (50% = detected)"
		vRank <- (vRank/2) + 50
		vRank <- c( vRank, 0)
	}

	plot( 1, 1,  main=paste( "Marker Gene Profile:      SampleID =", sampleID, label), type="n", 
			xlim=c(0.4,NG*1.2+0.3), ylim=c( -15,103), ylab=ylabel, xlab=xlabel, xaxt="n", yaxt="n",
			cex.lab=1.1, cex.axis=1.1, font.lab=2, font.axis=2)
	axis( 2, at=seq(0,100,20), font=2, cex.axis=1.2)

	# show the marker genes
	outX <- outY <- outPCH <- outCol <- vector()
	for ( j in 1:NG) {
		myGrp <- levels(grpFac)[ groupOrder[ j]]
		myColor <- groupColor[j]
		who <- which( markerDF$Group == myGrp)
		myGenes <- markerDF$GENE_ID[ who]
		myPCH <- ifelse( markerDF$Direction[ who] == "UP", 24, 6)
		# in a transcriptome, every gene should be there
		# in a proteome, only 'detected' genes should be seen
		if ( mode == "transcriptome") {
			where <- match( myGenes, genes, nomatch=0)
			use <- which( where > 0)
		} else {
			where <- match( myGenes, genes, nomatch=Ngenes+1)
			use <- which( where > 0)
		}
		outX <- c( outX, rep.int( j, length(use)))
		outY <- c( outY, vRank[where[use]])
		outPCH <- c( outPCH, myPCH[use])
		outCol <- c( outCol, rep.int(myColor,length(use)))
	}
	lines( c(-10, NG*2), c(50,50), lty=2, lwd=2, col='gray50')
	text( c(0.4,0.4), c(75,25), c( "Positive Weights", "Negative Weights"), srt=90, 
		col='gray50', cex=1.03, font=2)

	points( jitter(outX, factor=1.5), outY, pch=outPCH, col=outCol, bg=outCol, cex=pt.cex)
	axis( 1, at=1:NG, levels(grpFac)[groupOrder], font=2, cex.axis=1.2, las=if (NG>4) 3 else 1)

	legend( "topright", levels(grpFac)[groupOrder], fill=groupColor, bg='white', cex=1.15) 
	legend( "bottomright", c( "Positve Marker", "Negative Marker"), pch=c(24,6), 
		col="brown", pt.bg="brown", pt.cex=2, bg='white', cex=1.1) 

	ans <- calcMarkerGenesScore( tbl, markerDF, geneColumn=geneColumn, 
				intensityColumn=intensityColumn, mode=mode)
	ord <- match( names( ans), levels(grpFac)[groupOrder])

	text( c( 0.6, ord), -9, c( "SCORE:", as.character(round( ans, digits=2))), font=2, cex=1.15,
		srt=if(NG>4) 90 else 0)

	return( ans)
}


`calcMarkerGenesScore` <- function( x, markerDF, geneColumn="GENE_ID", intensityColumn="RPKM_M",
					mode=c("transcriptome","proteome"), sort=TRUE) {

	mode <- match.arg( mode)

	# arg 1 can be a data.frame or a file name
	if ( is.data.frame(x)) {
		tbl <- x
	} else {
		f <- x
		if ( length(f) > 1) warning( "'showMarkerGenes': Using only the first given filename")
		if ( ! file.exists( f[1])) {
			cat( "\n", mode, " file not found:  ", f[1])
			return(NULL)
		}
		tbl <- read.delim( f[1], as.is=T)
	}

	if ( ! all( c( geneColumn, intensityColumn) %in% colnames(tbl))) {
		cat( "\nSome needed columns not seen:  ", geneColumn, intensityColumn)
		cat( "\nFound:  ", colnames(tbl))
		return( NULL)
	}
	genes <- tbl[[ geneColumn]]
	isNG <- grep( "(ng)", genes, fixed=T)
	if ( length( isNG)) {
		tbl <- tbl[ -isNG, ]
		genes <- genes[ -isNG]
	}
	genes <- shortGeneName( genes, keep=1)
	markerDF$GENE_ID <- shortGeneName( markerDF$GENE_ID, keep=1)

	v <- tbl[[ intensityColumn]]
	generank <- as.rankPercentile( v)
	Ngenes <- length( genes)
	
	if ( mode == "proteome") {
		generank <- (generank/2) + 50
		generank <- c( generank, 0)
	}

	neededColumns <- c( "GENE_ID", "Group", "Direction")
	if ( ! all( neededColumns %in% colnames(markerDF))) {
		cat( "Expected 'MarkerGene' data frame columns not found.   Need: ", neededColumns)
		return(NULL)
	}

	# make sure most of the markers are seen
	markgenes <- markerDF$GENE_ID
	NMG <- nrow(markerDF)

	if ( mode == "transcriptome") {
		where <- match(markgenes, genes, nomatch=0)
	} else {
		where <- match(markgenes, genes, nomatch=Ngenes+1)
	}
	if ( sum( where == 0) > NMG*0.5) {
		cat( "Too many 'Marker Genes' not in transcriptome..  Check current species..")
		return( NULL)
	}
	markerDF$SCORE <- ifelse( toupper( markerDF$Direction) == "UP", 1, -1)
	markerDF$SCORE[ where == 0] <- 0
	markerFac <- factor( markerDF$Group)

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

	out <- scores
	if (sort) out <- sort( scores, decreasing=T)

	return( out)
}


`compareMarkerGenesScores` <- function( fileSet, sampleIDset, markerDF, 
				mode=c("transcriptome","proteome"), displayMode=c("absolute", "relative"), 
				groupOrder=NULL, geneColumn="GENE_ID", intensityColumn="RPKM_M",
				col=NULL, label="", plotMode=c("lines", "bars", "none"), 
				pch=22, pt.cex=2.5, lwd=4, label.cex=1, legend.cex=1.1, ...) {
	
	mode <- match.arg( mode)

	found <- ( file.exists(fileSet) | is.na(fileSet))
	if ( any( ! found)) {
		cat( "\nSome files not found: ", fileSet[!found])
		fileSet <- fileSet[ found]
		sampleIDset <- sampleIDset[ found]
	}
	NS <- length( sampleIDset)
	isNA <- which( is.na( fileSet))
	grpFac <- factor( markerDF$Group)
	NG <- nlevels( grpFac)
	if ( any( c( NS, NG) < 1)) stop( "Not enough samples or marker groups")

	m <- matrix( NA, nrow=NG, ncol=NS)
	rownames(m) <- levels(grpFac)
	colnames(m) <- sampleIDset

	for ( i in 1:NS) {
		s <- sampleIDset[i]
		f <- fileSet[i]
		if ( is.na(s) || is.na(f)) next
		tbl <- read.delim( f, as.is=T)
		ans <- calcMarkerGenesScore( tbl, markerDF, geneColumn=geneColumn, 
					intensityColumn=intensityColumn, mode=mode)
		if ( is.null( ans)) next
		ord <- order( names(ans))
		m[ , i] <- ans[ ord]
		cat( "\r", i, s, "\tBest: ", names(ans)[1], ans[1])
	}
	cat( "\n")

	if ( ! is.null( groupOrder)) {
		if ( length(groupOrder) != NG) stop( "'groupOrder' length does not match group count")
		m <- m[ groupOrder, ]
	} else {
		groupOrder <- 1:NG
	}
	
	if ( is.null(col)) {
		groupColor <- rainbow( NG, end=0.82)
	} else {
		groupColor <- rep( col, length.out=NG)
	}

	# keep the raw scores, or relative
	ylabel <- "Marker Gene Score   (Absolute Scale)"
	displayMode <- match.arg( displayMode)
	if ( displayMode == "relative") {
		for ( i in 1:nrow(m)) {
			mymed <- median( m[ i, ], na.rm=T)
			m[ i, ] <- m[ i, ] - mymed
		}
		ylabel <- "Marker Gene Score   (Relative to Group Medians)"
	}

	# spacing for the legend depends on counts and label length...
	plotMode <- match.arg( plotMode)
	if ( plotMode == "none") {
		return( m)
	}
	if ( plotMode == "bars") {
		gap <- 2.5
		bigX <- NS * (NG+gap) * (1 + ((max( nchar(levels(grpFac)))+6)/120))
		barplot( m, beside=T, col=groupColor, main=paste( "Marker Gene Comparison:    ", label),
			las=if(NS<5) 1 else 3, xlim=c(1,bigX), space=c(0,gap), ylab=ylabel, 
			cex.lab=1.1, cex.axis=1.1, font.lab=2, font.axis=2)
		legend( "topright", levels(grpFac)[groupOrder], fill=groupColor, bg='white', cex=legend.cex) 
	}
	if ( plotMode == "lines") {
		xLegendPad <- 1.3
		if ( NS < 6) xLegendPad <- 1.5
		ylim <- range( m, na.rm=T)
		plot( 1, 1, type="n", main=paste( "Marker Genes Scoring:   ", label), ylab=ylabel,
			xlim=c(0.5-(NS/10),NS*xLegendPad), ylim=ylim, font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1, 
			xaxt="n", xlab=NA, ...)
		axis( 1, at=1:NS, colnames(m), las=3, font=2, cex=1.05)
		lines( c(-10, NS*3), c(0,0), lty=2, lwd=2, col='gray50')

		for ( i in 1:NS) {
			v <- m[ ,i]
			ord <- order( v)
			points( rep.int(i,length(ord)), v[ord], bg=groupColor[ord], pch=pch, cex=pt.cex)
		}
		for( j in 1:NG) {
			lines( 1:NS, m[ j, ], col=groupColor[j], lwd=lwd)
			text( 0.9, m[j,1], rownames(m)[j], pos=2, cex=label.cex)
			text( NS+0.1, m[j,NS], rownames(m)[j], pos=4, cex=label.cex)
		}

		legend( "topright", rownames(m), fill=groupColor, bg='white', cex=legend.cex)
	}

	if ( length( isNA)) {
		m <- m[ , -isNA]
	}
	return( m)
}


`markerGeneOverlap` <- function( markerDF) {

	# try to summarize how commonly genes are UP markers in 2+ groups
	neededColumns <- c( "GENE_ID", "Group", "Direction")
	if ( ! all( neededColumns %in% colnames(markerDF))) {
		cat( "Expected 'MarkerGene' data frame columns not found.   Need: ", neededColumns)
		return(NULL)
	}

	markerDF <- subset( markerDF, Direction == "UP")
	upGenes <- markerDF$GENE_ID
	dupGenes <- upGenes[ duplicated( upGenes)]

	dupGroups <- unlist( sapply( dupGenes, function(x) {
			hits <- which( markerDF$GENE_ID == x)
			if (length(hits) < 2) return( NA)
			myGrps <- sort( unique( markerDF$Group[hits]))
			if ( length( myGrps) > 4) cat( "\nWarn: Very common gene: ", length(myGrps), "|", x, "|", myGrps)
			return( paste( myGrps, collapse=" = "))
		}))
	if ( all( is.na(dupGroups))) return( data.frame())

	# use the aveage number of genes per group to assess how common a duplication event is
	nGperG <- round( mean( tapply( 1:nrow(markerDF), factor(markerDF$Group), length)))

	dupTbl <- table( dupGroups)
	dupPcts <- round( as.numeric( dupTbl) * 100 / nGperG, digits=1)

	out <- data.frame( "Groups.Sharing.a.Gene"=names(dupTbl), "Gene.Count"=as.numeric(dupTbl), "Gene.Percentage"=dupPcts, 
			stringsAsFactors=F)
	ord <- order( dupPcts, decreasing=T)
	out <- out[ ord, ]
	drops <- which( out$Gene.Count < 2)
	if ( length(drops)) out <- out[ -drops, ]
	rownames(out) <- 1:nrow(out)
	
	out
}


