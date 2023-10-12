# markerGenes.R -- tools to evaluate transcriptomes based on a set of marker genes


`showMarkerGenes` <- function( x, markerDF, sampleID="", geneColumn="GENE_ID", intensityColumn="RPKM_M",
				mode=c("transcriptome","proteome"), groupOrder=NULL, col=NULL, label="", 
				pt.cex=1.0, legend.cex=1.0, names.cex=1.2, nFDRsimulations=100, ...) {
	
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

	xlabel <- if ( NG > 5) NA else "Marker Gene Group Name"
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
	axis( 1, at=1:NG, levels(grpFac)[groupOrder], font=2, cex.axis=names.cex, las=if (NG>5) 3 else 1)

	legend( "topright", levels(grpFac)[groupOrder], fill=groupColor, bg='white', cex=legend.cex) 
	legend( "bottomright", c( "Positve Marker", "Negative Marker"), pch=c(24,6), 
		col="brown", pt.bg="brown", pt.cex=2, bg='white', cex=legend.cex) 

	ans <- calcMarkerGenesScore( tbl, markerDF, geneColumn=geneColumn, 
				intensityColumn=intensityColumn, mode=mode, 
				nFDRsimulations=nFDRsimulations)
	if ( is.null(ans)) return(NULL)

	ord <- match( ans$Group, levels(grpFac)[groupOrder])

	text( c( 0.6, ord), -9, c( "SCORE:", as.character(round( ans$Score, digits=2))), font=2, cex=1.15,
		srt=if(NG>5) 90 else 0)
	if ( nFDRsimulations > 0 && NG <= 5) {
		text( c( 0.6, ord), -15, c( "P-value:", as.character(round( ans$FDR, digits=3))), font=2, cex=0.97)
	}
	dev.flush()

	return( ans)
}


`calcMarkerGenesScore` <- function( x, markerDF, geneColumn="GENE_ID", intensityColumn="RPKM_M",
					mode=c("transcriptome","proteome"), sort=TRUE, nFDRsimulations=100) {

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

	# make sure enough of the markers are seen
	markgenes <- markerDF$GENE_ID
	NMG <- nrow(markerDF)
	NGG <- nrow(tbl)
	minGenesToFind <- min( NMG, NGG) * 0.25

	if ( mode == "transcriptome") {
		noMatchValue <- 0
	} else {
		noMatchValue <- Ngenes + 1
	}
	where <- match(markgenes, genes, nomatch=noMatchValue)
	if ( sum( where != 0) < minGenesToFind) {
		cat( "\nToo many 'Marker Genes' not found in transcriptome..  Check current species..")
		cat( "\nN_Genes in data set: ", Ngenes)
		cat( "\nN_MarkerGenes:       ", NMG)
		cat( "\nN_MarkerGenes found: ", sum( where != 0))
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
			where <- match(markgenes, randGenes, nomatch=noMatchValue)
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
	if (sort) {
		ord <- order( scores, decreasing=T)
		out <- out[ ord, ]
		rownames(out) <- 1:NGRPS
	}

	return( out)
}


`compareMarkerGenesScores` <- function( fileSet, sampleIDset, markerDF, 
				mode=c("transcriptome","proteome"), displayMode=c("absolute", "relative"), 
				groupOrder=NULL, geneColumn="GENE_ID", intensityColumn="RPKM_M",
				col=NULL, label="", plotMode=c("lines", "bars", "none"), 
				pch=22, pt.cex=2.5, lwd=4, label.cex=1, legend.cex=1.1, nFDRsimulations=100, ...) {
	
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

	m <- pv <- matrix( NA, nrow=NG, ncol=NS)
	rownames(m) <- rownames(pv) <- levels(grpFac)
	colnames(m) <- colnames(pv) <- sampleIDset

	for ( i in 1:NS) {
		s <- sampleIDset[i]
		f <- fileSet[i]
		if ( is.na(s) || is.na(f)) next
		tbl <- read.delim( f, as.is=T)
		ans <- calcMarkerGenesScore( tbl, markerDF, geneColumn=geneColumn, 
					intensityColumn=intensityColumn, mode=mode,
					nFDRsimulations=nFDRsimulations)
		if ( is.null( ans)) next
		ord <- order( ans$Group)
		m[ , i] <- ans$Score[ ord]
		pv[ , i] <- ans$FDR[ ord]
		cat( "\r", i, s, "\tBest: ", ans$Group[1], ans$Score[1], ans$FDR[1])
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
		return( list( "score"=m, "p.value"=pv))
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
		ylim <- range( m, na.rm=T) * 1.05
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
		# show P-values
		pvTxt <- rep.int( "", NS)
		bestP <- apply( pv, 2, min, na.rm=T)
		bestS <- apply( m, 2, max, na.rm=T)
		pvTxt[ bestP < 0.05] <- "*"
		pvTxt[ bestP < 0.01] <- "**"
		pvTxt[ bestP < 0.001] <- "***"
		text( 1:NS, bestS, pvTxt, pos=3)

		legend( "topright", rownames(m), fill=groupColor, bg='white', cex=legend.cex)
	}

	if ( length( isNA)) {
		m <- m[ , -isNA]
		pv <- pv[ , -isNA]
	}
	return( list( "score"=m, "p.value"=pv))
}


`markerGeneOverlap` <- function( markerDF, min.overlap=2) {

	# try to summarize how commonly genes are UP markers in 2+ groups
	neededColumns <- c( "GENE_ID", "Group", "Direction")
	if ( ! all( neededColumns %in% colnames(markerDF))) {
		cat( "Expected 'MarkerGene' data frame columns not found.   Need: ", neededColumns)
		return(NULL)
	}

	# get the IDs we see 2+ times
	markerDF <- subset( markerDF, Direction == "UP")
	upGenes <- markerDF$GENE_ID
	dupGenes <- sort( unique( upGenes[ duplicated( upGenes)]))

	# see what sets of 2+ cell types each gene was in
	dupGroups <- sapply( dupGenes, function(x) {
			hits <- which( markerDF$GENE_ID == x)
			if (length(hits) < 2) return("")
			myGrps <- sort( unique( markerDF$Group[hits]))
			if ( length( myGrps) > 4) cat( "\nWarn: Very common gene: ", length(myGrps), "|", x, "|", myGrps)
			return( paste( myGrps, collapse=" = "))
		})
	if ( all( is.na(dupGroups))) return( data.frame())

	# use the aveage number of genes per group to assess how common a duplication event is
	nGperG <- round( mean( tapply( 1:nrow(markerDF), factor(markerDF$Group), length)))

	# make the counts are how often each 2+ cell type call was seen
	dupTbl <- table( dupGroups)
	dupPcts <- round( as.numeric( dupTbl) * 100 / nGperG, digits=1)

	# lastly, make that list of genes that was in each group of cell types
	dupGrpGeneList <- sapply( names(dupTbl), function(x) { wh <- which( dupGroups == x); 
					return( paste( sort( unique( dupGenes[wh])), collapse="; "))
				})

	out <- data.frame( "Groups.Sharing.a.Gene"=names(dupTbl), "Gene.Count"=as.numeric(dupTbl), 
			"Gene.Percentage"=dupPcts, "Gene.List"=dupGrpGeneList, stringsAsFactors=F)
	ord <- order( dupPcts, decreasing=T)
	out <- out[ ord, ]
	drops <- which( out$Gene.Count < min.overlap)
	if ( length(drops)) out <- out[ -drops, ]
	rownames(out) <- 1:nrow(out)
	out
}


`createMarkerGenes` <- function( folder, file.pattern=paste( getCurrentSpeciesFilePrefix(), "Meta.JOINED.txt$", sep="."),
				Ngenes=100, geneColumn="GENE_ID", sep="\t", UP.only=FALSE) {

	# given a folder that should hold DE results, find all the groups' results
	fset <- dir( folder, pattern=file.pattern, recursive=FALSE, full.name=T)
	if (! length( fset)) {
		cat( "\nWarning:  No DE files found that match the given file pattern: ", file.pattern)
		return( NULL)
	}

	out <- data.frame()
	for (f in fset) {
		tbl <- read.delim( f, as.is=TRUE, sep=sep)
		if ( ! geneColumn %in% colnames(tbl)) {
			cat( "\nWarning:  No column name matches the given gene field: ", geneColumn, " in file: ", f)
			next
		}
		grpName <- sub( paste( ".", file.pattern, sep=""), "", basename(f))
		Nuse <- min( Ngenes, nrow(tbl)/2)
		smlTbl <- tbl[ 1:Nuse, ]
		sml <- data.frame( "GENE_ID"=shortGeneName( smlTbl[[ geneColumn]], keep=1), "Group"=grpName, 
				"Direction"="UP", "Rank"=1:nrow(smlTbl), stringsAsFactors=F)
		out <- rbind( out, sml)
		if ( UP.only) next
		tbl <- tbl[ rev( 1:nrow(tbl)), ]
		smlTbl <- tbl[ 1:Nuse, ]
		sml <- data.frame( "GENE_ID"=shortGeneName( smlTbl[[ geneColumn]], keep=1), "Group"=grpName, 
				"Direction"="DOWN", "Rank"=1:nrow(smlTbl), stringsAsFactors=F)
		out <- rbind( out, sml)
		cat( "  ", grpName, nrow(out), "\n")
	}
	return( out)
}

