# transcriptDeconvolution.R -- higher level code to do transcript deconvolution 
#				uses lower level 'transcriptBlend' tool


`getTranscriptDeconvolutionTargetMatrix` <- function( speciesID=getCurrentSpecies(), target=NULL) {

	if ( speciesID != getCurrentSpecies()) {
		oldSpecies <- getCurrentSpecies()
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( oldSpecies))
	}

	# let's update this system to use the exact same data as from the cell type tools
	CellTypeSetup()
	targetM <- getCellTypeMatrix()
	return( targetM)
	
	# old method...
	myPrefix <- getCurrentSpeciesFilePrefix()
	CellTypeSetup()
	reference <- getCellTypeReference()
	datafile <- paste( myPrefix, reference, "TargetMatrix", sep=".")
	targetM <- NULL
	data( list=datafile, package="DuffyTools", envir=environment())
	if ( is.null( targetM)) cat( "\nFailed to load deconvolution target matrix.  Tried: ", datafile)
	return( targetM)
}


`buildTranscriptDeconvolutionTargetMatrix` <- function( dimensionIDset, transcript.files, de.files, nGenesPerDimension=NULL,
					speciesID=getCurrentSpecies(), dropGenePattern=NULL, dropProductPattern=NULL) {

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	geneMap <- subset( getCurrentGeneMap(), REAL_G=TRUE)

	# step 1:  make sure all files and sampleIDs are in agreement
	NS <- length(dimensionIDset)
	NTF <- length(transcript.files)
	NDEF <- length(de.files)
	if ( ! all( c( NTF, NDEF) == NS)) stop( "All arguments must be same length: dimensionIDset, transcript.files, de.files")
	for ( i in 1:NS) {
		sid <- dimensionIDset[i]
		tfid <- basename( transcript.files[i])
		defid <- basename( de.files[i])
		if ( ! grepl( sid, tfid, fixed=T)) stop( paste( "SampleID not in Transcript filename: ", sid, tfid))
		if ( ! grepl( sid, defid, fixed=T)) stop( paste( "SampleID not in DE ratios filename: ", sid, defid))
	}
	if ( ! all( grep( prefix, basename(transcript.files)) %in% 1:NS)) stop( "Expected species prefix not in all transcript files")
	if ( ! all( grep( prefix, basename(de.files)) %in% 1:NS)) stop( "Expected species prefix not in all DE ratio files")
	if ( ! all( file.exists( transcript.files))) stop( "Some transcript files not found")
	if ( ! all( file.exists( de.files))) stop( "Some DE ratio files not found")
	# OK, the files and IDs are self consistent

	# step 2:  get the reference transcriptomes for these samples
	cat( "\nLoading reference transcriptomes..\n")
	targetM <- expressionFileSetToMatrix( transcript.files, dimensionIDset, verbose=T)
	# force check that the genes we see are in this annotation
	drops <- which( ! ( rownames(targetM) %in% geneMap$GENE_ID))
	if ( length(drops)) {
		cat( "\nSome GeneIDs not in current species genome annotation.  Dropping: ", length( drops))
		targetM <- targetM[ -drops, ]
	}

	# step 3:  gather the genes that will be the marker gene set for each dimension
	if ( is.null( nGenesPerDimension)) {
		gpd <- round( nrow(geneMap) / NS / 2)
		nGenesPerDimension <- gpd
	}
	cat( "\nLoading marker genes..\n")
	markerGenes <- markerRanks <- markerTypes <- vector()
	for ( i in 1:NS) {
		f <- de.files[i]
		tmp <- read.delim( f, as.is=T)
		drops1 <- drops2 <- vector()
		if ( ! is.null( dropGenePattern)) {
			drops1 <- grep( dropGenePattern, tmp$GENE_ID)
		}
		if ( ! is.null( dropProductPattern)) {
			drops2 <- grep( dropProductPattern, tmp$PRODUCT)
		}
		drops <- sort( unique( c( drops1, drops2)))
		if ( length(drops)) {
			cat( "\nDropping unwanted genes: ", length(drops), "\n")
			tmp <- tmp[ -drops, ]
		}
		# grap up to 2N here, then trim out excess later
		Ng2 <- nGenesPerDimension * 2
		if ( Ng2 > nrow(tmp)) Ng2 <- nrow(tmp)
		markerGenes <- c( markerGenes, tmp$GENE_ID[1:Ng2])
		markerRanks <- c( markerRanks, 1:Ng2)
		markerTypes <- c( markerTypes, rep( dimensionIDset[i], Ng2))
		cat( "\r", i, dimensionIDset[i], Ng2, length(markerGenes))
	}

	# step 4:  keep just one copy of each gene, and track its type and rank
	ord <- order( markerGenes, markerRanks, markerTypes)
	best <- match( sort( unique( markerGenes)), markerGenes[ord])
	bestMarkerGenes <- markerGenes[ ord[ best]]
	bestMarkerRanks <- markerRanks[ ord[ best]]
	bestMarkerTypes <- markerTypes[ ord[ best]]

	# now enforce the N limit
	dropExcess <- vector()
	for (ct in sort( unique( bestMarkerTypes))) {
		now <- which( bestMarkerTypes == ct)
		if ( length(now) > nGenesPerDimension) {
			dropNow <- setdiff( now, now[1:nGenesPerDimension])
			dropExcess <- c( dropExcess, dropNow)
		}
	}
	if ( length( dropExcess)) {
		bestMarkerGenes <- bestMarkerGenes[ -dropExcess]
		bestMarkerRanks <- bestMarkerRanks[ -dropExcess]
		bestMarkerTypes <- bestMarkerTypes[ -dropExcess]
	}

	# now finalize with just those best
	# given this set as the proposed marker genes, make the target matrix we will use
	who <- match( bestMarkerGenes, rownames(targetM))
	targetM <- targetM[ who, ]
	out <- data.frame( "GENE_ID"=bestMarkerGenes, "Dimension"=bestMarkerTypes, 
				"Rank"=bestMarkerRanks, "PRODUCT"=gene2Product( bestMarkerGenes),
				targetM, stringsAsFactors=F)
	outfile <- paste( prefix, "Deconvolution.TargetMatrix.txt", sep=".")
	write.table( out, outfile, sep="\t", quote=F, row.names=F)
	outfile2 <- paste( prefix, "Deconvolution.TargetMatrix.rda", sep=".")
	save( targetM, file=outfile2)
	cat( "\nWrote new Deconvolution Target Matrix: ", nGenesPerDimension, "  \tN_Row: ", nrow(targetM))
	return()
}


`fileSet.TranscriptDeconvolution` <- function( files, fids, targetM=getTranscriptDeconvolutionTargetMatrix(), 
					geneColumn="GENE_ID", intensityColumn="RPKM_M",
					sep="\t", useLog=FALSE, minIntensity=0, 
					arrayFloorIntensity=NULL, dropLowVarianceGenes=NULL, geneUniverse=NULL,
					algorithm=c("port","default","plinear","LM","GenSA","steepDescent"),
					startFractions=NULL, plot=TRUE, plot.path=".", plot.col=NULL, 
					label="", dropRiboClearGenes=TRUE, verbose=TRUE) {

	if ( length(files) != length(fids)) {
		cat( "\nLength mismatch:  'files' and 'fids' must be same length")
		return(NULL)
	}

	# evaluate all files, and then decide how to visualize
	found <- file.exists(files)
	if ( sum(found) < length(files)) {
		cat( "\nSome transcript files not found: ", fids[ !found])
		files <- files[ found]
		fids <- fids[ found]
	}
	NS <- length( files)
	ND <- ncol(targetM)

	ans <- matrix( 0, nrow=ND, ncol=NS)
	colnames(ans) <- fids
	rownames(ans) <- colnames(targetM)

	algorithm <- match.arg( algorithm)
	if (algorithm == "GenSA") require( GenSA)
	if (algorithm == "LM") require( minpack.lm)

	# try to do in multi core
	myStarts <- NULL
	mcAns <- multicore.lapply( 1:NS, FUN=function(i) {
			f <- files[ i]
			tbl <- read.delim( f, sep=sep, as.is=T)
			if (i > 1) verbose <<- FALSE
			# allow usuing specific starting percentages
			# note that previous runs output 0 to 100, but the low level tool wants 0 to 1.0
			if ( ! is.null( startFractions)) {
				wh <- which( colnames(startFractions) == fids[i])[1]
				if ( ! is.na(wh)) {
					#cat( "  Given start %s..")
					myStarts <- as.numeric(startFractions[ ,wh]) / 100
				}
			}
			fitAns <- fit.transcriptBlend( tbl, targetM, geneColumn=geneColumn, 
					intensityColumn=intensityColumn, useLog=useLog, 
					minIntensity=minIntensity, 
					arrayFloorIntensity=arrayFloorIntensity, 
					dropLowVarianceGenes=dropLowVarianceGenes, 
					geneUniverse=geneUniverse, dropRiboClearGenes=dropRiboClearGenes, 
					algorithm=algorithm, startFractions=myStarts, verbose=verbose)
			if ( ! is.null(fitAns)) cat( "  ", i, fids[i], "RMS.Dev=", round(fitAns$RMS.Deviation,digits=3))
			return( fitAns)
		})

	rmsDev <- r2cod <- r2p <- pval <- vector( length=NS)
	for ( i in 1:NS) {
		# catch case of a single element affect lapply
		fitAns <- if ( NS == 1) mcAns else mcAns[[i]]
		if ( is.null(fitAns)) next
		ans[ , i] <- fitAns$BestFit
		rmsDev[i] <- fitAns$RMS.Deviation
		r2cod[i] <- fitAns$R2.CoD
		r2p[i] <- fitAns$R2.Pearson
		pval[i] <- fitAns$Pvalue
	}
	out2 <- data.frame( "SampleID"=fids, "RMS.Deviation"=rmsDev, "R2.CoD"=r2cod, "R2.Pearson"=r2p,
				"P.Value"=pval, stringsAsFactors=F)
	cat( "  Done.\n")

	# standardize the values to 100 %
	# NO, let's not do this here anymore.  Allow original values to be passed back
	pcts <- ans
	# for ( i in 1:NS) pcts[ , i] <- ans[ , i] * 100 / sum( ans[ , i])

	if (plot) {
		if ( label == "") label <- paste( "Sample =", fids[1], "  Algorithm =", algorithm, "  Log =", useLog)
		plotTranscriptProportions(pcts, col=plot.col, label=label)
		# new plot printing method does not need explicit extension
		plotFile <- paste( getCurrentSpeciesFilePrefix(), getCellTypeReference(), "DeconvProportions", algorithm, sep=".")
		if ( length(fids) == 1) plotFile <- paste( fids[1], plotFile, sep=".")
		plotFile <- file.path( plot.path, plotFile)
		printPlot( plotFile)
	}

	return( list( "BestFit"=pcts, "Statistics"=out2))
}


`matrix.TranscriptDeconvolution` <- function( m, targetM=getTranscriptDeconvolutionTargetMatrix(), 
					useLog=FALSE, minIntensity=0, 
					arrayFloorIntensity=NULL, dropLowVarianceGenes=NULL, geneUniverse=NULL,
					algorithm=c("port","default","plinear","LM","GenSA","steepDescent"),
					plot=TRUE, plot.path=".", plot.col=NULL, label="", dropRiboClearGenes=TRUE, verbose=TRUE) {

	NS <- ncol(m)
	ND <- ncol(targetM)
	geneIDs <- rownames(m)
	if ( is.null( geneIDs)) stop( "Matrix has no GeneID rownames..")

	ans <- matrix( 0, nrow=ND, ncol=NS)
	colnames(ans) <- fids <- colnames(m)
	rownames(ans) <- colnames(targetM)

	algorithm <- match.arg( algorithm)
	if (algorithm == "GenSA") require( GenSA)
	if (algorithm == "LM") require( minpack.lm)

	# try to do in multi core
	mcAns <- multicore.lapply( 1:NS, FUN=function(i) {
			v <- m[ ,i]
			tbl <- data.frame( "GENE_ID"=geneIDs, "INTENSITY"=v, stringsAsFactors=F)
			if (i > 1) verbose <<- FALSE
			fitAns <- fit.transcriptBlend( tbl, targetM, geneColumn="GENE_ID", 
					intensityColumn="INTENSITY", useLog=useLog, 
					minIntensity=minIntensity, 
					arrayFloorIntensity=arrayFloorIntensity, 
					dropLowVarianceGenes=dropLowVarianceGenes, 
					geneUniverse=geneUniverse, dropRiboClearGenes=dropRiboClearGenes, 
					algorithm=algorithm, verbose=verbose)
			if ( ! is.null(fitAns)) cat( "  ", i, fids[i], "RMS.Dev=", round(fitAns$RMS.Deviation,digits=3))
			return( fitAns)
		})

	rmsDev <- r2cod <- r2p <- pval <- vector( length=NS)
	for ( i in 1:NS) {
		# catch case of a single element affect lapply
		fitAns <- if ( NS == 1) mcAns else mcAns[[i]]
		if ( is.null(fitAns)) next
		ans[ , i] <- fitAns$BestFit
		rmsDev[i] <- fitAns$RMS.Deviation
		r2cod[i] <- fitAns$R2.CoD
		r2p[i] <- fitAns$R2.Pearson
		pval[i] <- fitAns$Pvalue
	}
	out2 <- data.frame( "SampleID"=fids, "RMS.Deviation"=rmsDev, "R2.CoD"=r2cod, "R2.Pearson"=r2p,
				"P.Value"=pval, stringsAsFactors=F)
	cat( "  Done.\n")

	# standardize the values to 100 %
	pcts <- ans
	for ( i in 1:NS) pcts[ , i] <- ans[ , i] * 100 / sum( ans[ , i])

	if (plot) {
		if ( label == "") label <- paste( "Sample =", colnames(m)[1], "  Algorithm =", algorithm, "  Log =", useLog)
		plotTranscriptProportions(pcts, col=plot.col, label=label)
		plotFile <- paste( getCurrentSpeciesFilePrefix(), getCellTypeReference(), "DeconvProportions", algorithm, sep=".")
		if ( length(fids) == 1) plotFile <- paste( fids[1], plotFile, sep=".")
		plotFile <- file.path( plot.path, plotFile)
		printPlot( plotFile)
	}

	return( list( "BestFit"=pcts, "Statistics"=out2))
}


`plotTranscriptProportions` <- function( pcts, mode=c("bars","auto","pie","lines"), 
					col=rainbow( nrow(pcts), end=0.82), label="", useLog=FALSE, min.value.show=0.2,
					pch=22, pt.cex=2.5, lwd=4, legend.cex=0.9, label.cex=0.9, cex.axis=1,
					text.rotation=if( ncol(pcts) < 9) 0 else 90, 
					gaps=NULL, significance.values=NULL, significance.scaling=FALSE, ...) {

	# evaluate all files, and then decide how to visualize
	NS <- ncol(pcts)
	ND <- nrow(pcts)
	mode <- match.arg( mode)
	fids <- colnames(pcts)

	# previously, the proportions were alway true percentages, that sum to 100%.
	# No longer a given.  Force it now
	for ( i in 1:NS) pcts[ ,i] <- pcts[ ,i] * 100 / sum( pcts[ ,i], na.rm=T)

	doPie <- (( mode == "pie") || (mode == "auto" && NS < 3))
	doBars <- (( mode == "bars") || (mode == "auto" && NS > 2))
	doLines <- (( mode == "lines") || (mode == "auto" && NS > 16))

	# visuals...
	if ( doPie) {
		par( mfcol=c(1,1))
		# set up with useful sizes for this number of pies
		maintext <- "Transcript Proportions:     "
		saveMAI <- par( 'mai')
		if (NS > 1) par( mfrow=c(1,2))
		if (NS > 2) par( mfrow=c(2,2))
		if (NS > 4) par( mfrow=c(2,3))
		if (NS > 6) par( mfrow=c(3,3))
		if (NS > 9) par( mfrow=c(4,4))
		if (NS > 3) maintext <- ""
		if (NS > 1 && all( saveMAI == c(1.02,0.82,0.82,0.42))) {
			par( mai=c( 0.15, 0.4, 0.6, 0.4))
			on.exit( par( 'mai'=saveMAI))
		}

		for ( i in 1:NS) {
			v <- pcts[ ,i]
			nams <- paste( rownames(pcts), as.percent( v, big.value=100, digits=1), sep="  ")
			toShow <- which( v > min.value.show)
			pie( v[toShow], labels=nams[toShow], col=col[toShow], main=paste( maintext, fids[i]), ...)
		}
		par( mfcol=c(1,1))
		dev.flush()
		return()
	}

	LAS <- if ( NS < 5) 1 else 3
	X_LIM_SCALE <- 1.3
	nextra <- length(gaps)

	if ( doLines) {
		# more than 2, so try colored points
		par( mfcol=c(1,1))
		xLimits <- c( 0, NS + 1.15 + ( 0.3 * (NS-1)))
		yLimits <- range( pcts)
		yLabel <- "Proportion per Component" 
		Log=""
		# perhaps we were given some P values, to weight how we draw things
		myCEX <- rep.int( pt.cex, ND)
		myLWD <- rep.int( lwd, ND)
		myTXTCEX <- rep.int( label.cex, ND)
		myTXTCOL <- rep.int( 'black', ND)
		if ( ! is.null( significance.values)) {
			whSig <- match( rownames(pcts), names(significance.values), nomatch=0)
			mySignif <- rep.int( 0.001, ND)
			mySignif[ whSig > 0] <- significance.values[whSig]
			# best significance keep the given weight/color scheme, worse get less
			if (significance.scaling) {
				whNow <- which( mySignif > 0.01)
				myCEX[whNow] <- myCEX[whNow] * 0.7
				myLWD[whNow] <- myLWD[whNow]  * 0.6
				myTXTCEX[whNow] <- myTXTCEX[whNow] * 0.8
				myTXTCOL[whNow] <- adjustColor( myTXTCOL[whNow], 0.2)
				whNow <- which( mySignif > 0.05)
				myCEX[whNow] <- myCEX[whNow] * 0.7
				myLWD[whNow] <- myLWD[whNow]  * 0.6
				myTXTCEX[whNow] <- myTXTCEX[whNow] * 0.8
				myTXTCOL[whNow] <- adjustColor( myTXTCOL[whNow], 0.2)
				whNow <- which( mySignif > 0.1)
				myCEX[whNow] <- myCEX[whNow] * 0.7
				myLWD[whNow] <- myLWD[whNow]  * 0.6
				myTXTCEX[whNow] <- myTXTCEX[whNow] * 0.8
				myTXTCOL[whNow] <- adjustColor( myTXTCOL[whNow], 0.2)
				whNow <- which( mySignif > 0.2)
				myCEX[whNow] <- myCEX[whNow] * 0.7
				myLWD[whNow] <- myLWD[whNow]  * 0.6
				myTXTCEX[whNow] <- myTXTCEX[whNow] * 0.8
				myTXTCOL[whNow] <- adjustColor( myTXTCOL[whNow], 0.2)
			}
		}

		if (useLog) {
			Log <- 'y'
			yLimits[1] <- min.value.show
			pcts[ pcts < min.value.show] <- min.value.show
			yLabel <- "Proportion per Component  (log scale)" 
		}
		plot( 1, 1, type="n", main=paste( "Transcriptome Proportions:   ", label), ylab=yLabel,
			xlim=xLimits, ylim=yLimits, font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1, 
			xaxt="n", xlab=NA, log=Log, ...)
		axis( 1, at=1:NS, colnames(pcts), las=LAS, font=2, cex.axis=cex.axis)
		for ( i in 1:NS) {
			v <- pcts[ ,i]
			toShow <- which( v >= min.value.show)
			ord <- order( v[toShow])
			points( rep.int(i,length(toShow)), v[toShow[ord]], bg=col[toShow[ord]], pch=pch, cex=myCEX[toShow[ord]])
		}
		# if gaps are wanted, that inserts breaks (NAs) in the lines
		nextra <- 0
		linePts <- 1:NS
		if ( ! is.null(gaps)) {
			gaps <- sort( unique( as.integer( gaps)))
			for (g in gaps) {
				if ( is.na(g)) next
				if ( g < 1) next
				if ( g+nextra >= ncol(pcts)) next
				pcts <- cbind( pcts[,1:(g+nextra)], rep.int(NA,nrow(pcts)), pcts[,(g+nextra+1):ncol(pcts)])
				linePts <- c( linePts[1:(g+nextra)], NA, linePts[(g+nextra+1):length(linePts)])
				nextra <- nextra + 1
			}
			NSpts <- ncol(pcts)
		}
		# try to draw the line in order of pct
		pctMeans <- apply( pcts, 1, mean, na.rm=T)
		lineOrd <- order( pctMeans)
		for( j in lineOrd) {
			if ( all(is.na(pcts[j,])) || max(pcts[j,],na.rm=T) < min.value.show) next
			lines( linePts, pcts[ j, ], col=col[j], lwd=myLWD[j])
			text( 0.9, pcts[j,1], rownames(pcts)[j], pos=2, cex=myTXTCEX[j], col=myTXTCOL[j])
			text( NS+0.1, pcts[j,NS+nextra], rownames(pcts)[j], pos=4, cex=myTXTCEX[j], col=myTXTCOL[j])
		}
		# the  line plot has no fixed order, so don't reverese the colors
		legend( "topright", rownames(pcts), fill=col, bg='white', cex=legend.cex)
		dev.flush()
		return()
	}
	if ( doBars) {
		# more than 2, so try a barchart
		par( mfcol=c(1,1))
		nextra <- 0
		barAts <- 1:NS
		if ( ! is.null(gaps)) {
			gaps <- sort( unique( as.integer( gaps)))
			for (g in gaps) {
				if ( is.na(g)) next
				if ( g < 1) next
				if ( g+nextra >= ncol(pcts)) next
				pcts <- cbind( pcts[,1:(g+nextra)], rep.int(NA,nrow(pcts)), pcts[,(g+nextra+1):ncol(pcts)])
				barAts <- c( barAts[1:g], (barAts[(g+1):length(barAts)] + 1))
				nextra <- nextra + 1
			}
		}
		xLimits <- c( 0, max(barAts) + 1.15 + ( 0.3 * (NS+length(gaps)-1)))
		ans <- barplot( pcts, col=col, main=paste( "Transcriptome Proportions:   ", label), 
			ylab="Proportion per Component", las=LAS, 
			xlim=xLimits, font.axis=2, font.lab=2, cex.axis=cex.axis, cex.lab=1.1, 
			legend=TRUE, args.legend=list( "cex"=legend.cex, "bg"='white'), ...)

		# calc the center of each band, to send back to caller
		bandCtr <- matrix( 0, nrow(pcts), ncol(pcts))
		rownames(bandCtr) <- rownames(pcts)
		colnames(bandCtr) <- colnames(pcts)
		for ( k in 1:ncol(pcts)) {
			myTops <- cumsum( pcts[,k])
			myDiffs <- diff( c( 0, myTops))
			bandCtr[,k] <- myTops - (myDiffs/2)
		}

		# draw some labels where we can
		samplesToLabel <- barAts
		if ( text.rotation == 0) {
			if ( NS > 16) samplesToLabel <- samplesToLabel[ seq( 1, NS, by=2)]
			if ( NS > 32) samplesToLabel <- samplesToLabel[ seq( 1, NS, by=3)]
			if ( NS > 48) samplesToLabel <- samplesToLabel[ seq( 1, NS, by=4)]
		}
		rotate0.sizeCut <- 3
		rotate90.sizeCut <- 12
		for ( i in samplesToLabel) {
			v <- pcts[ ,i]
			if ( all( is.na( v))) next
			y0 <- 0
			xExtra <- 0.01
			for ( j in 1:ND) {
				yNow <- v[j]
				if ( is.na( yNow)) next
				if ( text.rotation == 90 && yNow >= rotate90.sizeCut) {
					text( ans[i]+xExtra, (y0+yNow/2), rownames(pcts)[j], col=1, cex=label.cex, 
						srt=90)
					xExtra <- 0.2 - xExtra
				} 
				if ( text.rotation == 0 && yNow >= rotate0.sizeCut) {
					text( ans[i], (y0+yNow/2), rownames(pcts)[j], col=1, cex=label.cex, srt=0)
				}
				y0 <- y0 + yNow
			}
		}
		dev.flush()

		# package up a bit of context to send back
		barAts <- ans
		names(barAts) <- colnames(pcts)
		return( invisible( list( "bar.centers"=barAts, "band.centers"=bandCtr)))
	}
}


`compareTranscriptProportions` <- function( pcts, groups, levels=sort(unique(as.character(groups))), 
					col=rainbow(nrow(pcts),end=0.8),
					minPerGroup=3, test=t.test, plot=TRUE, plot.path=".", 
					plot.mode=c("bars","auto","pie","lines"), label="", useLog=FALSE,
					wt.fold=1, wt.pvalue=2, min.percent=0.1, gaps=NULL, 
					stats.color=NULL, significance.scaling=FALSE, ...) {

	NC <- ncol(pcts)
	NR <- nrow(pcts)

	if ( length(groups) != NC) stop( "Length of 'groups' must match columns of 'pcts'")
	grpFac <- factor( groups, levels=levels)
	grpPtrs <- tapply( 1:NC, grpFac, FUN=NULL)
	NG <- nlevels( grpFac)
	grpLabels <- base::levels(grpFac)
	if ( ! all( 1:NG %in% grpPtrs)) stop( "Some group levels have zero members")

	# previously, we were gauranteed that all proportions summed to 100%.  No longer always true.
	# But we need them to be for comparitive plotting.  So do it now with a warning
	pctSums <- apply( pcts, 2, sum, na.rm=T)
	if ( any( round(pctSums) != 100)) {
		cat( "\n  Note: normalizing Proportions to sum to 100% for plotting.")
		for ( j in 1:NC) pcts[ ,j] <- pcts[ ,j] * 100 / pctSums[j]
	}

	bigOut <- vector( mode="list")
	nOut <- 0

	# do all 2-way compares we can
	for ( j1 in 1:(NG-1)) {
		isGrp1 <- which( grpPtrs == j1)
		#if ( length( isGrp1) < minPerGroup) next
		label1 <- grpLabels[ j1]
		for ( j2 in (j1+1):NG) {
			isGrp2 <- which( grpPtrs == j2)
			#if ( length( isGrp2) < minPerGroup) next
			label2 <- grpLabels[ j2]

			# OK to do these 2 groups
			v1 <- v2 <- pval <- vector( length=NR)
			for ( i in 1:NR) {
				x <- pcts[ i, isGrp1]
				y <- pcts[ i, isGrp2]
				v1[i] <- mean( x, na.rm=T)
				v2[i] <- mean( y, na.rm=T)
				# values too close to zero can give great P-values.  Override that
				if ( all( c(x,y) < min.percent)) {
					pval[i] <- 1.0
				} else {
					# make sure we have enough values
					if ( (nNow <- length(x)) < minPerGroup) x <- c( x, jitter(sample(x,size=minPerGroup-nNow, replace=T)))
					if ( (nNow <- length(y)) < minPerGroup) y <- c( y, jitter(sample(y,size=minPerGroup-nNow, replace=T)))
					isXzero <- which( x < min.percent)
					if (Nzero <- length(isXzero)) x[ isXzero] <- runif( Nzero, min.percent/10, min.percent)
					isYzero <- which( y < min.percent)
					if (Nzero <- length(isYzero)) y[ isYzero] <- runif( Nzero, min.percent/10, min.percent)
					if ( is.function(test)) {
						ans <- test( x, y)
						pval[i] <- ans$p.value
					} else {
						pval[i] <- 1
					}
				}
			}
			# report the fold as Group2 over Group1
			fold <- log2( (v2+min.percent) / (v1+min.percent))
			pval[ is.nan(pval)] <- 1.0
			pval[ is.na(pval)] <- 1.0
			signif <- rep.int( "", length(pval))
			signif[ pval < 0.1] <- "."
			signif[ pval < 0.05] <- "*"
			signif[ pval < 0.01] <- "**"
			signif[ pval < 0.001] <- "***"
		
			out <- data.frame( "Component"=rownames(pcts), "N1"=length(isGrp1), "N2"=length(isGrp2), 
					"Avg1"=round(v1,digits=2), "Avg2"=round(v2,digits=2), 
					"Log2Fold"=round(fold,digits=3), "PI.value"=round( piValue(fold,pval), digits=3), 
					"P.value"=round(pval,digits=5), "Signif"=signif, stringsAsFactors=F)
			colnames(out)[2:3] <- paste( "N", c( label1, label2), sep="_")
			colnames(out)[4:5] <- paste( "Avg", c( label1, label2), sep="_")
			ord <- diffExpressRankOrder( out$Log2Fold, out$P.value, wt.fold=wt.fold, wt.pvalue=wt.pvalue)
			out <- out[ ord, ]
			rownames(out) <- 1:nrow(out)

			nOut <- nOut + 1
			bigOut[[nOut]] <- out
			names( bigOut)[nOut] <- paste( label1, label2, sep=".vs.")
		}
	}

	# we have the compares done, draw it too

	# try to pass info about how significant each was down to the plot tool, to weight/color accordingly
	comps <- bigOut[[1]]$Component
	nComps <- length(comps)
	bestPvals <- rep.int( 1, nComps)
	iList <- 0
	for ( i1 in 1:(NG-1)) for (i2 in (i1+1):NG) {
		iList <- iList + 1
		if ( i2 - i1 != 1) next  # not an adjacent pair
		# don't use stats where the caller put gaps in the plot
		if ( ! is.null(gaps) && (i1 %in% gaps)) next
		smlDF <- bigOut[[iList]]
		wh <- match( comps, smlDF$Component)
		thisP <- smlDF$P.value[wh]
		bestPvals <- pmin( bestPvals, thisP, na.rm=T)
	}
	names(bestPvals) <- comps

	if (plot) {
		plot.mode <- match.arg( plot.mode)
		pctsGrp <- t( apply( pcts, 1, function(x) tapply( x, grpFac, mean)))
		# show the group counts, if not too many
		if (NG <=12) colnames(pctsGrp) <- paste( colnames(pctsGrp), "\n(N=", as.numeric(table(as.numeric(grpPtrs))), ")", sep="")
		plotAns <- plotTranscriptProportions( pctsGrp, label=label, col=col, mode=plot.mode, gaps=gaps, useLog=useLog, 
							significance.values=bestPvals, significance.scaling=significance.scaling, ...)

		# perhaps try to show the significance? 
		# various code for some of the modes
		if ( plot.mode == "bars" && NG == 2) {
			myStats <- bigOut[[1]]
			myBarCenters <- plotAns$bar.centers
			myBandCenters <- plotAns$band.centers
			toShow <- which( myStats$Signif != "")
			if ( length(toShow)) {
				where <- match( myStats$Component[toShow], rownames(myBandCenters))
				textToShow <- myStats$Signif[toShow]
				upDownCall <- myStats$Log2Fold[toShow]
				myX <- ifelse( upDownCall > 0, myBarCenters[2] + 0.6, myBarCenters[1] - 0.6)
				myY <- ifelse( upDownCall > 0, myBandCenters[ where, 2], myBandCenters[ where, 1])
				statsColor <- col[where]
				if ( ! is.null(stats.color)) statsColor <- stats.color[1]
				text( myX, myY, textToShow, col=statsColor, cex=1.95)
			}
		}
		if ( plot.mode == "lines") {
			# we can show the stats between adjacent groups, step along the list of stats
			smallX <- if ( NG > 2) 0.15 else -0.1
			if ( NG > 8) smallX <- 0.28
			smallY <- diff( range( as.vector( pctsGrp), na.rm=T)) * 0.01
			if ( useLog) smallY <- 0
			iList <- 0
			for ( i1 in 1:(NG-1)) for (i2 in (i1+1):NG) {
				iList <- iList + 1
				if ( i2 - i1 != 1) next  # not an adjacent pair
				# don'thsow stats where the caller put gaps in the plot
				if ( ! is.null(gaps) && (i1 %in% gaps)) next
				myStats <- bigOut[[iList]]
				toShow <- which( myStats$Signif != "")
				if ( length(toShow)) {
					where <- match( myStats$Component[toShow], rownames(pctsGrp))
					textToShow <- myStats$Signif[toShow]
					upDownCall <- myStats$Log2Fold[toShow]
					myX <-rep.int( i2 - smallX, length(toShow))
					myY <-ifelse( upDownCall > 0, pctsGrp[where,i2]+smallY, pctsGrp[where,i2]-smallY)
					statsColor <- col[where]
					if ( ! is.null(stats.color)) statsColor <- stats.color[1]
					text( myX, myY, textToShow, col=statsColor, cex=2, pos=2, offset=0.15)
				}
			}
		}
	}

	return( if (nOut == 1) bigOut[[1]] else bigOut)
}


`mergeTranscriptProportions` <- function( pcts, groups, drops=NULL) {

	# generic way to reduce the number of cel types by merging similar types
	if ( length( groups) != nrow(pcts)) stop( "'groups' must be same length as 'nrow(pcts)'")

	out <- data.frame()
	tapply( 1:nrow(pcts), grpFac <- factor(groups), function(x) {
				sml <- pcts[ x, , drop=F]
				smlDF <- as.data.frame(sml[1, , drop=F])
				vals <- apply( sml, 2, sum, na.rm=T)
				smlDF[ 1, ] <- vals
				out <<- rbind( out, smlDF)
			})

	out <- as.matrix( out)
	rownames(out) <- levels(grpFac)
	if ( ! is.null(drops)) {
		whoDrop <- which( rownames(out) %in% drops)
		if ( length(whoDrop)) out <- out[ -whoDrop, ]
	}

	return( out)
}

