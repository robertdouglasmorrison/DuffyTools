# transcriptDeconvolution.R -- higher level code to do transcript deconvolution 
#				uses lower level 'transcriptBlend' tool


`getTranscriptDeconvolutionTargetMatrix` <- function( speciesID=getCurrentSpecies(), target=NULL) {

	# grapb the 'most appropriate' target
	if ( speciesID %in% PARASITE_SPECIES) {
		targetM <- getLifeCycleMatrix()
		# the life cycle data is all genes
		# while a typical target is the higher expressing subset of genes
		bigV <- apply( targetM, 1, max, na.rm=T)
		keep <- which( bigV >= 5.0)  # RPKM threshold of detection in PF
		targetM <- targetM[ keep, ]
		return( targetM)
	}

	if ( speciesID != getCurrentSpecies()) {
		oldSpecies <- getCurrentSpecies()
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( oldSpecies))
	}
	myPrefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( target)) target <- "ImmuneCell.TargetMatrix"

	f <- paste( myPrefix, target, sep=".")
	targetM <- NULL
	data( list=f, package="DuffyTools", envir=environment())
	if ( is.null( targetM)) cat( "\nFailed to load target data.  Tried: ", f)
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
					sep="\t", useLog=FALSE, normalize=FALSE, minIntensity=0, 
					arrayFloorIntensity=NULL, dropLowVarianceGenes=NULL,
					algorithm=c("port","default","plinear","LM","GenSA","steepDescent"),
					startFractions=NULL, plot=TRUE, plot.path=".", plot.col=NULL, 
					label="", verbose=TRUE) {

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
			# allow suing specific starting percentages
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
					normalize=normalize, minIntensity=minIntensity, 
					arrayFloorIntensity=arrayFloorIntensity, 
					dropLowVarianceGenes=dropLowVarianceGenes, 
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
	pcts <- ans
	for ( i in 1:NS) pcts[ , i] <- ans[ , i] * 100 / sum( ans[ , i])

	if (plot) {
		if ( label == "") label <- paste( "Sample =", fids[1], "  Algorithm =", algorithm, "  Log =", useLog)
		plotTranscriptProportions(pcts, col=plot.col, label=label)
		plotFile <- paste( "TranscriptProportions", algorithm, "png", sep=".")
		if ( length(fids) == 1) plotFile <- paste( fids[1], if (useLog) "YesLog" else "NoLog", plotFile, sep=".")
		plotFile <- file.path( plot.path, plotFile)
		dev.print( png, plotFile, width=900)
	}

	return( list( "BestFit"=pcts, "Statistics"=out2))
}


`matrix.TranscriptDeconvolution` <- function( m, targetM=getTranscriptDeconvolutionTargetMatrix(), 
					useLog=FALSE, normalize=FALSE, minIntensity=0, 
					arrayFloorIntensity=NULL, dropLowVarianceGenes=NULL,
					algorithm=c("port","default","plinear","LM","GenSA","steepDescent"),
					plot=TRUE, plot.path=".", plot.col=NULL, label="", verbose=TRUE) {

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
					normalize=normalize, minIntensity=minIntensity, 
					arrayFloorIntensity=arrayFloorIntensity, 
					dropLowVarianceGenes=dropLowVarianceGenes, 
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
		plotFile <- paste( "TranscriptProportions", algorithm, "png", sep=".")
		if ( ncol(m) == 1) plotFile <- paste( colnames(m)[1], if (useLog) "YesLog" else "NoLog", plotFile, sep=".")
		plotFile <- file.path( plot.path, plotFile)
		dev.print( png, plotFile, width=900)
	}

	return( list( "BestFit"=pcts, "Statistics"=out2))
}


`plotTranscriptProportions` <- function( pcts, mode=c("bars","auto","pie","lines"), 
					col=rainbow( nrow(pcts), end=0.82), label="", 
					pch=22, pt.cex=2.5, lwd=4, legend.cex=0.9, label.cex=0.9, 
					text.rotation=if( ncol(pcts) < 9) 0 else 90, 
					gaps=NULL, ...) {

	# evaluate all files, and then decide how to visualize
	NS <- ncol(pcts)
	ND <- nrow(pcts)
	mode <- match.arg( mode)
	fids <- colnames(pcts)

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
			toShow <- which( v > 0.25)
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
		ylim <- range( pcts)
		plot( 1, 1, type="n", main=paste( "Transcript Proportions:   ", label), 
			ylab="Proportion per Component", 
			xlim=c(0.2,(NS+nextra)*X_LIM_SCALE+0.5), ylim=ylim, font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1, 
			xaxt="n", xlab=NA, ...)
		axis( 1, at=1:NS, colnames(pcts), las=LAS, font=2, cex=1.05)
		for ( i in 1:NS) {
			v <- pcts[ ,i]
			toShow <- which( v > 0.25)
			ord <- order( v[toShow])
			points( rep.int(i,length(toShow)), v[toShow[ord]], bg=col[toShow[ord]], pch=pch, cex=pt.cex)
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

		for( j in 1:ND) {
			if ( max( pcts[ j, ], na.rm=T) < 2) next
			lines( linePts, pcts[ j, ], col=col[j], lwd=lwd)
			text( 0.9, pcts[j,1], rownames(pcts)[j], pos=2, cex=label.cex)
			text( NS+0.1, pcts[j,NS+nextra], rownames(pcts)[j], pos=4, cex=label.cex)
		}

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
		ans <- barplot( pcts, col=col, main=paste( "Transcript Proportions:   ", label), 
			ylab="Proportion per Component", las=LAS, 
			xlim=c(0.2,(NS+nextra+1)*X_LIM_SCALE+0.5), font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1, 
			legend=TRUE, args.legend=list( "cex"=legend.cex, "bg"='white'), ...)
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
		return()
	}
}


`compareTranscriptProportions` <- function( pcts, groups, levels=sort(unique(as.character(groups))), 
					minPerGroup=3, test=t.test, plot=TRUE, plot.path=".", 
					plot.mode=c("bars","auto","pie","lines"), label="", 
					wt.fold=1, wt.pvalue=2, min.percent=0.1, ...) {

	NC <- ncol(pcts)
	NR <- nrow(pcts)

	if ( length(groups) != NC) stop( "Length of 'groups' must match columns of 'pcts'")
	grpFac <- factor( groups, levels=levels)
	grpPtrs <- tapply( 1:NC, grpFac, FUN=NULL)
	NG <- nlevels( grpFac)
	grpLabels <- base::levels(grpFac)
	if ( ! all( 1:NG %in% grpPtrs)) stop( "Some group levels have zero members")

	bigOut <- vector( mode="list")
	nOut <- 0

	# do all 2-way compares we can
	for ( j1 in 1:(NG-1)) {
		isGrp1 <- which( grpPtrs == j1)
		if ( length( isGrp1) < minPerGroup) next
		label1 <- grpLabels[ j1]
		for ( j2 in (j1+1):NG) {
			isGrp2 <- which( grpPtrs == j2)
			if ( length( isGrp2) < minPerGroup) next
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
		
			out <- data.frame( "Component"=rownames(pcts), "N1"=length(isGrp1), "N2"=length(isGrp2), 
					"Avg1"=round(v1,digits=2), "Avg2"=round(v2,digits=2), 
					"Log2Fold"=round(fold,digits=3), "P.value"=round(pval,digits=5), 
					"PI.value"=round( piValue(fold,pval), digits=3), stringsAsFactors=F)
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
	if (plot) {
		plot.mode <- match.arg( plot.mode)
		pctsGrp <- t( apply( pcts, 1, function(x) tapply( x, grpFac, mean)))
		# show the group counts, if not too many
		if (NG <=12) colnames(pctsGrp) <- paste( colnames(pctsGrp), "\n(N=", as.numeric(table(as.numeric(grpPtrs))), ")", sep="")
		plotTranscriptProportions( pctsGrp, label=label, mode=plot.mode, ...)
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

