# eosinophilTools.R -- try to capture info from eosinophil marker genes



# extract expression values from a variety of object types
`eosinophilExpression` <- function( x, value.mode=c("absolute","relative"), sep="\t") {

	eosinMarkerGenes <- c(  "FRRS1", "ADORA3", "ALOX15", "CCL23", "AFF2", 
				"CEBPE", "AOC1", "TFF3", "PRSS33", "CCL15", 
				"KHDRBS3", "VSTM1", "IDO1", "CD24", "OLIG2")
	eosinRPKMvalues <- c(    61.4777, 100.8599, 284.4524, 94.6696, 20.9992, 
				128.5605, 24.3029, 32.6940, 86.2324, 8.9697, 
				6.6057, 128.8078, 206.893, 88.6310, 21.5762)

	isMAT <- isDF <- isVEC <- FALSE
	out <- NULL

	if ( is.character(x)) {
		x  <- read.delim( x[1], as.is=T, sep=sep)
	}

	if ( is.matrix(x)) {
		geneNames <- shortGeneName( rownames(x), keep=1)
		values <- x
		isMAT <- TRUE
	} else if ( is.data.frame( x)) {
		gidColumn <- which( toupper(colnames(x)) %in% c( "GENE_ID", "GENEID", "GENE"))[1]
		if ( is.na( gidColumn)) {
			cat( "\nUnable to find 'GENE_ID' column in data frame.  \nEncountered: ", colnames(x))
			return( NULL)
		}
		intenColumn <- which( toupper(colnames(x)) %in% c( "RPKM_M", "INTENSITY", "TPM_M"))[1]
		if ( is.na( intenColumn)) {
			cat( "\nUnable to find 'RPKM_M' column in data frame.  \nEncountered: ", colnames(x))
			return( NULL)
		}
		geneNames <- shortGeneName( x[[gidColumn]], keep=1) 
		values <- as.numeric( x[[intenColumn]])
		isDF <- TRUE
	} else {
		geneNames <- shortGeneName( names(x), keep=1)
		values <- as.numeric(x)
		isVEC <- TRUE
	}

	where <- match( eosinMarkerGenes, geneNames, nomatch=NA)
	if ( all( is.na( where))) {
		cat( "\nWarning:  failed to find any named eosinophil genes..")
		return( NULL)
	}

	value.mode <- match.arg( value.mode)
	if ( isVEC || isDF) {
		out <- values[ where]
		names(out) <- eosinMarkerGenes
		if ( value.mode == "relative") out <- out / eosinRPKMvalues
	} else if ( isMAT) {
		out <- values[ where, , drop=F]
		rownames(out) <- eosinMarkerGenes
		if ( value.mode == "relative") for ( j in 1:ncol(out)) out[,j] <- out[,j] / eosinRPKMvalues
	}
	return( out)
}


`eosinophilExpressionComparison` <- function( x, groups=names(x), levels=sort(unique(groups)), 
						value.mode=c("absolute","relative"), col=NULL,
						sep="\t", main="", legend="topright", AVG.FUN=mean, ...) {

	isMAT <- isVEC <- FALSE
	out <- NULL
	value.mode <- match.arg( value.mode)
	mainText <- paste( "Esosinophil Marker Genes: ", main)

	# we may be given a list of filenames, or a matrix
	if ( is.character(x)) {
		# for files, extract the gene values from each
		files <- x
		N <- length(files)
		groups <- rep( groups, length.out=N)
		smlAns <- lapply( files, eosinophilExpression, value.mode=value.mode, sep=sep)
		drops <- sapply( smlAns, is.null)
		smlAns <- smlAns[ ! drops]
		files <- files[ ! drops]
		groups <- groups[ ! drops]
		if ( ! length(smlAns)) return(NULL)
		N <- length(smlAns)
		# visit each small set, to get the final set of genes
		allGenes <- names( smlAns[[1]])
		for ( j in 2:N) allGenes <- intersect( allGenes, names(smlAns[[j]]))
		NG <- length( allGenes)
		if ( ! NG) {
			cat( "\nNo Eosinophil genes in intersect of both data sets.")
			return(NULL)
		}
		allGenes <- sort( allGenes)
		m <- matrix( NA, nrow=NG, ncol=N)
		rownames(m) <- allGenes
		colnames(m) <- basename( files)
		for ( j in 1:N) {
			wh <- match( allGenes, names(smlAns[[j]]))
			m[ ,j] <- smlAns[[j]][wh]
		}
	} else if ( is.matrix(x)) {
		m <- eosinophilExpression( x, value.mode=value.mode)
		N <- ncol(m)
		NG <- nrow(m)
		groups <- rep( groups, length.out=N)
	}

	# now reduce by group
	grpFac <- factor(groups, levels=levels)
	uniqGrps <- levels(grpFac)
	NGRP <- length( uniqGrps)
	legend.label.suffix <- paste( " (N=", tapply(1:N,grpFac,length),")", sep="")
	if ( is.null(col)) {
		col <- rainbow( NGRP, end=0.75)
	} else {
		col <- rep( col, length.out=NGRP)
	}
	m2 <- matrix( NA, nrow=NG, ncol=NGRP)
	rownames(m2) <- rownames(m)
	colnames(m2) <- uniqGrps
	for ( i in 1:NG) m2[ i, ] <- tapply( m[ i, ], grpFac, AVG.FUN, na.rm=T)

	# now compute compares
	g1 <- g2 <- fc <- pv <- avg1 <- avg2 <- vector()
	nout <- 0
	avgAns <- apply( m2, 2, AVG.FUN, na.rm=T)
	for (i1 in 1:(NGRP-1)) for (i2 in (i1+1):NGRP) {
		xAns <- m2[ ,i1]
		yAns <- m2[ ,i2]
		testAns <- wilcox.test( xAns, yAns, paired=T)
		xMean <- avgAns[i1]
		yMean <- avgAns[i2]
		thisFC <- log2( (yMean+1) / (xMean+1))
		nout <- nout + 1
		g1[nout] <- uniqGrps[i1]
		g2[nout] <- uniqGrps[i2]
		avg1[nout] <- xMean
		avg2[nout] <- yMean
		fc[nout] <- thisFC
		pv[nout] <- testAns$p.value
	}

	# reformat to plot the results
	plotAns <- barplot( t(m2), beside=TRUE, col=col, las=3, ylab="Gene Expression", xlab=NA, 
				xlim=c(0,NG*(NGRP+1.75)), ylim=c(0,max(m2,na.rm=T)*1.05), main=mainText, ...)

	for ( j in 1:NGRP) lines( c(0.5,NG*(NGRP+1)+0.5), rep.int(avgAns[j],2), col=col[j], lty=3, lwd=2)

	if ( ! is.na(legend)) {
		if ( is.null(legend) || legend == "") legend <- "topright"
		graphics::legend( legend, paste( uniqGrps, legend.label.suffix), fill=col, bg='white')
	}

	out <- data.frame( "Group1"=g1, "Group2"=g2, "Average1"=avg1, "Average2"=avg2, "Log2Fold"=fc, "P.Value"=pv, stringsAsFactors=F)

	# try to show some P-values?
	mostSig <- which.min( out$P.Value)
	xLeftTic <- NG * (NGRP+1) + 1.5
	xRightTic <- xLeftTic * 1.02
	yTic1 <- out$Average1[mostSig]
	yTic2 <- out$Average2[mostSig]
	pShow <- out$P.Value[mostSig]
	pShowText <- if ( pShow > 0.0001) formatC( pShow, format="f", digits=4) else formatC( pShow, format="e", digits=2) 
	lines( c( xLeftTic, xRightTic, xRightTic, xLeftTic), c( yTic1, yTic1, yTic2, yTic2), lwd=1, col=1, lty=1)
	text( xRightTic, (yTic1+yTic2)/2, paste( "P =", pShowText), col=1, pos=4, cex=0.9)

	return( out)
}

