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


`eosinophilExpressionComparison` <- function( x, y, value.mode=c("absolute","relative"), col=c('tan','grey70'), 
						labels=c("x","y"), sep="\t", main=NULL, legend="topright", 
						...) {

	isMAT <- isVEC <- FALSE
	out <- NULL
	value.mode <- match.arg( value.mode)

	if ( is.null(main)) main <- "Esosinophil Marker Gene Expression:"

	xAns <- eosinophilExpression( x, value.mode=value.mode, sep=sep)
	if ( is.null( xAns)) return(NULL)
	yAns <- eosinophilExpression( y, value.mode=value.mode, sep=sep)
	if ( is.null( yAns)) return(NULL)

	# try to compare and display in a sensible manner
	if ( is.vector(xAns) && is.vector(yAns)) {
		xGenes <- names(xAns)
		yGenes <- names(yAns)
		xyGenes <- sort( intersect( xGenes, yGenes))
		NG <- length(xyGenes)
		if ( ! NG) {
			cat( "\nNo Eosinophil genes in intersect of both data sets.")
			return(NULL)
		}
		xAns <- xAns[ match( xyGenes, xGenes)]
		yAns <- yAns[ match( xyGenes, yGenes)]
		legend.label.suffix <- ""
	} else if ( is.matrix(xAns) && is.matrix(yAns)) {
		xGenes <- rownames(xAns)
		yGenes <- rownames(yAns)
		xyGenes <- sort( intersect( xGenes, yGenes))
		NG <- length(xyGenes)
		if ( ! NG) {
			cat( "\nNo Eosinophil genes in intersect of both data sets.")
			return(NULL)
		}
		xAnsM <- xAns[ match( xyGenes, xGenes), ]
		yAnsM <- yAns[ match( xyGenes, yGenes), ]
		xAns <- apply( xAnsM, 1, mean)
		yAns <- apply( yAnsM, 1, mean)
		legend.label.suffix <- paste( " (N=", c(ncol(xAnsM),ncol(yAnsM)), ")", sep="")
	} else {
		cat( "\nError: both x and y must be of the same type")
		return(NULL)
	}

	# compare the 2 sets of values
	testAns <- wilcox.test( xAns, yAns, paired=T)
	xMean <- mean(xAns)
	yMean <- mean(yAns)
	fc <- log2( (yMean+1) / (xMean+1))

	# reformat to plot the results
	m <- matrix( c(xAns,yAns), nrow=length(xyGenes), ncol=2)
	rownames(m) <- xyGenes
	colnames(m) <- labels
	plotAns <- barplot( t(m), beside=TRUE, col=col, las=3, ylab="Gene Expression", xlab=NA, 
				ylim=c(0,max(m)*1.05), main=main, ...)

	lines( c(0.5,NG*3+0.5), c(xMean,xMean), col=col[1], lty=2, lwd=2)
	lines( c(0.5,NG*3+0.5), c(yMean,yMean), col=col[2], lty=2, lwd=2)

	if ( ! is.na(legend)) {
		if ( is.null(legend) || legend == "") legend <- "topright"
		graphics::legend( legend, paste(labels,legend.label.suffix), fill=col, bg='white')
	}

	out <- list( "x.avg"=xMean, "y.avg"=yMean, "Log2Fold"=fc, "p.value"=testAns$p.value)
	return( out)
}

