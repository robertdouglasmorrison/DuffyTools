# pcaTools.R -- make a PCA plots from various data

fileSet.PCAplot <- function( files, fids, geneColumn = "GENE_ID", 
				intensityColumn = "RPKM_M", sep = "\\t", useLog = FALSE, 
				d1=1, d2=2, main="", col=rainbow( length(fids), end=0.8), 
				pt.cex=2, pch=21, label.cex=0.9, na.mode=c("drop","zero"), verbose=FALSE, 
				plotOrder=1:length(files), ...) {

	m <- expressionFileSetToMatrix( files, fids, geneColumn=geneColumn, 
					intensityColumn=intensityColumn, verbose=verbose)

	if ( useLog) m <- log2( m + 1)
	na.mode <- match.arg( na.mode)

	out <- matrix.PCAplot( m, main=main, col=col, d1=d1, d2=d2, pt.cex=pt.cex, 
				label.cex=label.cex, pch=pch, na.mode=na.mode, 
				plotOrder=plotOrder, ...)

	return( out)
}


matrix.PCAplot <- function( m, main="", col=rainbow( ncol(m), end=0.8), d1=1, d2=2,
				pt.cex=2, pch=21, label.cex=0.9, na.mode=c("drop","zero"), 
				plotOrder=1:ncol(m), ...) {

	# drop any NA rows
	hasNA <- apply( m, 1, function(x) any( is.na( x)))
	if ( any( hasNA)) {
		na.mode <- match.arg( na.mode)
		if ( na.mode == "drop") {
			m <- m[ !hasNA, ]
			cat( "\nDropping rows with NA values: ", sum( hasNA))
		} else {
			nNA <- sum( is.na(m))
			m[ is.na(m)] <- 0
			cat( "\nZeroing entries with NA values: ", nNA)
		}
	}

	# do the PCA
	cat( "\nPCA using ", nrow(m), " observations per sample\n")
	pca <- prcomp( m)

	xValues <- pca$rotation[,d1]
	yValues <- pca$rotation[,d2]

	mainText <- paste( "PCA Plot:    ", main)
	xLim <- range( xValues)
	xLim <- xLim + diff(xLim)*c(-0.15,0.15)
	yLim <- range( yValues)
	yLim <- yLim + diff(yLim)*c(-0.1,0.1)

	# do a tiny jitter in case perfect overlap
	xShow <- jitter( xValues, factor=1.0)
	yShow <- jitter( yValues, factor=1.0)
	ord <- plotOrder

	plot( xShow[ord], yShow[ord], pch=pch, cex=pt.cex, bg=col[ord], main=mainText, xlim=xLim, ylim=yLim, 
			xlab=paste("Principal Component",d1), ylab=paste("Principal Component",d2),
			...)

	lbls <- colnames(m)
	require( plotrix)
	if ( any( label.cex > 0)) thigmophobe.labels( xShow, yShow, lab=lbls, cex=label.cex)
	dev.flush()

	out <- data.frame( "SampleID"=colnames(m), "PC_X.Value"=xValues, 
				"PC_Y.Value"=yValues, "Color"=col, 
				row.names=1:ncol(m), stringsAsFactors=F)
	return( out)
}

