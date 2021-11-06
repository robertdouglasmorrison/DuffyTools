# expressionHeatmap.R -- turn abundance into a log 2 fold change heatmap

# trying to provide a few different ways to generate better heatmaps


# using the 'heatmap.plus' package...

`expressionHeatmap` <- function( m, Ngenes=NULL, average.FUN=median, minIntensity=0, small.offset=1, 
				heatColors=redgreen(75), ...) {


	require( gplots)
	require( heatmap.plus)

	mv <- expressionMatrixToMvalue( m, average.FUN=average.FUN, minIntensity=minIntensity, small.offset=small.offset)

	difs <- apply( mv, MARGIN=1, function(x) diff( range( x, na.rm=T)))
	difs[ is.na(difs)] <- 0

	drops <- which( difs <= 0)
	if ( length( drops) > 0) {
		cat( "\nDropping rows with zero variance: ", length( drops))
		mv <- mv[ -drops, ]
		difs <- difs[ -drops]
	}

	if ( ! is.null( Ngenes) && nrow(mv) > Ngenes) {
		cat( "\nTrimming to highest variance genes: ", Ngenes)
		ord <- order( difs, decreasing=T)
		who <- ord[ 1:Ngenes]
		mv <- mv[ sort(who), ]
		difs <- difs[ sort(who)]
	}

	# lastly, the colors will be divided linearly, but there may be outliers at the edges
	limits <- quantile( as.vector(mv), c(0.02, 0.98))
	# force these to be symmetric
	useLimit <- min( abs( limits))
	limits[1] <- -useLimit
	limits[2] <- useLimit
	cat( "\nClipping extreme M-values to limits: ", limits)
	mv[ mv < limits[1]] <- limits[1]
	mv[ mv > limits[2]] <- limits[2]

	# plot that heat
	if ( ! is.null(Ngenes) && Ngenes > 1000) cat( "\nCalling heatmap()..")
	ans <- heatmap( mv, heatColors=heatColors, ...)
	return( invisible(ans))
}


# using the 'RColorBrewer' and 'gplots' packages

`expressionHeatmap2` <- function( m, mode=c("fold.change","expression","log2.expression","rank"), average.FUN=median, 
				minIntensity=0, small.offset=1, crop.log2fold=NULL, min.difference=NULL,
				palette="YlOrRd", n.palette.colors=9, rev.palette=FALSE, n.heatmap.colors=50, 
				Rowv=TRUE, Colv=TRUE, dendrogram=c("both","row","column","none"), 
				scale=c("none","row","column"), trace=c("column","row","both","none"), margins=c(8,5), 
				key.xlab="", keysize=1.0, density.info=c("histogram","density","none"), key.title="Color Key",
				...) {

	require( gplots)
	require( RColorBrewer)

	checkX11( width=8, height=12)

	# we are always given gene expression data, perhaps transform it
	mode <- match.arg( mode)
	mUse <- m
	if ( mode == "fold.change") {
		mUse <- expressionMatrixToMvalue( m, average.FUN=average.FUN, minIntensity=minIntensity, small.offset=small.offset)
		if ( ! is.null( crop.log2fold)) {
			cropValue <- as.numeric( crop.log2fold)
			tooBig <- which( abs( mUse) > crop.log2fold)
			if ( length(tooBig)) {
				cat( "\nCropping", length(tooBig), "excessive Log2Fold values to: ", cropValue)
				mUse[ mUse > cropValue] <- cropValue
				mUse[ mUse < -cropValue] <- -cropValue
			}
		}
		if (key.xlab == "") key.xlab <- "Log2 Fold Change"
	} else if ( mode == "rank") {
		for ( i in 1:ncol(m)) mUse[ , i] <- rank( m[ , i])
		if (key.xlab == "") key.xlab <- "Expression Rank"
	} else if ( mode == "log2.expression") {
		for ( i in 1:ncol(m)) mUse[ , i] <- log2( m[ , i] + small.offset)
		mUse[ is.na( mUse)] <- 0
		mUse[ is.nan( mUse)] <- 0
		mUse[ is.infinite( mUse)] <- 0
		if (key.xlab == "") key.xlab <- "Log2 Expression"
	} else {
		# default is to leave the expression untransformed
		if (key.xlab == "") key.xlab <- "Expression Units"
	}
	
	# all expression transformations are done.  Now see if we should drop genes with too little difference
	if ( ! is.null( min.difference)) {
		minDifCut <- as.numeric( min.difference)
		mDiff <- apply( mUse, 1, function(x) diff( range( x, na.rm=T)))
		drops <- which( mDiff < minDifCut)
		if ( length( drops)) {
			cat( "\nDropping gene rows with too little difference.  Cutoff =", minDifCut)
			cat( "\nN_Rows_Dropped: ", length(drops))
			mUse <- mUse[ -drops, ]
			cat( "\nN_Rows_Kept:    ", nrow(mUse))
			if ( nrow( mUse) < 3) stop( "Too many gene rows removed by 'min.difference' setting...")
		}
	}
	
	# perhaps use color brewer methods/tools for heatmap coloring
	# 'palette' can be a palette name, or a vector of colors
	if ( is.character( palette) && length( palette) == 1 && nchar(palette) > 1) {
		palette <- brewer.pal( n.palette.colors, palette)
		if (rev.palette) palette <- rev( palette)
	} 
	if ( is.null(palette) || is.na( palette) || nchar(palette) < 2) {
		palette <- heatMapColors( n.palette.colors, "red-white-blue", plotRamp=F, rampExponent = 0.75)
	}
	heatmapcolors <- palette
	if ( length(palette) < n.heatmap.colors) {
		morecols <- colorRampPalette( palette)
		heatmapcolors <- morecols( n.heatmap.colors)
	} 
	
	# ready to make the heatmap
	dendrogram <- match.arg( dendrogram)
	scale <- match.arg( scale)
	trace <- match.arg( trace)
	density.info <- match.arg( density.info)
	heatAns <- heatmap.2( mUse, col=heatmapcolors, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, 
				scale=scale, trace=trace, margins=margins, key.xlab=key.xlab, density.info=density.info, 
				keysize=keysize, key.title=key.title, ...)

	return( invisible( heatAns))
}
