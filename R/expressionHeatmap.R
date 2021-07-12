# expressionHeatmap.R -- turn abundance into a log 2 fold change heatmap


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

