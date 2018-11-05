
# heatMapColors.R - make a green to red color ramp

heatMapColors <- function( nColors, inflexionPoint=0.5, plotRamp=TRUE) {

	# make a color map...  we are dividing R & G by transcriptome...
	redFraction <- inflexionPoint
	nbig <- round( nColors * redFraction)
	nsmall <- nColors - nbig
	zsmall <- rep( 0, times=nsmall)
	zbig <- rep( 0, times=nbig)
	rampup <- seq( 0.02,1,length.out=nbig)^0.5
	rampdown <- seq( 1,0.02,length.out=nsmall)^0.5
	colmap <- rgb( red=c( zsmall, rampup), green=c( rampdown, zbig), blue=0)
	colHgt <- c( -rampdown, rampup)

	# do the heatmap...
	if (plotRamp) plot( 1:nColors, colHgt, pch=21, cex=8, bg=colmap)

	return( colmap)
}
