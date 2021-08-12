
# heatMapColors.R - make a green to red color ramp

heatMapColors <- function( nColors, palette=c("red-black-green", "red-white-blue", "red-white-green"),
				inflexionPoint=0.5, rampExponent=0.5, plotRamp=TRUE, pt.cex=30/sqrt(nColors)) {

	# make a color map...  for use by heatmaps, etc.
	palette <- match.arg( palette)

	if ( nColors < 3) nColors <- 5
	upFraction <- inflexionPoint
	nUp <- round( nColors * upFraction)
	nDown <- nColors - nUp
	smallStep <- max( round( 1/nColors, digits=2), 0.001)
	rampUp <- seq( smallStep, 1, length.out=nUp) ^ rampExponent
	rampUp[ nUp] <- 1.0
	rampDown <- seq( 1, smallStep, length.out=nDown) ^ rampExponent
	rampDown[ nDown] <- smallStep
	zeroDown <- rep.int( 0, nDown)
	zeroUp <- rep.int( 0, nUp)
	oneDown <- rep.int( 1, nDown)
	oneUp <- rep.int( 1, nUp)

	# different methods, based on the palette
	if ( palette == "red-black-green") {
		colmap <- rgb( red=c( zeroDown, rampUp), green=c( rampDown, zeroUp), blue=rep.int(0,nColors))
	} else if (palette == "red-white-blue") {
		colmap <- rgb( red=c( rev(rampDown), oneUp), green=c( rev(rampDown), rev(rampUp)), blue=c( oneDown, rev(rampUp)))
	} else if (palette == "red-white-green") {
		colmap <- rgb( red=c( rev(rampDown), oneUp), green=c( oneDown, rev(rampUp)), blue=c( rev(rampDown), rev(rampUp)))
	}
	colHgt <- c( -rampDown, rampUp)

	# do the heatmap...
	if (plotRamp) plot( 1:nColors, colHgt, pch=21, cex=pt.cex, bg=colmap)

	return( colmap)
}
