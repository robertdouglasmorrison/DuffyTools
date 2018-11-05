# colorTools.R  - add flexibility to R colors


rainbow <- function (n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1,
				cosineShiftMagnitude=0.03) {

	if ((n <- as.integer(n[1L])) < 1) return( character())
	if (start == end || any(c(start, end) < 0) || any(c(start, end) > 1)) 
            stop("'start' and 'end' must be distinct and in [0, 1].")

	# start with the standard intervals like builtin rainbow
	x <- seq.int(start, ifelse(start > end, 1, 0) + end, length.out = n) %% 1

	# use a slightly non-linear ramp, to get less red, green, dark blue in the final ramp
	if ( cosineShiftMagnitude != 0) {
		# those spots are at 0.00, 0.33, 0.67, so make a cosine adjustment in those intervals
		modX <- cos( 3 * pi * (x %% (1/3)))
		# any term that are exactly on the red/green/blue center points get no adjustment
		modX[ abs(modX) >= 1] <- 0
		# the visual effect is stronger at the blue end of the spectrum, so taper it a bit
		cosineFactor <- rep.int( cosineShiftMagnitude, length(x))
		isBluish <- which( x > 0.50)
		cosineFactor[ isBluish] <- cosineFactor[ isBluish] * (1-x[isBluish])
		modX <- modX * cosineFactor
		x <- (x + modX) %% 1
	}

	hsv(h=x, s, v, alpha)
}


heatMapColors <- function( nColors, inflexionPoint=0.5, plotRamp=TRUE) {

	# make a color map...  we are dividing R & G by transcriptome...
	redFraction <- inflexionPoint
	nbig <- round( nColors * redFraction)
	nsmall <- nColors - nbig
	zsmall <- rep( 0, times=nsmall)
	zbig <- rep( 0, times=nbig)
	rampup <- seq( 0.0001,1,length.out=nbig)^0.5
	rampdown <- seq( 1,0.0001,length.out=nsmall)^0.5
	colmap <- rgb( red=c( zsmall, rampup), green=c( rampdown, zbig), blue=0)
	colHgt <- c( -rampdown, rampup)

	# do the heatmap...
	if (plotRamp) plot( 1:nColors, colHgt, pch=21, cex=8, bg=colmap)

	return( colmap)
}


`adjustColor` <- function( col, adjust=0) {
	x <- col2rgb(col); 
	dx <- if ( adjust > 0) (255 - x) else x;
	newx <- x + (dx * adjust)
	newx <- ifelse( newx > 255, 255, newx);
	newx <- ifelse( newx < 0, 0, newx);
	return( rgb( t(newx), max=255))
}


`adjustColorSet` <- function( colors, max.adjust=0.5) {

	out <- colors
	if ( any( duplicated( colors))) {
		tapply( 1:length(colors), INDEX=factor(colors), FUN=function(x) {
				if ( length(x) < 2) return()
				n <- length(x)
				delta <- max.adjust / n
				for ( i in 2:n) {
					out[x[i]] <<- adjustColor( out[x[1]], adjust=(delta*(i-1)))
				}
				return()
			})
	}
	return( out)
}


`colorBySpecies` <- function( speciesSet, palette=c('cyan', 'springgreen', 'gold', 'hotpink'),
				intergenic=NULL, intergenic.color='brown') {

	# given a vector of species IDs, return some color choices for each
	spFac <- factor( speciesSet)
	nSP <- nlevels( spFac)
	
	if ( nSP > length( palette)) palette <- rep( palette, length.out=nSP)

	colorSet <- palette[ as.numeric( spFac)]
	legendColors <- palette[ 1:nSP]
	names(legendColors) <- levels(spFac)

	if ( ! is.null( intergenic)) {
		if ( is.logical( intergenic)) intergenic <- which( intergenic)
		colorSet[ intergenic] <- intergenic.color
		legendColors <- c( legendColors, "intergenic"=intergenic.color)
	}

	out <- list( 'colors'=colorSet, 'legend.colors'=legendColors)
	return(out)
}
