# localityPlot.R

# draw where the given list of genes are physically distributed on the genome

`localityPlot` <- function( geneSet, mode=c("separate", "combined"), seqMap=NULL, color=2, label="") {


	smap <- getCurrentSeqMap()
	gmap <- getCurrentGeneMap()

	if ( ! is.null( seqMap)) smap <- seqMap

	# plot the chromos...
	nChromo <- nrow( smap)

	# allow multiple sets at once...
	if ( typeof( geneSet) != "list") {
		geneSet <- list( "geneSet"=geneSet)
	}
	if ( length(color) != length(geneSet)) color <- rep( color, length.out=length( geneSet))

	if ( match.arg(mode) == "separate") {

	smallX <- 1
	bigX <- max( smap$LENGTH)
	padding <- 250000
	xlimits <- c( smallX-padding, bigX)
	ylimits <- c( 0, nChromo)

	plot( 1, 1, type="n", main=paste( "Distribution of Gene Locations", label, sep="\n"),
		xlab="Gene location on chromosome", ylab="Chromosome Number", 
		xlim=xlimits, ylim=ylimits)
	for ( i in 1:nChromo) {
		lines( x=c(1,smap$LENGTH[i]), y=c(i,i), type="l", lwd=3)
		text( x=0, y=i, label=smap$SEQ_ID[i], pos=2, cex=0.9)
	}

	for( j in 1:length( geneSet)) {
		oneSet <- geneSet[[j]]
		oneColor <- color[j]
		gmapPtr <- base::match( oneSet, gmap$GENE_ID, nomatch=0)
		gmapPtr <- gmapPtr[ gmapPtr > 0]
		chromo <- base::match( gmap$SEQ_ID[ gmapPtr], smap$SEQ_ID)
		gbeg <- gmap$POSITION[ gmapPtr]
		gend <- gmap$END[ gmapPtr]
		rect( gbeg, chromo-0.2, gend, chromo+0.2, col=oneColor, border=oneColor)
	}
	legend( "bottomright", names( geneSet), col=color, fill=color)
	}

	if ( match.arg(mode) == "combined") {

	smallX <- 1
	bigX <- 1000000
	padding <- 100000
	xlimits <- c( smallX-padding, bigX)
	yRange <- 0

	# do the density first, to get y scaling
	densityCurves <- gLocs <- vector( mode="list")
	for( j in 1:length( geneSet)) {
		oneSet <- geneSet[[j]]
		gmapPtr <- base::match( oneSet, gmap$GENE_ID, nomatch=0)
		gmapPtr <- gmapPtr[ gmapPtr > 0]
		chromo <- base::match( gmap$SEQ_ID[ gmapPtr], smap$SEQ_ID)
		gbeg <- gmap$POSITION[ gmapPtr]
		gend <- gmap$END[ gmapPtr]
		# convert to combined location...
		myXscale <- bigX / smap$LENGTH[chromo]
		gbeg <- gbeg * myXscale
		gend <- gend * myXscale
		densitySet <- base::unlist( base::mapply( FUN=`:`, gbeg, gend, USE.NAMES=FALSE))
		den <- density( densitySet, kernel="epanechnikov", n=2^10, cut=0)
		densityCurves[[j]] <- den
		gLocs[[j]] <- list( "beg"=gbeg, "end"=gend)
		yRange <- range( c( den$y, yRange))
	}
	lowY <- ( yRange[2] * -0.2)
	yhgt <- lowY * -0.1
	ylimits <- c( lowY, yRange[2]*1.15)

	plot( 1, 0, type="n", main=paste( "Distribution of Gene Locations (all chromosomes combined)", 
		label, sep="\n"), xlab="Gene location on 'combined chromosome'", ylab="Density", 
		xlim=xlimits, ylim=ylimits)
	lines( x=c(1,bigX), y=c(lowY,lowY), type="l", lwd=3)
	text( x=0, y=lowY, label="combined", pos=2, cex=0.9)

	# draw them all... and save to redraw from L to R
	xl <- yl <- xh <- yh <- cc <- vector()
	for( j in 1:length( geneSet)) {
		oneColor <- color[j]
		den <- densityCurves[[j]]
		thisLoc <- gLocs[[j]]
		rect( thisLoc$beg, lowY-yhgt, thisLoc$end, lowY+yhgt, col=oneColor, border=oneColor)
		lines( den$x, den$y, col=oneColor, lwd=2)
		xl <- c( xl, thisLoc$beg)
		yl <- c( yl, lowY-yhgt)
		xh <- c( xh, thisLoc$end)
		yh <- c( yh, lowY+yhgt)
		cc <- c( cc, oneColor)
	}
	ord <- order( xl)
	rect( xl[ord], yl[ord], xh[ord], yh[ord], col=cc[ord], border=cc[ord])

	legend( "topright", names( geneSet), col=color, fill=color)
	}

	return()
}
