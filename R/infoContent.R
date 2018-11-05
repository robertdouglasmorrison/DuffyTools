# infoContent.R

`loadInfoContent` <- function( seqID="Pf3D7_01") {

	#setCurrentSpecies( speciesID)
	#prefix <- getCurrentSpeciesFilePrefix()
	#infoFile <- paste( prefix, "InfoContent", sep=".")
	infoFile <- paste( seqID, "InfoContent", sep=".")
	if ( exists( "InfoPath")) {
		infoFile <- file.path( InfoPath, paste( infoFile, "rda", sep="."))
		load( infoFile)
		infoContent <<- infoContent
	} else {
		data( list=infoFile)
	}
	return()
}


# add a info content curve to the current existing plot made by whoever...

`plotInfoContent` <- function( seqID, col=1, lwd=1, yMin=0, side=4, legendLocation=NULL) { 

	species <- getSpeciesFromSeqID( seqID)

	needLoad <- TRUE
	if ( exists( "infoContent") && infoContent$SeqID == seqID) needLoad <- FALSE

	if (needLoad) loadInfoContent( seqID)

	# get the info values for this entire chromosome
	v <- infoContent$Freq

	# use the current plot limits to see what part we want
	usr <- par("usr")
	xl <- round(usr[1]);   xh <- round(usr[2])
	yl <- round(usr[3]);   yh <- round(usr[4])
	# R stretches the plot limits by 4% so try to reverse that
	nX <- xh - xl + 1
	#xl <- xl + round(nX*0.05)
	#xh <- xh - round(nX*0.05)
	yh <- yh * 0.96
	nX <- xh - xl + 1
	myX <- xl:xh

	# the plot limits may extend beyond the chromosome ends...
	myY <- rep( NA, times=nX)
	fromXl <- max( 1, xl)
	toXh <- min( xh, length( v))
	dl <- fromXl - xl
	nX <- toXh - fromXl + 1
	myY[ (dl+1) : (dl+nX)] <- v[ fromXl:toXh]

	myYbig <- max( myY, 10, na.rm=T)
	myYbig <- log2( myYbig)
	myY <- log2( myY)
	isLog <- TRUE

	scaleFactor <- (yh) / myYbig
	myY <- (myY * scaleFactor) + yMin

	lines( myX, myY, col=col, lwd=lwd)

	# legend ?
	if ( ! is.null(legendLocation)) {
		myAt <- pretty( c( 1, myYbig), n=5)
		if (isLog) {
			myAt <- c( 0, myAt)
			myLabels <- formatC( 2^myAt, format="f", digits=0)
		} else {
			myLabels <- formatC( myAt, format="f", digits=myDigits)
		}
		axis( side=side, at=((myAt*scaleFactor)+yMin), labels=myLabels)
		mtext( "DNA Frequency  --  inverse( Information Content )", side=side, line=3)
		legend( legendLocation, c( "DNA Frequency", expression( '(Information)'^-1)), col=c(col,NA),
			lwd=lwd*2, bg="white")
	}
	return()
}


`getGeneInfoContent` <- function( geneID) {

	out <- rep( NA, times=length(geneID))
	curSeq <- ""

	for ( ig in 1:length(geneID)) {
		thisG <- geneID[ig]

		# get the exonic bases for this gene
		emap <- subset( getCurrentExonMap(), GENE_ID == thisG)
		if ( nrow( emap) < 1) emap <- subset( getCurrentGeneMap(), GENE_ID == thisG)
		if ( nrow( emap) < 1) {
			cat( "\nGene not found: ", thisG)
			next
		}

		# we may need to set up....
		seqID <- emap$SEQ_ID[1]
		if ( ig == 1) {
			species <- getSpeciesFromSeqID( seqID)
			needLoad <- TRUE
			if ( exists( "infoContent") && infoContent$Info$Species == species) needLoad <- FALSE
			if (needLoad) loadInfoContent( species)
		}

		if ( seqID != curSeq) {
			allInfo <- infoContent[[ seqID]]
			curSeq <- seqID
		}

		myV <- vector()
		for (ie in 1:nrow(emap)) {
			myX <- emap$POSITION[ie] : emap$END[ie]
			myV <- base::append( myV, allInfo[ myX])
		}
		out[ig] <- mean( myV, na.rm=TRUE)
	}
	names(out) <- geneID
	return( out)
}


makeInfoContentTable <- function( speciesID=getCurrentSpecies()) {

	setCurrentSpecies( speciesID)
	loadInfoContent( speciesID)

	gmap <- getCurrentGeneMap()
	infoScores <- getGeneInfoContent( gmap$GENE_ID)

	out <- data.frame( "GENE_ID"=gmap$GENE_ID, "PRODUCT"=gmap$PRODUCT, "INFO_SCORE"=infoScores,
			stringsAsFactors=FALSE)
	ord <- base::order( out$INFO_SCORE)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}


# add a SNP_Pvalue content curve to the current existing plot made by whoever...

`plotSNP_Pvalue` <- function( seqID, col=1, lwd=1, yMin=0, side=4, legendLocation=NULL) { 

	species <- getSpeciesFromSeqID( seqID)

	needLoad <- TRUE
	if ( exists( "infoContent") && infoContent$SeqID == seqID) needLoad <- FALSE

	if (needLoad) loadInfoContent( seqID)

	# get the info values for this entire chromosome
	if ( ! ( "SNP_Pvalue" %in% names( infoContent))) return()
	v <- infoContent$SNP_Pvalue

	# use the current plot limits to see what part we want
	usr <- par("usr")
	xl <- round(usr[1]);   xh <- round(usr[2])
	yl <- round(usr[3]);   yh <- round(usr[4])
	# R stretches the plot limits by 4% so try to reverse that
	nX <- xh - xl + 1
	#xl <- xl + round(nX*0.05)
	#xh <- xh - round(nX*0.05)
	yh <- yh * 0.86
	nX <- xh - xl + 1
	myX <- xl:xh

	# the plot limits may extend beyond the chromosome ends...
	myY <- rep( NA, times=nX)
	fromXl <- max( 1, xl)
	toXh <- min( xh, length( v))
	dl <- fromXl - xl
	nX <- toXh - fromXl + 1
	myY[ (dl+1) : (dl+nX)] <- v[ fromXl:toXh]

	# at this point turn Pvalues into values > 1
	myY[ myY < 1e-200] <- 1e-200
	myY <- -log10(myY)

	myYbig <- max( myY, 10, na.rm=T)
	isLog <- FALSE

	scaleFactor <- (yh) / myYbig
	myY <- (myY * scaleFactor) + yMin

	lines( myX, myY, col=col, lwd=lwd)

	# legend ?
	if ( ! is.null(legendLocation)) {
		myAt <- pretty( c( 1, myYbig), n=5)
		if (isLog) {
			myAt <- c( 0, myAt)
			myLabels <- formatC( 2^myAt, format="f", digits=0)
		} else {
			myLabels <- formatC( 10^(-myAt), format="e", digits=0)
		}
		axis( side=side, at=((myAt*scaleFactor)+yMin), labels=myLabels)
		mtext( "Genomic Base Probability", side=side, line=3)
		legend( legendLocation, c( "Genomic Base Probability", "(SNP_Pvalue)"), col=c(col,NA),
			lwd=lwd*2, bg="white")
	}
	return()
}