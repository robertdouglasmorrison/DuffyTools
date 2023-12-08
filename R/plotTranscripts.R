# plotTranscripts.R  -- assorted plots for comparison of transcripts and ratios


`plotRatios` <- function( file, geneColumn="GENE_ID", value1Column="RPKM_1_M", value2Column="RPKM_2_M", cex=1,
			units="RPKM", offset=1, keepIntergenics=FALSE, label="Plot", plotType=c("Scatter","MA"),
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=1, 
			marker.pos=NULL, minYmax=NULL, maxYmax=NULL,
			sep="\t", sym.asp=TRUE, lab1="", lab2="", hideZero=FALSE, ...) {

	if ( is.character(file)) {
		tmp <- read.delim( file, as.is=T, sep=sep)
		cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tmp))
	} else if (is.data.frame(file)) {
		tmp <- file
	} else {
		stop( "Argument 'file' must be a character string or a data frame.")
	}

	if ( !( all( c( geneColumn, value1Column, value2Column) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, value1Column, value2Column, 
				"  Found: ", colnames(tmp))
		return()
	}

	# extract the parts we want
	genes <- tmp[[ geneColumn]]
	int1 <- as.numeric( tmp[[ value1Column]])
	int2 <- as.numeric( tmp[[ value2Column]])


	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		# gmap <- getCurrentGeneMap()
		# dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		# drops2 <- which( genes %in% dropableGenes)
		# drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		int1 <- int1[ -drops]
		int2 <- int2[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}

	plotType <- base::match.arg( plotType)
	if ( plotType == "MA") {
		ans <- makeMAplot( genes, int1, int2, offset=offset, units=units, label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, marker.pos=marker.pos,
			minYmax=minYmax, maxYmax=maxYmax, hideZero=hideZero, ...)
	} else {
		ans <- makeScatterplot( genes, int1, int2, offset=offset, units1=units, label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, marker.pos=marker.pos,
			minYmax=minYmax, maxYmax=maxYmax, sym.asp=sym.asp, lab1=lab1, lab2=lab2, hideZero=hideZero, ...)
	}

	return( invisible( ans))
}


`plotTranscripts` <- function( file1, file2, fid1=basename(file1), fid2=basename(file2), 
			gene1Column="GENE_ID", gene2Column=gene1Column, 
			value1Column="RPKM_M", value2Column=value1Column, cex=1,
			units1="RPKM", units2=units1, offset=1, keepIntergenics=FALSE, label="Plot", 
			plotType=c("Scatter","MA"),
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=1, 
			marker.pos=NULL, minYmax=NULL, maxYmax=NULL, showAllGenes=FALSE, 
			sep1="\t", sep2=sep1, sym.asp=TRUE, hideZero=FALSE, show.cor=TRUE, ...) {

	tmp <- read.delim( file1, as.is=T, sep=sep1)
	cat( "\nRead file: ", file1, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( gene1Column, value1Column) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", gene1Column, value1Column, "  Found: ", colnames(tmp))
		return()
	}
	genes1 <- tmp[[ gene1Column]]
	genes1 <- shortGeneName( genes1, keep=1)
	int1 <- as.numeric( tmp[[ value1Column]])

	tmp <- read.delim( file2, as.is=T, sep=sep2)
	cat( "\nRead file: ", file2, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( gene2Column, value2Column) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", gene2Column, value2Column, "  Found: ", colnames(tmp))
		return()
	}
	genes2 <- tmp[[ gene2Column]]
	genes2 <- shortGeneName( genes2, keep=1)
	int2 <- as.numeric( tmp[[ value2Column]])

	# combine and resolve
	if (showAllGenes) {
		both <- union( genes1, genes2)
		newint1 <- newint2 <- rep.int( 0, length(both))
		wh1 <- base::match( genes1, both)
		newint1[wh1] <- int1
		wh2 <- base::match( genes2, both)
		newint2[wh2] <- int2
		int1 <- newint1
		int2 <- newint2
		genes <- both
	} else {
		both <- intersect( genes1, genes2)
		wh1 <- base::match( both, genes1)
		int1 <- int1[ wh1]
		wh2 <- base::match( both, genes2)
		int2 <- int2[ wh2]
		genes <- both
	}

	# put into MA order
	v <- log2( (int1+1)/(int2+1))
	ord <- base::order( v, decreasing=T)
	genes <- genes[ord]
	int1 <- int1[ord]
	int2 <- int2[ord]

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		# gmap <- getCurrentGeneMap()
		# dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		# drops2 <- which( genes %in% dropableGenes)
		# drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		int1 <- int1[ -drops]
		int2 <- int2[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}
	cat( "\nN_Genes in common: ",  length(genes), "\n")

	plotType <- match.arg( plotType)
	if ( plotType == "MA") {
		ans <- makeMAplot( genes, int1, int2, offset=offset, units=units1, label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, marker.pos=marker.pos,
			minYmax=minYmax, maxYmax=maxYmax, hideZero=hideZero, ...)
	} else {
		ans <- makeScatterplot( genes, int1, int2, offset=offset, units1=units1, units2=units2, 
			label=label, cex=cex,
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, marker.pos=marker.pos,
			minYmax=minYmax, maxYmax=maxYmax, sym.asp=sym.asp, lab1=fid1, lab2=fid2, 
			hideZero=hideZero, show.cor=show.cor, ...)
	}

	return( invisible( ans))
}


`plotRatiosTwoFiles` <- function( file1, file2, fid1=basename(file1), fid2=basename(file2), 
			geneColumn="GENE_ID", valueColumn="LOG2FOLD_M", cex=1,
			units="Log2 Fold", offset=0, keepIntergenics=FALSE, label="Plot", 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=1, 
			marker.pos=NULL, minYmax=NULL,maxYmax=NULL,  
			sep="\t", sym.asp=TRUE, hideZero=FALSE, show.cor=TRUE, ...) {

	tmp <- read.delim( file1, as.is=T, sep=sep)
	cat( "\nRead file: ", file1, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( geneColumn, valueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, valueColumn, "  Found: ", colnames(tmp))
		return()
	}
	# extract the parts we want
	genes1 <- tmp[[ geneColumn]]
	genes1 <- shortGeneName( genes1, keep=1)
	int1 <- as.numeric( tmp[[ valueColumn]])

	tmp <- read.delim( file2, as.is=T, sep=sep)
	cat( "\nRead file: ", file2, "\nN_Genes: ", nrow(tmp), "\n")
	if ( !( all( c( geneColumn, valueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, valueColumn, "  Found: ", colnames(tmp))
		return()
	}
	# extract the parts we want
	genes2 <- tmp[[ geneColumn]]
	genes2 <- shortGeneName( genes2, keep=1)
	int2 <- as.numeric( tmp[[ valueColumn]])

	# combine and resolve
	both <- intersect( genes1, genes2)
	wh1 <- base::match( both, genes1)
	genes1 <- genes1[ wh1]
	int1 <- int1[ wh1]
	wh2 <- base::match( both, genes2)
	genes2 <- genes2[ wh2]
	int2 <- int2[ wh2]
	genes <- both

	ord <- base::order( int1, decreasing=T)
	genes <- genes[ord]
	int1 <- int1[ord]
	int2 <- int2[ord]

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		# gmap <- getCurrentGeneMap()
		# dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		# drops2 <- which( genes %in% dropableGenes)
		# drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		int1 <- int1[ -drops]
		int2 <- int2[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}
	cat( "\nN_Genes in common: ",  length(genes), "\n")

	ans <- makeScatterplot( genes, int1, int2, offset=offset, units1=units, label=label, cex=cex,
		marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex,
		marker.labels=marker.labels, marker.pch=marker.pch, marker.pos=marker.pos,
		minYmax=minYmax, maxYmax=maxYmax, sym.asp=sym.asp, lab1=fid1, lab2=fid2, useLog="", 
		hideZero=hideZero, show.cor=show.cor, ...)
	
	return( invisible( ans))
}


`plotFoldChange` <- function( file, geneColumn="GENE_ID", foldColumn="LOG2FOLD", pvalueColumn="PVALUE", cex=1,
			keepIntergenics=FALSE, label="Plot", plotType=c("Volcano"),
			pch=21, col=c('blue','red', 'black'), cut.fold=1, cut.pvalue=0.05, shortNames=TRUE,
			signif.labels=TRUE, signif.label.cex=1, 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.label.cex=1, 
			marker.pch=21, marker.pos=NULL, sep="\t", min.intensity=0, intensityColumn="RPKM_1", 
			left.label=NULL, right.label=NULL, forceXmax=NULL, forceYmax=NULL, ...) {

	if ( is.character(file)) {
		tmp <- read.delim( file, as.is=T, sep=sep)
		cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tmp))
	} else if (is.data.frame(file)) {
		tmp <- file
	} else {
		stop( "Argument 'file' must be a character string or a data frame.")
	}

	if ( !( all( c( geneColumn, foldColumn, pvalueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, foldColumn, pvalueColumn, 
				"\n  \tFound: ", colnames(tmp))
		return()
	}

	# extract the parts we want
	genes <- tmp[[ geneColumn]]
	if (shortNames) genes <- shortGeneName( genes, keep=1)
	fold <- as.numeric( tmp[[ foldColumn]])
	pval <- as.numeric( tmp[[ pvalueColumn]])

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		# gmap <- getCurrentGeneMap()
		# dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		# drops2 <- which( genes %in% dropableGenes)
		# drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		fold <- fold[ -drops]
		pval <- pval[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}

	# allow the removal of very low expression genes
	if ( min.intensity > 0) {
		inten <- tmp[[ intensityColumn]]
		drops <- which( inten < min.intensity)
		if ( length(drops)) {
			genes <- genes[ -drops]
			fold <- fold[ -drops]
			pval <- pval[ -drops]
			cat( "\nAfter dropping low intensity: ", length(genes))
		}
	}


	plotType <- base::match.arg( plotType)
	if ( plotType == "Volcano") {
		ans <- makeVolcanoPlot( genes, fold, pval, label=label, cex=cex,
			pch=pch, col=col, cut.fold=cut.fold, cut.pvalue=cut.pvalue,
			signif.labels=signif.labels, signif.label.cex=signif.label.cex, 
			marker.genes=marker.genes, marker.col=marker.col, marker.cex=marker.cex, marker.label.cex=marker.label.cex,
			marker.labels=marker.labels, marker.pch=marker.pch, marker.pos=marker.pos, 
			left.label=left.label, right.label=right.label, forceXmax=forceXmax, forceYmax=forceYmax, ...)
	}

	return( invisible( ans))
}


`makeMAplot` <- function( genes, int1, int2, offset=1, units="Intensity", label="MA plot", 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.label.cex=1,
			marker.pch=1, marker.pos=NULL, minYmax=NULL, maxYmax=NULL, cex=1, hideZero=FALSE, ...) {

	if (hideZero) {
		zeros <- which( int1 == 0 | int2 == 0)
		if( length(zeros) > 0) {
			int1 <- int1[ -zeros]
			int2 <- int2[ -zeros]
			genes <- genes[ -zeros]
		}
	}

	# allow a linear shift of background intensity
	int1 <- int1 + offset
	int2 <- int2 + offset

	# A is average, M is fold change
	a <- log( sqrt(int1 * int2), 2)
	m <- log( (int1 / int2), 2)
	myYrange <- range(m) * 1.2
	if ( !is.null(minYmax)) {
		if (myYrange[2] < minYmax) myYrange[2] <- minYmax
		if (myYrange[1] > (-minYmax)) myYrange[1] <- (-minYmax)
	}
	if ( !is.null(maxYmax)) {
		if (myYrange[2] > maxYmax) myYrange[2] <- maxYmax
		if (myYrange[1] < (-maxYmax)) myYrange[1] <- (-maxYmax)
	}

	plot ( a, m, main=label, xlab=paste( "A:  log_2( Average",units,")"), 
		ylab=paste("M:  log_2( Ratio",units,")"), ylim=myYrange, cex=cex, font.axis=2, font.lab=2, ...)

#	addConfidenceLines( a, m, conf=0.95)

	# label selected genes, assumes the genes 'up in set1' are in the first half...
	if ( length( marker.genes) > 0) {
		if ( marker.genes[1] == "identify") {
			identify( a, m, shortGeneName(genes, keep=1), col=marker.col, cex=marker.cex)
		} else {
		who1 <- base::match( marker.genes, genes, nomatch=0)
		who2 <- base::match( alias2Gene(marker.genes), genes, nomatch=0)
		who <- pmax( who1, who2)
		marker.genes <- marker.genes[ who > 0]
		if ( length( marker.col) > 1) marker.col <- marker.col[ who > 0]
		who <- who[ who > 0]
		if ( length(who) > 0) {
			# put name above for first half, below for second half...
			if ( is.null( marker.pos)) {
				pos <- ifelse( m[who] > 0, 3, 1)
			} else { 
				pos <- marker.pos
			}
			points( a[who], m[who], col=marker.col, bg=marker.col, pch=marker.pch, cex=marker.cex)
			if (marker.labels) text( a[who], m[who], marker.genes, pos=pos, col=marker.col, cex=marker.label.cex)
		}}
	}
	dev.flush()

	Rpearson <- cor( int1, int2, use="complete", method="pearson")
	Rspearman <- cor( int1, int2, use="complete", method="spearman")

	return( list( "x"=a, "y"=m, "id"=genes, "Pearson_R"=Rpearson, "Spearman_Rho"=Rspearman))
}


`makeScatterplot` <- function( genes, int1, int2, offset=1, units1="Intensity", 
			units2=units1, label="MA plot", 
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.label.cex=1,
			marker.pch=1, marker.pos=NULL, minYmax=NULL, maxYmax=NULL, lab1="", lab2="", cex=1, useLog="xy", 
			hideZero=FALSE, sym.asp=TRUE, show.cor=TRUE, ...) {

	# we are plotting file 1 on Y axis, file 2 on X axis

	if (hideZero) {
		zeros <- which( int1 == 0 | int2 == 0)
		if( length(zeros) > 0) {
			int1 <- int1[ -zeros]
			int2 <- int2[ -zeros]
			genes <- genes[ -zeros]
		}
	}

	# allow a linear shift of background intensity
	int1 <- int1 + offset
	int2 <- int2 + offset

	myRange1 <- range( int1)
	myRange2 <- range( int2)
	if ( !is.null(minYmax)) {
		if (myRange1[2] < minYmax) myRange1[2] <- minYmax
		if (myRange2[2] < minYmax) myRange2[2] <- minYmax
	}
	if ( !is.null(maxYmax)) {
		if (myRange1[2] > maxYmax) myRange1[2] <- maxYmax
		if (myRange2[2] > maxYmax) myRange2[2] <- maxYmax
	}
	if (sym.asp) myRange1 <- myRange2 <- range( c( myRange1, myRange2))

	plot ( int2, int1, main=label, log=useLog, xlab=paste( lab2, "    (",units2, ")", sep=""), 
		ylab=paste( lab1, "    (", units1, ")", sep=""), xlim=myRange2, ylim=myRange1, 
		cex=cex, font.axis=2, font.lab=2, xaxt="n", yaxt="n", ...)

	xAts <- c( 0, axTicks(1)+offset)
	xAtLbls <- c( 0, axTicks(1))
	axis( 1, at=xAts, label=xAtLbls, font=2, ...)
	yAts <- c( 0, axTicks(2)+offset)
	yAtLbls <- c( 0, axTicks(2))
	axis( 2, at=yAts, label=yAtLbls, font=2, ...)

	Rpearson <- cor( log10(int2), log10(int1), use="complete", method="pearson")
	Rspearman <- cor( log10(int2), log10(int1), use="complete", method="spearman")
	if (show.cor) {
		cat( "\nPearson's R =", Rpearson, "\nN_Zeros: ", sum(int1 <= offset), sum(int2 <= offset))
		if (useLog == "xy") {
			try( abline( reg=lsfit( log10(int2), log10(int1)), col=2, lwd=3, untf=F))
		} else {
			abline( reg=lsfit( int2, int1), col=2, lwd=3)
		}
		legend( 'topleft', paste( "Pearson's R = ", formatC(Rpearson, format='f', digits=2), " "), bg='white')
	}
	legend( 'bottomright', paste( "N_Genes in common = ", length(genes)), bg='white')

	# label selected genes, assumes the genes 'up in set1' are in the first half...
	if ( length( marker.genes) > 0) {
		if ( marker.genes[1] == "identify") {
			identify( int2, int1, shortGeneName(genes, keep=1), col=marker.col, cex=marker.cex)
		} else {
		who1 <- base::match( marker.genes, genes, nomatch=0)
		who2 <- base::match( alias2Gene(marker.genes), genes, nomatch=0)
		who <- pmax( who1, who2)
		marker.genes <- marker.genes[ who > 0]
		if ( length( marker.col) > 1) marker.col <- marker.col[ who > 0]
		who <- who[ who > 0]
		if ( length(who) > 0) {
			# put name left for first half, right for second half...
			if ( is.null( marker.pos)) {
				pos <- ifelse( int1[who] > int2[who], 2, 4)
			} else {
				pos <- marker.pos
			}
			points( int2[who], int1[who], col=marker.col, bg=marker.col, pch=marker.pch, cex=marker.cex)
			if (marker.labels) text( int2[who], int1[who], marker.genes, pos=pos, col=marker.col, cex=marker.label.cex)
		}}
	}
	dev.flush()

	return( list( "x"=int2, "y"=int1, "id"=genes, "Pearson_R"=Rpearson, "Spearman_Rho"=Rspearman))
}


`makeVolcanoPlot` <- function( genes, fold, pval, label="Volcano Plot: ", 
			pch=21, col=c('blue','red','black'), cut.fold=1, cut.pvalue=0.05,
			clip.fold=10, clip.pvalue=1e-10, forceXmax=NULL, forceYmax=NULL,
			signif.labels=TRUE, signif.label.cex=1,
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.label.cex=1,
			marker.pch=21, marker.pos=NULL, cex=1, left.label=NULL, right.label=NULL, ...) {

	# we are plotting -log10(pval) on Y axis, Fold on X axis

	# prevent absurd numbers from dominating
	fold[ fold > clip.fold] <- clip.fold
	fold[ fold < -clip.fold] <- -clip.fold
	pval[ pval < clip.pvalue] <- clip.pvalue

	# add extra room on X for the labels
	# put some minimums in, so no change is clear
	bigX <- max( 1, abs( fold), na.rm=T)
	myRangeX <- range( c( -1, 1, fold), na.rm=T) * 1.15
	if ( ! is.null( forceXmax)) {
		xMax <- as.numeric(forceXmax)
		myRangeX[1] <- -xMax
		myRangeX[2] <- xMax
		fold[ fold > xMax] <- xMax
		fold[ fold < -xMax] <- -xMax
	}

	y <- -( log10( pval))
	myRangeY <- c( 0, max( 1, y, na.rm=F)*1.05)
	if ( ! is.null( forceYmax)) {
		yMax <- as.numeric(forceYmax)
		myRangeY[2] <- yMax
		y[ y > yMax] <- yMax
	}

	# use color to show when we exceed the wanted fold change
	N <- length( fold)
	myCol <- rep.int( col[1], N)

	# change the logic a bit, to separate the significance from the marker genes
	# first color by P-value
	isSignif <- which( abs(fold) >= cut.fold & pval <= cut.pvalue)
	myCol[ isSignif] <- col[2]

	if ( ! is.null( marker.genes)) {
		# given explicit genes, just show those
		whereMarker <- match( marker.genes, genes, nomatch=0)
		if ( all( whereMarker == 0)) whereMarker <- pmatch( marker.genes, genes, nomatch=0)
		myCol[ whereMarker] <- marker.col
	}

	plot ( fold, y, type="p", main=label, xlab="Log2 Fold Change",
		ylab="-Log10 P", xlim=myRangeX, ylim=myRangeY, 
		pch=pch, col=myCol, bg=myCol, cex=cex, font.axis=2, font.lab=2, ...)

	# label selected genes, assumes the genes 'up in set1' are in the first half...
	# add labels when the P-value is great enough
	# now do the significance and the markers as 2 separate passes
	if (signif.labels && length(isSignif)) {
		myPos <- ifelse( fold[isSignif] > 0, 4, 2)
		if ( length(isSignif) > 2) {
			require( plotrix)
			myPos <- thigmophobe( fold[isSignif], y[isSignif])
			# but force these to take right sides
			myPos[ fold[isSignif] < 0 & myPos == 4] <- 2
			myPos[ fold[isSignif] > 0 & myPos == 2] <- 4
		}
		text( fold[isSignif], y[isSignif], genes[isSignif], pos=myPos, cex=signif.label.cex)
	}

	if ( length( marker.genes) > 0) {
		if ( marker.genes[1] == "identify") {
			identify( fold, y, shortGeneName(genes, keep=1), col=marker.col, cex=marker.cex)
		} else {
			who1 <- match( marker.genes, genes, nomatch=0)
			who2 <- match( alias2Gene( marker.genes), genes, nomatch=0)
			who <- pmax( who1, who2)
			marker.genes <- marker.genes[ who > 0]
			if ( length( marker.col) > 1) marker.col <- marker.col[ who > 0]
			who <- who[ who > 0]
			if ( length(who) > 0) {
				# put name left for down reg, right for up
				if ( is.null( marker.pos)) {
					pos <- ifelse( fold[who] < 0, 2, 4)
				} else {
					pos <- marker.pos
				}
				points( fold[who], y[who], col=marker.col, bg=marker.col, pch=marker.pch, cex=marker.cex)
				if (marker.labels) text( fold[who], y[who], marker.genes, pos=pos, col=col[3], cex=marker.label.cex)
			}
		}
	}

	# optional labels to remind which group is which

	#if ( !is.null(left.label)) text( myRangeX[1]*.75, myRangeY[2]*0.025, paste( "UP in '", left.label, "'", sep=""), cex=1, font=2)
	if ( !is.null(right.label)) text( myRangeX[1]*.75, myRangeY[2]*0.025, paste( "DOWN in '", right.label, "'", sep=""), cex=1, font=2)
	if ( !is.null(right.label)) text( myRangeX[2]*.75, myRangeY[2]*0.025, paste( "UP in '", right.label, "'", sep=""), cex=1, font=2)
	dev.flush()

	return( list( "x"=fold, "y"=y, "id"=genes ))
}


addConfidenceLines <- function( int1, int2, conf=0.95) {

	# get the regression line model, and its fitted values with std errs.
	model <- lm( int2 ~ int1)
	preds <- predict( model, se.fit=TRUE)

	# 'Working-Hotelling multiplier'
	W <- sqrt( 2 * qf( conf, 2, (length(int1) - 2)))

	# upper and lower conf bands
	dw <- W * preds$se.fit
	bands <- data.frame( "lowBand"=preds$fit - dw, "highBand"=preds$fit + dw)

	# the regression line
	abline( model$coeficents[1], model$coeficents[2], untf=TRUE, col=2, lty=1)

	ordX <- base::order(int1)
	lines( int1[ordX], bands$lowBand[ordX], lty=2, col=2)
	lines( int1[ordX], bands$highBand[ordX], lty=2, col=2)
	return()
}
