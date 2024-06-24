# vargeneTools.R -- work with the PF var gene domains


`getVargeneDomainMap` <- function() {

	vargeneDomainMap <- NULL
	data( "vargeneDomainMap", package="DuffyTools", envir=environment())
	if ( is.null( vargeneDomainMap)) return( data.frame())
	return( vargeneDomainMap)
}


`getVar2csaDomainMap` <- function( strain=NULL) {

	data( "var2csaDomainConstructs", envir=environment())
	myVar2csa <- var2csa
	if ( ! is.null( strain)) {
		myVar2csa <- subset( myVar2csa, STRAIN == strain)
	}
	return(myVar2csa)
}


`getVargeneDomains` <- function( genes) {

	# see if any of our domain sets cover this gene
	dmap3D7 <- subset.data.frame( getVargeneDomainMap(), GENE_NAME %in% genes)
	if ( nrow(dmap3D7)) {
		# make sure all possibly needed columns are the same, regardless of source
		dmap <- dmap3D7
		dmap$GENE_ID <- dmap$GENE_NAME
		dmap$AA_START <- dmap$AA_STOP <- NA
		dmap$DNA_START <- dmap$POSITION
		dmap$DNA_STOP <- dmap$END
		return( dmap)
	}
		
	dmapVSA <- subset.data.frame( getVSAdomainMap(), GENE_ID %in% genes)
	if ( nrow(dmapVSA)) {
		# make sure all possibly needed columns are the same, regardless of source
		dmap <- dmapVSA
		dmap$STRAND <- "+"
		return( dmap)
	}
		
	dmapJOS <- subset.data.frame( getJOSdomainMap(), GENE_ID %in% genes)
	if ( nrow(dmapJOS)) {
		# make sure all possibly needed columns are the same, regardless of source
		dmap <- dmapJOS
		dmap$STRAND <- "+"
		return( dmap)
	}
	return( data.frame())
}


`findVar2csaDomains` <- function( aaSeq, minScorePerAA=2, one.time=TRUE, plot.Y=NULL, rect.border=1, rect.fill='goldenrod',
				rect.halfHeight=1, rect.lwd=1, text.font=2, text.cex=1, text.X.repeat=NULL) {

	require( Biostrings)
	require( pwalign)
	data( "var2csaDomainConstructs", envir=environment())
	data( BLOSUM62)

	# use newer idea: find one best, mask, and repeat
	aaSeqNow <- aaSeq
	domDF <- data.frame()

	myVar2csa <- var2csa
	refLens <- nchar( myVar2csa$REF_SEQ)

	repeat {
		# find the best next one domain
		domAns1 <- pairwiseAlignment( myVar2csa$REF_SEQ, aaSeqNow, type="global-local", scoreOnly=T, 
					substitutionMatrix=BLOSUM62)
		domScores <- domAns1 / refLens
		best <- which.max( domScores)
		bestScore <- domScores[best]
		if ( bestScore < minScorePerAA) break
		# find its extent in the full length construct
		bestRefSeq <- myVar2csa$REF_SEQ[best]
		domAns2 <- pairwiseAlignment( bestRefSeq, aaSeqNow, type="global-local", scoreOnly=F, 
					substitutionMatrix=BLOSUM62)
		from <- start( subject( domAns2))
		to <- from + width( subject( domAns2)) - 1
		thisLen <- to - from + 1
		thisSubSeq <- substr( aaSeqNow, from, to)
		#cat( "\nDebug Init: ", from, to, thisSubSeq)

		# small chance of hitting mask residues on the flanks
		# with a very small chance of spanning into a second non-masked area too
		leftPatt <- "^[^X]{1,9}XXXXXXXXXX"
		if ( grepl( leftPatt, thisSubSeq)) substr( thisSubSeq, 1, 10) <- "XXXXXXXXXX"
		rightPatt <- "XXXXXXXXXX[^x]{1,9}$"
		if ( grepl( rightPatt, thisSubSeq)) substr( thisSubSeq, thisLen-9, thisLen) <- "XXXXXXXXXX"

		if ( substr( thisSubSeq, 1,1) == "X") {
			maskL <- sub( "(^X+)(.+)", "\\1", thisSubSeq)
			nMaskL <- nchar( maskL)
			from <- from + nMaskL
			thisSubSeq <- sub( maskL, "", thisSubSeq)
			thisLen <- nchar(thisSubSeq)
		}
		if ( substr( thisSubSeq, thisLen,thisLen) == "X") {
			maskR <- sub( "(^[^X]+)(X+$)", "\\2", thisSubSeq)
			nMaskR <- nchar( maskR)
			to <- to - nMaskR
		}
		thisLen <- to - from + 1
		if ( thisLen < 10) break
		thisSubSeq <- substr( aaSeqNow, from, to)
		#cat( "\nDebug Final: ", from, to, thisSubSeq)

		# if there still any 'XXX' in our hit, than what we grabbed completely spans a
		# domain we already masked out.   Call this an impossible outcome and quit
		if (grepl( "XXXXXX", thisSubSeq)) break

		# build the new small record and add it to the growing set
		smlDF <- myVar2csa[ best, ]
		smlDF$QUERY_START <- from
		smlDF$QUERY_STOP <- to
		smlDF$QUERY_DOMAIN_ID <- paste( smlDF$DOMAIN_ID, " (", smlDF$STRAIN, ")", sep="")
		smlDF$SCORE_PER_AA <- round( bestScore, digits=3)
		smlDF$EDIT_DIST <- adist( x=thisSubSeq, y=bestRefSeq)[1,1]
		smlDF$QUERY_SEQ <- thisSubSeq
		domDF <- rbind( domDF, smlDF)

		# remove that domain from the set of all domains
		if ( one.time) {
			thisDOMID <- smlDF$DOMAIN_ID
			myVar2csa <- myVar2csa[ -which( myVar2csa$DOMAIN_ID == thisDOMID), ]
			if ( ! nrow( myVar2csa)) break
			refLens <- nchar( myVar2csa$REF_SEQ)
		}

		# mask out that region so it can't get used again, and go find the next best
		maskStr <- paste( rep.int("X",thisLen),collapse="")
		substr( aaSeqNow, from, to) <- maskStr
	}
	if ( ! nrow(domDF)) return( domDF)

	# final order is by Starts in this query construct
	out <- domDF
	out <- out[ order( out$QUERY_START), ]
	rownames(out) <- 1:nrow(out)

	if ( is.null( plot.Y)) return( out)

	# draw them
	yLocs <- as.numeric( plot.Y)
	from <- out$QUERY_START
	to <- out$QUERY_STOP
	midpts <- (from + to) / 2
	txtstr <- out$QUERY_DOMAIN_ID
	for ( y in yLocs) {
		rect( from-0.3, y-rect.halfHeight, to+0.3, y+rect.halfHeight, border=rect.border, 
			col=rect.fill, lwd=rect.lwd)
		# labels go either at the mid-points, of repeatead along
		if ( is.null( text.X.repeat)) {
			text( midpts, y, txtstr, cex=text.cex, font=text.font)
		} else {
			# repeat the domain along the boxes
			for (k in 1:nrow(out)) {
				xAts <- seq( from[k]+12, to[k]-12, by=text.X.repeat[1])
				if ( length( xAts) < 2) xAts <- midpts[k]
				text( xAts, y, txtstr[k], cex=text.cex, font=text.font)
			}
		}
	}
	return( out)
}


`findVsaDomains` <- function( aaSeq, minScorePerAA=1, maskVSApattern=NULL, keepVSApattern=NULL, 
				plot.Y=NULL, rect.border=1, rect.fill='goldenrod', text.column="DOMAIN_ID",
				rect.halfHeight=1, rect.lwd=1, text.font=2, text.cex=1, text.X.repeat=NULL,
				verbose=FALSE) {

	require( Biostrings)
	require( pwalign)
	data( "VSA.DomainConstructs", envir=environment())
	data( BLOSUM62)

	# use newer idea: find one best, mask, and repeat
	aaSeqNow <- aaSeq
	domDF <- data.frame()

	myVsa <- vsa
	if ( ! is.null( maskVSApattern)) {
		toDrop <- vector()
		for ( mask in maskVSApattern) toDrop <- c( toDrop, grep( mask, myVsa$GENE_ID))
		toDrop <- sort( unique( toDrop))
		dropGenes <- sort( unique( myVsa$GENE_ID[toDrop]))
		if (verbose) {
			cat( "\nMasking VSA genes to exclude from search: ", length(dropGenes), "\n")
			print( head( dropGenes))
		}
		myVsa <- myVsa[ -toDrop, ]
	}
	if ( ! is.null( keepVSApattern)) {
		toKeep <- vector()
		for ( mask in keepVSApattern) toKeep <- c( toKeep, grep( mask, myVsa$GENE_ID))
		toKeep <- sort( unique( toKeep))
		keepGenes <- sort( unique( myVsa$GENE_ID[toKeep]))
		if (verbose) {
			cat( "\nLimiting VSA genes to include in search: ", length(keepGenes), "\n")
			print( head( keepGenes))
		}
		myVsa <- myVsa[ toKeep, ]
	}

	refLens <- nchar( myVsa$REF_SEQ)

	repeat {
		# find the best next one domain
		domAns1 <- pairwiseAlignment( myVsa$REF_SEQ, aaSeqNow, type="global-local", scoreOnly=T, 
					substitutionMatrix=BLOSUM62)
		domScores <- domAns1 / refLens
		best <- which.max( domScores)
		bestScore <- domScores[best]
		if ( bestScore < minScorePerAA) break
		# find its extent in the full length construct
		domAns2 <- pairwiseAlignment( myVsa$REF_SEQ[best], aaSeqNow, type="global-local", scoreOnly=F, 
					substitutionMatrix=BLOSUM62)
		from <- start( subject( domAns2))
		to <- from + width( subject( domAns2)) - 1
		thisLen <- to - from + 1
		thisSubSeq <- substr( aaSeqNow, from, to)
		#cat( "\nDebug Init: ", from, to, thisSubSeq)

		# small chance of hitting mask residues on the flanks
		# with a very small chance of spanning into a second non-masked area too
		leftPatt <- "^[^X]{1,9}XXXXXXXXXX"
		if ( grepl( leftPatt, thisSubSeq)) substr( thisSubSeq, 1, 10) <- "XXXXXXXXXX"
		rightPatt <- "XXXXXXXXXX[^x]{1,9}$"
		if ( grepl( rightPatt, thisSubSeq)) substr( thisSubSeq, thisLen-9, thisLen) <- "XXXXXXXXXX"
		if ( substr( thisSubSeq, 1,1) == "X") {
			maskL <- sub( "(^X+)(.+)", "\\1", thisSubSeq)
			nMaskL <- nchar( maskL)
			from <- from + nMaskL
			thisSubSeq <- sub( maskL, "", thisSubSeq)
			thisLen <- nchar(thisSubSeq)
		}
		if ( substr( thisSubSeq, thisLen,thisLen) == "X") {
			maskR <- sub( "(^[^X]+)(X+$)", "\\2", thisSubSeq)
			nMaskR <- nchar( maskR)
			to <- to - nMaskR
		}
		thisLen <- to - from + 1
		if ( thisLen < 10) break
		thisSubSeq <- substr( aaSeqNow, from, to)
		#cat( "\nDebug Final: ", from, to, thisSubSeq)

		# if there still any 'XXX' in our hit, than what we grabbed completely spans a
		# domain we already masked out.   Call this an impossible outcome and quit
		if (grepl( "XXXXXX", thisSubSeq)) break

		# build the new small record and add it to the growing set
		smlDF <- myVsa[ best, ]
		smlDF$QUERY_START <- from
		smlDF$QUERY_STOP <- to
		smlDF$QUERY_DOMAIN_ID <- paste( smlDF$DOMAIN_ID, " (", smlDF$STRAIN, ":", smlDF$GENE_ID, ")", sep="")
		smlDF$SCORE_PER_AA <- round( bestScore, digits=3)
		smlDF$QUERY_SEQ <- thisSubSeq
		smlDF$EDIT_DIST <- adist( thisSubSeq, myVsa$REF_SEQ[best])[1]
		# drop unwanted, and reorder?
		drops <- match( c( "GENE_ID", "STRAIN"), colnames(smlDF), nomatch=0)
		drops <- setdiff( drops, 0)
		if ( length( drops)) smlDF <- smlDF[ , -drops, drop=F]
		ord <- c( 1, 9, 10, 2, 3, 12, 4:8, 11)
		smlDF <- smlDF[ ,ord]
		# join this to the rest
		domDF <- rbind( domDF, smlDF)
		if (verbose) cat( "  Found", nrow(domDF), smlDF$QUERY_DOMAIN_ID[1], smlDF$SCORE_PER_AA[1])

		# mask out that region so it can't get used again, and go find the next best
		maskStr <- paste( rep.int("X",thisLen),collapse="")
		substr( aaSeqNow, from, to) <- maskStr
	}
	if ( ! nrow(domDF)) return( domDF)

	# final order is by Starts in this query construct
	out <- domDF
	out <- out[ order( out$QUERY_START), ]
	rownames(out) <- 1:nrow(out)

	if ( is.null( plot.Y)) return( out)

	# draw them
	yLocs <- as.numeric( plot.Y)
	from <- out$QUERY_START
	to <- out$QUERY_STOP
	midpts <- (from + to) / 2
	txtstr <- out[[ text.column]]
	for ( y in yLocs) {
		rect( from-0.3, y-rect.halfHeight, to+0.3, y+rect.halfHeight, border=rect.border, 
			col=rect.fill, lwd=rect.lwd)
		# labels go either at the mid-points, of repeatead along
		if ( is.null( text.X.repeat)) {
			text( midpts, y, txtstr, cex=text.cex, font=text.font)
		} else {
			# repeat the domain along the boxes
			for (k in 1:nrow(out)) {
				xAts <- seq( from[k]+12, to[k]-12, by=text.X.repeat[1])
				if ( length( xAts) < 2) xAts <- midpts[k]
				text( xAts, y, txtstr[k], cex=text.cex, font=text.font)
			}
		}
	}
	return( out)
}


`vargeneProfile` <- function( m, gene="PF3D7_0906000", name=gene, verbose=F, doPlot=T) {

	# we need the matrix of all gene expression
	NC <- ncol(m)
	NG <- nrow(m)
	if ( "GENE_ID" %in% colnames(m)) {
		mGenes <- m$GENE_ID 
		expressColumns <- 3:NC
	} else {
		mGenes <- rownames(m)
		expressColumns <- 1:NC
	}
	NS <- length( expressColumns)

	# we will profile the relative of all var genes, so gather all we need about them
	vgmap <- getVargeneDomainMap()
	vargenes <- sort( unique( vgmap$GENE_NAME))
	NVG <- length(vargenes)
	vargroups <- vgmap$GENE_GROUP[ match( vargenes, vgmap$GENE_NAME)]
	varlevels <- sort(unique(vargroups))
	NVGRPS <- length(varlevels)
	varcolors <- rainbow( NVGRPS, end=0.77)[ match( vargroups, varlevels)]
	vargene.express <- as.matrix( m[ match( vargenes, mGenes), expressColumns])
	rownames(vargene.express) <- vargenes
	vgRange <- range(vargene.express)

	# the gene of interest may be a geneID, or a vector of values
	if ( is.character( gene)) {
		who <- match( alias2Gene( gene[1]), mGenes, nomatch=0)
		if ( who == 0) {
			cat( "\nGiven gene not in matrix of expression:  ", gene[1])
			return( NULL)
		}
		gene.express <- as.numeric( m[ who, expressColumns])
		geneProd <- gene2Product( alias2Gene(gene))
		geneText <- name
		axisTextSuffix <- "Expression  (RPKM)"
	} else {
		# we were given a vector of values...
		gene.express <- gene
		if ( length(gene.express)  != NS) {
			cat( "\nUnequal size of data vector and expression matrix: ", length(gene), NS)
			return(NULL)
		}
		gene <- name
		geneProd <- ""
		geneText <- name
		axisTextSuffix <- "Value"
	}
	sids <- names(gene.express) <- colnames(m)[ expressColumns]
	sidSuffix <- substr( sids, nc <- nchar(sids), nc)
	canTestSM <- all( sidSuffix %in% c( "A", "L","M", "S"))

	if ( doPlot) {
		par( mai=c( 1,1,0.8,0.9))
		mainText <- paste( "How '",name, "' (", gene, ") Impacts Var Gene Expression Repertroire", sep="") 
		if ( geneProd != "") mainText <- paste( mainText, "\n", geneProd)
		plot( 1,1, type="n", main=mainText,
			xlab=paste( "Pragyan's MALI Isolates, in Order of Increasing Expression of ", name),
			ylab="Var Gene Expression  (RPKM)",
			xlim=c(1,NS+2.2), ylim=vgRange*c(0.5,2), log='y', 
			font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1)
	}

	# visit each sample in order of given values
	geneOrd <- order( gene.express)
	if (verbose) {
		cat( "\n", name, " Data:  ", gene.express[geneOrd])
		cat( "\n\nSampleIDs: ", sids[geneOrd])
	}

	bigX <- bigY <- bigGrp <- vector()
	for (i in 1:NS) {
		who <- geneOrd[i]
		myVG <- vargene.express[ ,who]
		if ( doPlot) {
			vgOrd <- order( myVG)
			#points( rep.int(i, NVG), myVG[vgOrd], pch=21, col=varcolors[vgOrd], bg=varcolors[vgOrd], cex=1.6)
			points( jitter(rep.int(i,NVG),amount=0.05), myVG, pch=21, col=varcolors, bg=varcolors, cex=1.4)
			if (showVAR2CSA <- TRUE) {
				isVar2csa <-match("PF3D7_1200600", rownames(vargene.express))
				points( i, myVG[ isVar2csa], pch=21, col=varcolors[isVar2csa], bg=varcolors[isVar2csa], cex=1.5)
			}
		}
		bigX <- c( bigX, rep.int( i, NVG))
		bigY <- c( bigY, myVG)
		bigGrp <- c( bigGrp, vargroups)
	}
	if (doPlot) text( 1:NS, apply(vargene.express,2,max)[geneOrd], sidSuffix[geneOrd], cex=1.5, pos=3, font=2)

	# can we make a Severe vs Mild call?
	if ( canTestSM) {
		smSet <- which( sidSuffix[geneOrd] %in% c("S","L"))
		mildSet <- which( sidSuffix[geneOrd] %in% c("A","M"))
		sm.pval <- wilcox.test( smSet, mildSet)$p.value
	}

	# show the gene level on the plot...
	if (doPlot) {
		geneRange <- range( gene.express)
		geneValues <- sort( gene.express)
		geneScaled <- ((geneValues - geneRange[1])/diff(geneRange)) * diff(vgRange) + vgRange[1]
		lines( 1:NS, geneScaled, lty=2, lwd=2, col=1)
		geneTicks <- c(1,2,5,10,20,50,100,200,500,1000,2000,5000)
		geneAt <- ((geneTicks - geneRange[1])/diff(geneRange)) * diff(vgRange) + vgRange[1]
		axis( side=4, at=geneAt, geneTicks, font=2, cex=1.1)
		mtext( paste( name, " ", axisTextSuffix), side=4, line=3, font=2, cex=1.1)
	}

	# measure all the var gene group trends w.r.t. this gene
	bigDF <- data.frame( "X"=bigX, "Y"=bigY, "Group"=bigGrp, stringsAsFactors=F)
	bigDF$Y <- log10( bigDF$Y)
	bigAns <- lm( Y ~ X, data=bigDF)
	ans <- by( bigDF, factor(bigGrp), function(mydf) {
				return( lm( Y ~ X, data=mydf))
			})

	grpSlope <- grpPval <- vector( length=NVGRPS)
	for ( i in 1:NVGRPS) {
		thisLM <- ans[[i]]
		deltaLMans <- lmSlopeDifference( bigAns, thisLM)
		grpSlope[i] <- deltaLMans$difference
		grpPval[i] <- deltaLMans$p.value
		if (doPlot) {
			abline( reg=thisLM, untf=F, col=varcolors[match(varlevels[i],vargroups)], lty=2, lwd=3)
		}
	}

	if (doPlot) {
		legend( 'topright', varlevels, pch=21, pt.bg=rainbow(NVGRPS,end=0.77), pt.cex=1.9, 
			bg='white', cex=1.2, title="Var Group")
		legend( 'bottomright', name, lwd=2, lty=2, bg='white', cex=1.2)
		if (canTestSM) legend( "top", paste( "Severe vs Mild:  P =", format(sm.pval, format="f", digits=3)),
			bg='white', cex=1.1)
	}

	out <- data.frame( "VarGroup"=varlevels, "DeltaSlope"=grpSlope, "P.value"=grpPval, stringsAsFactors=F)
	ord <- diffExpressRankOrder( abs(out$DeltaSlope), out$P.value, wt.fold=3)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)
	return(out)
}


`findVsaCassettes` <- function( aaSeq, minScorePerAA=2, maskVSAgenes=NULL, plot.Y=NULL, rect.border=1, rect.fill='goldenrod',
				rect.halfHeight=1, rect.lwd=1, text.font=2, text.cex=1, text.X.repeat=NULL) {

	require( Biostrings)
	require( pwalign)
	data( "VSA.CassetteConstructs", envir=environment())
	data( BLOSUM62)

	# use newer idea: find one best, mask, and repeat
	aaSeqNow <- aaSeq
	dcDF <- data.frame()

	myDC <- dc
	if ( ! is.null( maskVSAgenes)) {
		toDrop <- vector()
		for ( mask in maskVSAgenes) toDrop <- c( toDrop, grep( mask, myDC$GENE_ID))
		toDrop <- sort( unique( toDrop))
		dropGenes <- sort( unique( myDC$GENE_ID[toDrop]))
		cat( "\nMasking VSA genes to exclude from search: ", length(dropGenes), "\n")
		print( dropGenes)
		myDC <- myDC[ -toDrop, ]
	}

	refLens <- nchar( myDC$REF_SEQ)

	repeat {
		# find the best next one domain
		domAns1 <- pairwiseAlignment( myDC$REF_SEQ, aaSeqNow, type="global-local", scoreOnly=T, 
					substitutionMatrix=BLOSUM62)
		domScores <- domAns1 / refLens
		best <- which.max( domScores)
		bestScore <- domScores[best]
		if ( bestScore < minScorePerAA) break
		# find its extent in the full length construct
		domAns2 <- pairwiseAlignment( myDC$REF_SEQ[best], aaSeqNow, type="global-local", scoreOnly=F, 
					substitutionMatrix=BLOSUM62)
		from <- start( subject( domAns2))
		to <- from + width( subject( domAns2)) - 1
		thisLen <- to - from + 1
		thisSubSeq <- substr( aaSeqNow, from, to)

		# small chance of hitting mask residues on the flanks
		if ( substr( thisSubSeq, 1,1) == "X") {
			maskL <- sub( "(^X+)(.+)", "\\1", thisSubSeq)
			nMaskL <- nchar( maskL)
			from <- from + nMaskL
			thisSubSeq <- sub( maskL, "", thisSubSeq)
			thisLen <- nchar(thisSubSeq)
		}
		if ( substr( thisSubSeq, thisLen,thisLen) == "X") {
			maskR <- sub( "(^[^X]+)(X+$)", "\\2", thisSubSeq)
			nMaskR <- nchar( maskR)
			to <- to - nMaskR
		}
		thisLen <- to - from + 1
		thisSubSeq <- substr( aaSeqNow, from, to)

		# build the new small record and add it to the growing set
		smlDF <- myDC[ best, ]
		smlDF$QUERY_START <- from
		smlDF$QUERY_STOP <- to
		smlDF$QUERY_CASSETTE <- paste( smlDF$CASSETTE, " (", smlDF$STRAIN, ":", smlDF$GENE_ID, ")", sep="")
		smlDF$SCORE_PER_AA <- bestScore
		smlDF$QUERY_SEQ <- thisSubSeq
		smlDF$EDIT_DIST <- adist( thisSubSeq, myDC$REF_SEQ[best])[1]
		# drop unwanted, and reorder?
		drops <- match( c( "GENE_ID", "STRAIN"), colnames(smlDF), nomatch=0)
		drops <- setdiff( drops, 0)
		if ( length( drops)) smlDF <- smlDF[ , -drops, drop=F]
		#cat( "\nDebug: ", dim(smlDF), colnames(smlDF))
		ord <- c( 1, 9, 10, 2, 3, 12, 4:8, 11)
		smlDF <- smlDF[ ,ord]
		# join this to the rest
		dcDF <- rbind( dcDF, smlDF)

		# mask out that region so it can't get used again, and go find the next best
		maskStr <- paste( rep.int("X",thisLen),collapse="")
		substr( aaSeqNow, from, to) <- maskStr
	}
	if ( ! nrow(dcDF)) return( dcDF)

	# final order is by Starts in this query construct
	out <- dcDF
	out <- out[ order( out$QUERY_START), ]
	rownames(out) <- 1:nrow(out)

	if ( is.null( plot.Y)) return( out)

	# draw them
	yLocs <- as.numeric( plot.Y)
	from <- out$QUERY_START
	to <- out$QUERY_STOP
	midpts <- (from + to) / 2
	txtstr <- out$QUERY_CASSETTE
	for ( y in yLocs) {
		rect( from-0.3, y-rect.halfHeight, to+0.3, y+rect.halfHeight, border=rect.border, 
			col=rect.fill, lwd=rect.lwd)
		# labels go either at the mid-points, of repeatead along
		if ( is.null( text.X.repeat)) {
			text( midpts, y, txtstr, cex=text.cex, font=text.font)
		} else {
			# repeat the domain along the boxes
			for (k in 1:nrow(out)) {
				xAts <- seq( from[k]+12, to[k]-12, by=text.X.repeat[1])
				if ( length( xAts) < 2) xAts <- midpts[k]
				text( xAts, y, txtstr[k], cex=text.cex, font=text.font)
			}
		}
	}
	return( out)
}

