# consensusRank.R

# find the overall ranking from multiple datasets


`consensusRank` <- function( fnames, fids=basename(fnames), colorSet=2:(length(fnames)+1), 
			geneColumn="GENE_ID", valueColumn="LOG2FOLD", sep="\t",
			speciesID=getCurrentSpecies(), keepNonGenes=FALSE, weights=1,
			altGeneMap=NULL, annotationColumns="PRODUCT", label="", cex=1.0, useLog=FALSE) {

	# pad all the needed args to the full length
	NF <- length(fnames)
	colorSet <- rep( colorSet, length.out=NF)
	geneColumn <- rep( geneColumn, length.out=NF)
	valueColumn <- rep( valueColumn, length.out=NF)
	weights <- round( rep( weights, length.out=NF))

	# make sure we can read them all
	# and get all the needed parts into commonly named sets
	dfList <- vector( mode="list", length=NF)
	for ( i in 1:NF) {
		thisF <- fnames[i]
		thisName <- geneColumn[i]
		thisValue <- valueColumn[i]

		if ( ! file.exists(thisF)) {
			cat( "\nFile not found: ", thisF, "   ...skipping...")
			next
		}
		tmp <- read.delim( thisF, as.is=T, sep=sep)
		if ( ! (thisName %in% colnames(tmp))) {
			cat( "\nGene name column not found: ", thisName, "   ...skipping...")
			next
		}
		if ( ! (thisValue %in% colnames(tmp))) {
			cat( "\nNumeric value column not found: ", thisValue, "   ...skipping...")
			next
		}

		mydf <- data.frame( "GENE_ID"=tmp[[ thisName]], "VALUE"=as.numeric( tmp[[ thisValue]]), 
				"RANK"=1:nrow(tmp), stringsAsFactors=TRUE)
		rownames(mydf) <- 1:nrow(tmp)
		dfList[[ i ]] <- mydf
	}
	names( dfList) <- basename( fnames)
	cat( "\nN_Files: ", length(dfList))
	bigNgenes <- 0
	for( i in 1:NF) {
		cat( "\n", names(dfList)[i], "\tN_Genes: ", nrow( dfList[[i]]))
		bigNgenes <- max( bigNgenes, nrow(dfList[[i]]))
	}
	if ( bigNgenes < 2) stop( "no Genes to Rank...")



	# trim out intergenics
	setCurrentSpecies( speciesID)
	gmap <- getCurrentGeneMap()
	nongenes <- vector()
	if ( ! keepNonGenes) {
		nongenes <- subset( gmap, REAL_G == FALSE)$GENE_ID
		for( i in 1:NF) {
			tmp <- dfList[[i]]
			if ( is.null(tmp)) next
			ngs <- which( tmp$GENE_ID %in% nongenes)
			if ( length(ngs) > 0) {
				tmp <- tmp[ -ngs, ]
				tmp$RANK <- 1:nrow(tmp)
				dfList[[ i ]] <- tmp
			}
		}
	}
	if ( ! is.null( altGeneMap)) {
		gmap <- altGeneMap
	}

	# find the common genes
	commonG <- gmap$GENE_ID
	for( i in 1:NF) {
		tmp <-dfList[[ i]]
		if ( is.null(tmp)) next
		commonG <- intersect( commonG, tmp$GENE_ID)
	}
	NG <- length( commonG)
	cat( "\nN_Genes in common to all files: ", NG)

	# trim all files down to those genes, keeping their order intact
	for( i in 1:NF) {
		tmp <-dfList[[ i]]
		if ( is.null(tmp)) next
		found <- ( tmp$GENE_ID %in% commonG)
		if ( sum(found) < nrow(tmp)) {
			tmp <- tmp[ found, ]
			tmp$RANK <- 1:nrow(tmp)
			dfList[[i]] <- tmp
		}
	}

	# we are now ready to record every gene's rank in every dataset, and summarize
	allRanks <- data.frame()
	for( i in 1:NF) {
		tmp <-dfList[[ i]]
		if ( is.null(tmp)) next

		cat( "\nSorting: ", names(dfList)[i], " by: ", valueColumn[i])
		ord <- base::order( tmp$VALUE, decreasing=T)
		tmp <- tmp[ ord, ]
		tmp$RANK <- 1:nrow(tmp)
	
		# allow different integer weights
		cat( "\n", names( dfList)[i], "\tWeight= ", weights[i])
		for( k in 1:weights[i]) allRanks <- rbind( allRanks, tmp)
	}
	gFactor <- factor( allRanks$GENE_ID)
	avgRank <- tapply( allRanks$RANK, INDEX=gFactor, FUN=sqrtmean)
	avgValue <- tapply( allRanks$VALUE, INDEX=gFactor, FUN=mean)

	out <- data.frame( "GENE_ID"=levels(gFactor), "VALUE"=avgValue, "RANK"=avgRank, stringsAsFactors=FALSE)
	ord <- base::order( out$RANK, -out$VALUE)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	if ( length( annotationColumns) > 0) {
		colptrs <- base::match( annotationColumns, colnames(gmap), nomatch=0)
		colptrs <- colptrs[ colptrs > 0]
		annotationColumns <- colnames(gmap)[ colptrs]
		
		where <- base::match( out$GENE_ID, gmap$GENE_ID)
		extras <- as.data.frame( gmap[ where, colptrs])
		colnames( extras) <- annotationColumns
		out <- cbind( out, extras)
	}

	# supplemental plot...
	yLim <- range( allRanks$VALUE)
	minY <- yLim[1]
	if (useLog) {
		minY <- min( allRanks$VALUE[ allRanks$VALUE > 0])
		yLim[1] <- minY
	}

	plot( 1,1, type="n", xlim=c(1,NG), ylim=yLim, xlab="Consensus Average Rank", ylab="Consensus Average Value",
		main=paste( "Consensus Rank Plot:", label, sep="\n"), font.axis=2, font.lab=2,
		log=if(useLog) "y" else "" )
	outGenes <- out$GENE_ID
	for( i in 1:NF) {
		tmp <- dfList[[i]]
		if ( is.null(tmp)) next
		xloc <- base::match( tmp$GENE_ID, outGenes)
		yval <- tmp$VALUE
		if (useLog) yval <- ifelse( yval > 0, yval, minY)
		points( xloc, yval, bg=colorSet[i], pch=21, col=1, cex=cex)
	}
	yval <- out$VALUE
	if (useLog) yval <- ifelse( yval > 0, yval, minY)
	points( x=1:nrow(out), y=yval, bg=1, pch=21, col=1, cex=cex)
	legend( "topright", c(fids,"Consensus"), pch=21, col=1, pt.bg=c(colorSet,1))

	colnames(out)[1:2] <- c( geneColumn[1], valueColumn[1])
	return( out)
}
