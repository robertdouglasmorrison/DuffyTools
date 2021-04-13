# expressionMatrixToDEfiles.R -- convert one matrix of expression data into a folder of DE files
#				 where each group is compared against the mean of all groups

`expressionMatrixToDEfiles` <- function( x, groups=colnames(x), folder=".", offset=1.0, AVG.FUN=sqrtmean,
					sep="\t", units="RKPM") {

	gids <- rownames(x)
	if ( is.null(gids) || all( gids == 1:nrow(x))) stop( "expression matrix must have gene rownames")
	NR <- nrow(x)
	groupIDs <- colnames(x)
	NC <- ncol(x)
	prefix <- getCurrentSpeciesFilePrefix()
	suffix <- if ( sep == ",") "csv" else "txt"

	if ( ! file.exists( folder)) dir.create( folder, recursive=T)

	# see if we reduce replicates by group
	if ( ! all( colnames(x) == groups)) {
		groupFac <- factor( groups)
		groupIDs <- levels( groupFac)
		NC <- nlevels(groupFac)
		cat( "\nReducing", ncol(x), "samples down to", NC, "groups..")
		x2 <- matrix( NA, NR, NC)
		colnames(x2) <- groupIDs
		rownames(x2) <- gids
		for ( i in 1:NR) x2[ i, ] <- tapply( x[i,], groupFac, FUN=AVG.FUN)
		x <- x2
	}

	# make one expression vector that is the average of all columns, to be the comparitor for each
	gAvg <- apply( x, MARGIN=1, AVG.FUN)

	# we will build a DE result file for each group, by doing all Fold Change calls first,
	# and then use those Mvalues to make the Pvalues
	cat( "\nCalculating Fold Change..")
	fcM <- matrix( NA, NR, NC)
	for( j in 1:NC) {
		myV <- x[ ,j]
		myFC <- log2( (myV+offset) / (gAvg+offset))
		fcM[ , j] <- myFC
	}
	cat( "\nCalculating P-values..")
	pvM <- matrix( NA, NR, NC)
	for ( i in 1:NR) {
		# wathc for values that break T test
		otherV <- fcM[ i, ]
		if ( all( abs( otherV) < 0.05)) next
		for( j in 1:NC) {
			myFC <- fcM[i,j]
			useV <- otherV[-j]
			if ( all( abs( useV) < 0.05)) next
			# small chance of not enough data for T test...
			pvM[ i, j] <- sparse.t.test( useV, mu=myFC)$p.value
		}
	}

	# make the results
	cat( "\nMaking result DE files..")
	prods <- gene2ProductAllSpecies(gids)
	celltypes <- gene2CellType(gids)
	for( j in 1:NC) {
		thisGrp <- groupIDs[j]
		thisRPKM <- x[ ,j]
		thisFC <- fcM[ ,j]
		thisPV <- pvM[ ,j]
		thisPV[ is.na(thisPV)] <- 1
		thisDF <- data.frame( "GENE_ID"=gids, "PRODUCT"=prods, "CellType"=celltypes, 
					"LOG2FOLD"=round(thisFC,digits=4), "PVALUE"=thisPV, 
					"RPKM"=round(thisRPKM,digits=2), "AVG_RPKM"=round(gAvg,digits=2), 
					stringsAsFactors=F)
		colnames(thisDF)[6:7] <- paste( c(thisGrp,"AVG"), units, sep="_")
		ord <- diffExpressRankOrder( thisFC, thisPV)
		thisDF <- thisDF[ ord, ]
		rownames(thisDF) <- 1:nrow(thisDF)
		outfile <- file.path( folder, paste( thisGrp, prefix, "Ratio", suffix, sep="."))
		write.table( thisDF, outfile, sep=sep, quote=(suffix == "csv"), row.names=F)
		cat( "\rWrote: ", j, thisGrp, outfile)
	}
	cat( "\nDone.\n")
}
