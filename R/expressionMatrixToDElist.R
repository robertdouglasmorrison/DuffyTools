# expressionMatrixToDElist.R -- convert one matrix of expression data into a list of DE tables
#				Like that used for meta rankings, etc.

`expressionMatrixToDElist` <- function( x, groups=colnames(x), minimumValue=NULL, AVG.FUN=mean) {

	deList <- vector( mode="list")

	groupIDset <- sort( unique( groups))

	if ( is.null( minimumValue)) {
		minimumValue <- min( x[ x > 0])
	}
	useX <- x + minimumValue

	# we will build a data frame for each group's set of columns
	for( i in 1:length( groupIDset)) {

		thisGroup <- groupIDset[i]
		myCols <- which( groups == thisGroup)
		otherCols <- which( groups != thisGroup)
		tmpG1 <- tmpG2 <- tmpV1 <- tmpV2 <- tmpF <- vector()

		for ( j in myCols) {
			tmpV1 <- c( tmpV1, useX[ ,j])
			tmpG1 <- c( tmpG1, rownames(x))
		}
		gFac <- factor( tmpG1)
		# if more than one column per group, average
		if ( length( myCols) > 1) {
			gFac <- factor( tmpG1)
			setV1 <- tapply( tmpV1, gFac, c)
			tmpV1 <- tapply( tmpV1, gFac, FUN=AVG.FUN, na.rm=T)
			tmpG1 <- levels(gFac)
		} else {
			setV1 <- as.list( tmpV1[ order(tmpG1)])
		}

		for ( j in otherCols) {
			tmpV2 <- c( tmpV2, useX[ ,j])
			tmpG2 <- c( tmpG2, rownames(x))
		}
		# if more than one column per group, average
		if ( length( otherCols) > 1) {
			gFac <- factor( tmpG2)
			setV2 <- tapply( tmpV2, gFac, c)
			tmpV2 <- tapply( tmpV2, gFac, FUN=AVG.FUN, na.rm=T)
			tmpG2 <- levels(gFac)
		} else {
			setV2 <- as.list( tmpV2[ order(tmpG2)])
		}

		tmpF <- log2( tmpV1/tmpV2)

		# make P-values too
		tmpP <- rep.int( 1, nlevels(gFac))
		for ( k in 1:length( tmpP)) {
			tmpP[k] <- wilcox.test( setV1[[k]], setV2[[k]], exact=F)$p.value
		}

		# order like any other DE dataset
		ord <- diffExpressRankOrder( tmpF, tmpP)

		thisDF <- data.frame( "GENE_ID"=tmpG1[ord], "LOG2FOLD"=tmpF[ord], "PVALUE"=tmpP[ord], 
					"VALUE_1"=tmpV1[ord], "VALUE_2"=tmpV2[ord], stringsAsFactors=F)
		rownames(thisDF) <- 1:nrow(thisDF)

		deList[[i]] <- thisDF
		names( deList)[i] <- thisGroup
	}
	return( deList)
}
