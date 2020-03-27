# samTools.R

# EdgeR -- differential expression analysis of digital gene expression data

#  	Mark Robinson, Davis McCarthy, Yunshun Chen, Gordon Smyth


EdgeR.DiffExpress <- function( fnames, fids, m=NULL, groupSet, targetGroup=sort(groupSet)[1], geneColumn="GENE_ID", 
			intensityColumn="READS_M", keepIntergenics=FALSE, missingGenes="fill", extraColumn=NULL,
			average.FUN=sqrtmean, minimumRPKM=1, wt.folds=1, wt.pvalues=1, ...) {

	# turn the set of transcript files into one matrix
	if ( is.null(m)) {
		m <- expressionFileSetToMatrix( fnames=fnames, fids=fids, geneColumn=geneColumn,
				intensityColumn=intensityColumn, missingGenes=missingGenes,
				keepIntergenics=keepIntergenics)
		if ( ! is.null( extraColumn)) {
			mExtra <- expressionFileSetToMatrix( fnames=fnames, fids=fids, geneColumn=geneColumn,
					intensityColumn=extraColumn, missingGenes=missingGenes,
					keepIntergenics=keepIntergenics)
		}
	}

	# drop non-genes...
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", rownames(m), fixed=T)
		if ( length(drops) > 0) {
			m <- m[ -drops, ]
			if ( ! is.null( extraColumn)) mExtra <- mExtra[ -drops, ]
		}
		nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
		drops <- which( rownames(m) %in% nonGenes)
		if ( length(drops) > 0) {
			m <- m[ -drops, ]
			if ( ! is.null( extraColumn)) mExtra <- mExtra[ -drops, ]
		}
	}

	# the column flags for EdgeR are to end up as 1,2, where 1 is the one we want...
	cl <- rep( 2, times=length(groupSet))
	cl[ groupSet == targetGroup] <- 1
	myGrpNames <- rep( targetGroup, times=length(groupSet))
	myGrpNames[ cl == 2] <- paste( "Not", targetGroup, sep=" ")
	if ( ! is.null( extraColumn)) {
		clExtra <- cl
	}

	# EdgeR will choke if less than two samples per group
	canDoDispersion <- TRUE
	dispersion <- "auto"
	if( sum( cl == 1) < 2) {
		canDoDispersion <- FALSE
	}
	if( sum( cl == 2) < 2) {
		canDoDispersion <- FALSE
	}

	grpFac <- factor( myGrpNames)
	grpNames <- levels(grpFac)

	mUse <- m

	require( edgeR)
	NG <- nrow( m)

	# call EdgeR,  it's a multi step process
	ans <- DGEList( counts=mUse, group=cl)
	ans <- calcNormFactors( ans)
	if ( canDoDispersion) {
		ans <- estimateCommonDisp( ans)
		ans <- estimateTagwiseDisp( ans)
	} else {
		dispersion <- 0.05
	}
	et <- exactTest( ans, pair=2:1, dispersion=dispersion)
	edgeRout <- topTags( et, n=NG)$table

	# get what e want/need out
	fout <- edgeRout[ , 1]
	pout <- edgeRout[ , 3]
	qout <- edgeRout[ , 4]
	gnames <- rownames( edgeRout)
	where <- match( gnames, rownames(m), nomatch=0)
	mTmp <- m[ where, ]

	# calc the average for each group, and put it into 'this group' order
	avgM <- t( apply( mTmp, MARGIN=1, function(x) tapply(x, grpFac, FUN=average.FUN)))
	if ( colnames(avgM)[1] != targetGroup) avgM <- avgM[ ,c(2,1)]

	# the fold change is based on read counts, not RPKM, so we need a different way to 
	# prevent divide by zero and falsely exagerated FC
	minimumREADS <- minimumRPKM * 10
	lowReads <- apply( avgM, 1, min)
	needRedo <- which( lowReads < (minimumREADS*5))
	if ( length( needRedo)) {
		newFold <- log2( (avgM[ needRedo, 1]+minimumREADS) / (avgM[ needRedo, 2]+minimumREADS))
		fout[ needRedo] <- newFold
	}

	# now we can assess PI values
	piout <- piValue( fout, pout)

	gprod <- gene2ProductAllSpecies( gnames)
	if ( any( gprod == "")) {
		tmp <- read.delim( fnames[1], as.is=T)
		gnams <- tmp[[ geneColumn]]
		gpros <- tmp[[ grep( "product", tolower(colnames(tmp)))]]
		needs <- which( gprod == "")
		where <- match( gnames[needs], gnams, nomatch=0)
		gprod[ needs[ where > 0]] <- gpros[ where]
	}

	# round to sensible digits of resolution
	fout <- round( fout, digits=4)
	avgM <- round( avgM, digits=2)
	piout <- round( piout, digits=3)

	out <- data.frame( gnames, gprod, fout, pout, qout, piout, avgM, stringsAsFactors=F)
	colnames(out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", "FDR", "PIVALUE", colnames(avgM))

	if ( ! is.null( extraColumn)) {
		avgExtra1 <- apply( mExtra[ , which( cl == 1)], MARGIN=1, FUN=average.FUN)
		avgExtra2 <- apply( mExtra[ , which( cl == 2)], MARGIN=1, FUN=average.FUN)
		out$AVG_EXTRA1 <- avgExtra1
		out$AVG_EXTRA2 <- avgExtra2
	}
	
	ord <- diffExpressRankOrder( out$LOG2FOLD, out$PVALUE, wt.folds, wt.pvalues)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}

