# DESeqTools.R

# DESeq -- differential expression analysis of digital gene expression data

#  	Simon Anders, Wolfgang Huber


DESeq.DiffExpress <- function( fnames, fids, m=NULL, groupSet, targetGroup=sort(groupSet)[1], geneColumn="GENE_ID", 
			intensityColumn="READS_M", keepIntergenics=FALSE, 
			minimumRPKM=1, missingGenes="fill", extraColumn=NULL,
			average.FUN=sqrtmean, wt.folds=1, wt.pvalues=1, fitType="local", 
			adjust.lowReadCounts=TRUE, ...) {

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

	# DESeq breaks if every gene has at least one zero count!!
	# try to catch and prevent...
	smallV <- apply( m, MARGIN=1, min)

	if ( sum( smallV == 0) > nrow(m)*0.9) {
		cat( "\nCatch and compensate for 'Every gene contains at least one zero' DESeq bug...")
		# find who is not zero
		bigV <- apply( m, MARGIN=1, max)
		twoPlus <- which( bigV > 1)
		avgCnt <- mean( bigV[twoPlus])
		rowTweak <- which( bigV >= avgCnt)
		for( j in rowTweak) {
			# these are the biggest intensity rows, change the zeros to one
			# this should best preserve the biggest differences in expression...
			whoZ <- (m[ j, ] == 0)
			m[ j, whoZ] <- 1
		}
		cat( "\nApplied a 1 to ", length(rowTweak), " rows with minimum read count of ", avgCnt)
	}

	# Fold Change is based on read counts, not RPKM, so make a proxy for minimum read depth
	minimumREADS <- minimumRPKM * 10

	# there are now 2 versions of DESeq -- see which we have
	current <- installed.packages()
	packages <- current[ ,"Package"]
	hasDESeq1 <- hasDESeq2 <- FALSE
	if ( "DESeq2" %in% packages) hasDESeq2 <- TRUE
	if ( "DESeq" %in% packages) hasDESeq1 <- TRUE

	# and see which is more appropriate -- the newer DEseq2 is terrible when no replicates
	# make this call based on the 'target' vs non-target group counts, not just on the groupIDs passed in
	mode <- ""
	grpLabel <- rep( 'NonTarget', times=length(groupSet))
	grpLabel[ groupSet == targetGroup] <- targetGroup
	grpCounts <- table( grpLabel)
	nSolos <- sum( grpCounts == 1)
	if ( hasDESeq1 && all( grpCounts == 1)) mode <- "DESeq"
	if ( hasDESeq1 && nSolos > 0) mode <- "DESeq"
	if ( hasDESeq2 && mode == "") mode <- "DESeq2"
	if ( hasDESeq1 && mode == "") mode <- "DESeq"
	if ( mode == "") stop( "No version of DESeq installed!")
	if ( mode == "DESeq2") {
		cat( "\nUsing newer 'DESeq2'\n")
	} else {
		cat( "\nUsing older 'DESeq'\n")
	}
	require( mode, character.only=TRUE)

	if ( mode == "DESeq") {
		# the column flags for DESeq are to end up as 1,2, where 1 is the one we want...
		mine <- 1
		other <- 2
		cl <- rep( other, times=length(groupSet))
		cl[ groupSet == targetGroup] <- mine
		myGrpNames <- rep( targetGroup, times=length(groupSet))
		myGrpNames[ cl == other] <- notTargetGroup <- paste( "Not", targetGroup, sep=" ")
		if ( ! is.null( extraColumn)) {
			clExtra <- cl
		}
	
		# the dispersion has to be estimated differently if there aren't replicates
		specialDispersion <- FALSE
		if ( any( table(cl) < 2)) specialDispersion <- TRUE
	
		mUse <- matrix( as.integer(m), nrow=nrow(m), ncol=ncol(m))
		colnames(mUse) <- colnames(m)
		rownames(mUse) <- rownames(m)
	
		NG <- nrow( m)
	
		# call DESeq,  it's a multi step process
		ans <- newCountDataSet( countData=mUse, conditions=cl)
		ans <- estimateSizeFactors( ans)
		if (specialDispersion) {
			ans <- estimateDispersions( ans, method="blind", sharingMode="fit-only", fitType=fitType)
		} else {
			ans <- estimateDispersions( ans, fitType=fitType)
		}
		deseqOut <- nbinomTest( ans, other, mine)
	
		# DESeq throws away rows that were all zero, and may do divide by zero
		gnames <- rownames(mUse)
		deseqGenes <- deseqOut[ ,1]
		gPtr <- match( gnames, deseqGenes, nomatch=0)
		gprod <- gene2ProductAllSpecies( gnames)
		if ( any( gprod == "")) {
			tmp <- read.delim( fnames[1], as.is=T)
			gnams <- tmp[[ geneColumn]]
			gpros <- tmp[[ grep( "product", tolower(colnames(tmp)))]]
			needs <- which( gprod == "")
			where <- match( gnames[needs], gnams, nomatch=0)
			gprod[ needs[ where > 0]] <- gpros[ where]
		}
	
		v1 <- v2 <- foldOut <- rep.int( 0, length(gnames))
		v1[ gPtr > 0] <- deseqOut$baseMeanA[ gPtr]
		v2[ gPtr > 0] <- deseqOut$baseMeanB[ gPtr]
		# use our floor of minimum reads to prevent divide by zero
		foldOut <- log2( (v2+minimumREADS) / (v1+minimumREADS))
		pvalOut <- rep.int( 1, length(gnames))
		pvalOut[ gPtr > 0] <- deseqOut$pval[ gPtr]
		pvalOut[ is.na( pvalOut)] <- 1
		pivalOut <- piValue( foldOut, pvalOut)
	
		# round to sensible digits of resolution
		foldOut <- round( foldOut, digits=4)
		pivalOut <- round( pivalOut, digits=4)
		v1 <- round( v1, digits=4)
		v2 <- round( v2, digits=4)

		out <- data.frame( gnames, gprod, foldOut, pvalOut, pivalOut, v2, v1,
				stringsAsFactors=F)
		colnames(out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", "PIVALUE",
				targetGroup, notTargetGroup)

	}  # old DESeq package
	
	if ( mode == "DESeq2") {
		# the column flags for DESeq2 are to end up as 1,2, where 2 is the one we want...
		mine <- 2
		other <- 1
		cl <- rep( other, times=length(groupSet))
		cl[ groupSet == targetGroup] <- mine
		myGrpNames <- rep( targetGroup, times=length(groupSet))
		myGrpNames[ cl == other] <- notTargetGroup <- paste( "Not", targetGroup, sep=" ")
		if ( ! is.null( extraColumn)) {
			clExtra <- cl
		}
		colData <- data.frame( "group"=as.character(cl), "groupName"=myGrpNames)
	
		# counts as integers
		mUse <- matrix( as.integer(m), nrow=nrow(m), ncol=ncol(m))
		colnames(mUse) <- colnames(m)
		rownames(mUse) <- rownames(m)

		dds <- DESeqDataSetFromMatrix( countData=mUse, colData=colData, design= ~ group)
		ddsAns <- DESeq( dds, quiet=TRUE, fitType=fitType)
		res <- results( ddsAns, independentFiltering=FALSE, cooksCutoff=FALSE)
		resDF <- as.data.frame( res)

		# extract what we want
		gnames <- rownames(resDF)
		gprod <- gene2ProductAllSpecies( gnames)
		foldOut <- resDF$log2FoldChange
		pvalOut <- resDF$pvalue

		# calc the average per group, using the normalized count data
		mNorm <- counts( ddsAns, normalized=TRUE)
		v2 <- apply( mNorm[ , which( cl == mine), drop=FALSE], MARGIN=1, FUN=average.FUN)
		v1 <- apply( mNorm[ , which( cl == other), drop=FALSE], MARGIN=1, FUN=average.FUN)

		# the fold change reported when read count is near zero is unrealistic
		# use a graduated correction when read counts go too low
		if (adjust.lowReadCounts) {
			adjustAns <- lowReadCountAdjustment( foldOut, pvalOut, v2, v1)
			foldOut <- adjustAns$fold.change
			pvalOut <- adjustAns$p.value
		}

		# now we can assess the PI values
		pivalOut <- piValue( foldOut, pvalOut)

		# round to sensible digits of resolution
		foldOut <- round( foldOut, digits=4)
		pivalOut <- round( pivalOut, digits=4)
		v1 <- round( v1, digits=2)
		v2 <- round( v2, digits=2)

		out <- data.frame( gnames, gprod, foldOut, pvalOut, pivalOut, v2, v1,
				stringsAsFactors=F)
		colnames(out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", "PIVALUE",
				targetGroup, notTargetGroup)
	}

	if ( ! is.null( extraColumn)) {
		avgExtra1 <- apply( mExtra[ , which( cl == mine)], MARGIN=1, FUN=average.FUN)
		avgExtra2 <- apply( mExtra[ , which( cl == other)], MARGIN=1, FUN=average.FUN)
		out$AVG_EXTRA1 <- avgExtra1
		out$AVG_EXTRA2 <- avgExtra2
	}

	# final order by both fold and P-value
	ord <- diffExpressRankOrder( out$LOG2FOLD, out$PVALUE, wt.folds, wt.pvalues)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}


`lowReadCountAdjustment` <- function( fold, pval, reads1, reads2, minimumRPKM=1) {

	# the fold change reported by tools like DESeq and EdgeR when read counts are near zero 
	# are unrealistic.  Often genes that are below the threshold of what we consider expressed
	# report extremely high fold change and extremely low P-values.  Mostly since the variance of
	# a bunch of zeros is artificially low.

	# try to compensate

	# step 1a:  we always add a small linear offset in the log2 fold calculation to prevent divide by zero
	minimumREADS <- minimumRPKM * 10
	conservativeFold <- log2( (reads1+minimumREADS) / (reads2+minimumREADS))
	
	# step1b: whichever fold change is more conservative, use that one
	foldOut <- ifelse( abs(fold) < abs(conservativeFold), fold, conservativeFold)
	
	# Step 2a: see how much the fold change got reduced, and adjust the p-value accordingly.
	foldRatio <- abs( foldOut / fold)
	foldRatio[ foldRatio > 1.0] <- 1.0

	# step 2b:  apply adjustment that backs off the P-value smoothly
	pvalOut <- pval ^ sqrt(foldRatio)

	return( list( "fold.change"=foldOut, "p.value"=pvalOut))
}

