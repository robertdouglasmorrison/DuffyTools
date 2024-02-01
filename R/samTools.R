# samTools.R

# Statistical Analysis of Microarrays SAM


SAM.DiffExpress <- function( fnames, fids, m=NULL, 
			groupSet, targetGroup=sort(groupSet)[1], geneColumn="GENE_ID", 
			intensityColumn="INTENSITY", keepIntergenics=FALSE, 
			minimumIntensity=NULL, missingGenes="fill", extraColumn=NULL,
			average.FUN=sqrtmean,
			wt.folds=1, wt.pvalues=1, wt.dists=1, useLog=TRUE, needJitter=FALSE,
			...) {

	# turn the set of transcript files into one matrix
	if ( is.null( m)) {
		cat( "\nLoading transcriptomes..")
		m <- expressionFileSetToMatrix( fnames=fnames, fids=fids, geneColumn=geneColumn,
				intensityColumn=intensityColumn, missingGenes=missingGenes,
				keepIntergenics=keepIntergenics)
		if ( ! is.null( extraColumn)) {
			mExtra <- expressionFileSetToMatrix( fnames=fnames, fids=fids, geneColumn=geneColumn,
					intensityColumn=extraColumn, missingGenes=missingGenes,
					keepIntergenics=keepIntergenics)
		}
		cat( "  Done.\n")
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

	# the column flags for SAM are to end up as 0,1, where 1 is the one we want...
	# and the result column names are set from the input group names
	cl <- rep( 0, times=length(groupSet))
	cl[ groupSet == targetGroup] <- 1
	myGrpNames <- rep( targetGroup, times=length(groupSet))
	myGrpNames[ cl == 0] <- paste( "Not", targetGroup, sep=".")
	uniqueGroupNames <- unique( groupSet)
	otherGroups <- setdiff( uniqueGroupNames, targetGroup)
	if ( length(otherGroups) == 1) myGrpNames[ cl == 0] <- otherGroups
	
	if ( ! is.null( extraColumn)) {
		clExtra <- cl
	}

	# SAM will choke if less than two samples per group
	if( sum( cl == 1) < 3) {

		# so let's extend with more columns simulated around the one
		mycolumn <- which( cl == 1)[1]
		nNeed <- 3 - sum( cl == 1)
		nWas <- ncol(m)
		cat( "\nSimulating extra transcripts for under-sampled condition: ", colnames(m)[mycolumn])

		mnew <- matrix( 0, nrow=nrow(m), ncol=nWas+nNeed)
		rownames(mnew) <- rownames(m)
		colnames(mnew) <- c( colnames(m), rep( "", times=nNeed))
		for (k in 1:ncol(m)) mnew[ ,k] <- m[ ,k]
		mnew[ , (nWas+1):(nWas+nNeed)] <- t( sapply( m[ ,mycolumn], function(x) {
				if ( is.na(x) || x <= 0) return( rep( 0, times=nNeed))
				useSD <- if ( x > 1) sqrt(x)/3 else x/3
				return( abs( rnorm( nNeed, mean=x, sd=useSD)))
			}))
		cl <- c( cl, rep( 1, times=nNeed))
		myGrpNames <- c( myGrpNames, rep( targetGroup,  times=nNeed))
		colnames( mnew)[(nWas+1):(nWas+nNeed)] <- paste( colnames(m)[mycolumn], 1:nNeed, sep="_")
		m <- mnew
	}
	if( sum( cl == 0) < 3) {

		# so let's extend with 2 more columns simulated around the one
		mycolumn <- which( cl == 0)[1]
		nNeed <- 3 - sum( cl == 0)
		nWas <- ncol(m)
		cat( "\nSimulating extra transcripts for under-sampled condition: ", colnames(m)[mycolumn])

		mnew <- matrix( 0, nrow=nrow(m), ncol=nWas+nNeed)
		rownames(mnew) <- rownames(m)
		colnames(mnew) <- c( colnames(m), rep( "", times=nNeed))
		for (k in 1:ncol(m)) mnew[ ,k] <- m[ ,k]
		mnew[ , (nWas+1):(nWas+nNeed)] <- t( sapply( m[ ,mycolumn], function(x) {
				if ( is.na(x) || x <= 0) return( rep( 0, times=nNeed))
				useSD <- if ( x > 1) sqrt(x)/3 else x/3
				return( abs( rnorm( nNeed, mean=x, sd=useSD)))
			}))
		cl <- c( cl, rep( 0, times=nNeed))
		otherName <- paste( "Not", targetGroup, sep=" ")
		if ( length(otherGroups) == 1) otherName <- otherGroups
		myGrpNames <- c( myGrpNames, rep( otherName, times=nNeed))
		colnames( mnew)[(nWas+1):(nWas+nNeed)] <- paste( colnames(m)[mycolumn], 1:nNeed, sep="_")
		m <- mnew
	}

	grpFac <- factor( myGrpNames)
	grpNames <- levels(grpFac)

	# apply the linear offset s0 if given...
	if ( ! is.null( minimumIntensity)) m <- m + minimumIntensity

	# some things like proteomic data are very discrete, perhaps breaking SAM.  Jitter to fix
	if ( needJitter) {
		cat( "\nApplying 'jitter()' to intensities..")
		for( j in 1:ncol(m)) m[,j] <- jitter( m[,j])
		cat( "  Done.\n")
	}

	require( siggenes)

	# call SAM
	samM <- m
	if (useLog) samM <- log2(m)

	samOut <- sam( samM, cl, use.dm=FALSE, R.unlog=useLog, ...)

	# now extract and repackage the useful bits...
	distActual <- samOut@d
	distExpect <- samOut@d.bar
	falseCnt <- samOut@vec.false
	pval <- samOut@p.value
	qval <- samOut@q.value
	myfold <- samOut@fold
	myfold <- log2( myfold)
	gnames <- names( distActual)

	# catch bad data...
	missing <- which( is.na(distActual))
	distExpectOut <- rep( 0, times=length(distActual))
	distExpectOut[ ! is.na(distActual)] <- distExpect
	distActual[ missing] <- 0
	deltaDist <- distActual - distExpectOut

	# calc the average for each group, and put it into 'this group' order
	avgM <- t( apply( m, MARGIN=1, function(x) tapply(x, grpFac, FUN=average.FUN, na.rm=T)))
	if ( colnames(avgM)[1] != targetGroup) avgM <- avgM[ ,c(2,1)]
	# remove that linear offset from the averages...
	if ( ! is.null( minimumIntensity)) avgM <- avgM - minimumIntensity
	# and never let them go negative
	avgM <- pmax( avgM, 0)

	gprod <- gene2ProductAllSpecies( gnames)
	if ( all(gprod == "")) gprod <- gene2ProductAllSpecies( shortGeneName(gnames, keep=1))
	if ( any( gprod == "") && is.null(m)) {
		tmp <- read.delim( fnames[1], as.is=T)
		gnams <- tmp[[ geneColumn]]
		gpros <- tmp[[ grep( "product", tolower(colnames(tmp)))]]
		needs <- which( gprod == "")
		where <- match( gnames[needs], gnams, nomatch=0)
		gprod[ needs[ where > 0]] <- gpros[ where]
	}

	# round to sensible digits of resolution
	myfold <- round( myfold, digits=4)
	distActual <- round( distActual, digits=4)
	avgM <- round( avgM, digits=2)

	out <- data.frame( gnames, gprod, myfold, pval, qval, distActual, avgM, stringsAsFactors=F)
	colnames(out) <- c( "GENE_ID", "PRODUCT", "LOG2FOLD", "PVALUE", 
			"QVALUE", "DISTANCE", colnames(avgM))

	if ( ! is.null( extraColumn)) {
		avgExtra1 <- apply( mExtra[ , which( cl == 1)], MARGIN=1, FUN=average.FUN)
		avgExtra2 <- apply( mExtra[ , which( cl == 0)], MARGIN=1, FUN=average.FUN)
		out$AVG_EXTRA1 <- avgExtra1
		out$AVG_EXTRA2 <- avgExtra2
	}
	
	ord <- diffExpressDistanceRankOrder( out$LOG2FOLD, out$PVALUE, out$DISTANCE, 
			wt.folds, wt.pvalues, wt.dists)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	return( out)
}
