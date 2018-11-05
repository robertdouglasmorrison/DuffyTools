# rankProduct.R -- test for differential expression

# from FEFS Letters:  573(2004) 83-92  Rainer Breiling, et.al.


rankProductDiffExpress <- function( fnames, fids, groupSet, targetGroup=groupSet[1], m=NULL,
			geneColumn="GENE_ID", intensityColumn="INTENSITY", 
			offset=0, keepIntergenics=FALSE, average.FUN=sqrtmean, extraColumn=NULL,
			poolSet=rep( 1, length(fnames)), nSimulations=100,
			missingGenes=c("drop", "fill")) {

	missingGenes <- match.arg( missingGenes)

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

	if ( ! keepIntergenics) {
		didDrop <- FALSE
		drops <- grep( "(ng)", rownames(m), fixed=T)
		if ( length(drops) > 0) {
			m <- m[ -drops, ]
			didDrop <- TRUE
		}
		nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
		drops <- which( rownames(m) %in% nonGenes)
		if ( length(drops) > 0) {
			m <- m[ -drops, ]
			didDrop <- TRUE
		}
		if (didDrop) cat( "\nDropped NonGenes..   N_TrueGenes: ", nrow(m))
	}

	allGenes <- rownames(m)
	nFiles <- ncol(m)
	nGenes <- length( allGenes)
	cat( "\nN_Files:  ", nFiles,"      N_Genes:  ", nGenes)

	# force a lower threshold to the intensities, to avoid divide by zero, or extreme fold changes.
	if ( offset > 0) {
		cat( "\nApplying linear offset to all intensities: ", offset)
		m <- m + offset
	}

	cat( "\n\nOrganizing...")
	allGroups <- groupSet
	allPools <- poolSet
	grpNames <- sort( unique( allGroups))
	if ( length( grpNames) < 2) stop( "'RankProduct' needs at least 2 groups of files")
	if ( ! (targetGroup %in% grpNames)) stop( paste("Not a given RankProduct group: ", targetGroup, 
				"   choices: ", paste(grpNames, collapse=", "))) 
	allInten <- m
	rownames( allInten) <- allGenes
	colnames( allInten) <- allGroups
	if ( ! is.null( extraColumn)) {
		allExtra <- mExtra
		rownames( allExtra) <- allGenes
		colnames( allExtra) <- allGroups
	}
	cat( "\n");  print( table( allGroups, dnn=NULL))

	# make the average gene intensity for each group
	avgInten1 <- apply( allInten[ , which( allGroups == targetGroup), drop=F], MARGIN=1, FUN=average.FUN, na.rm=TRUE)
	avgInten2 <- apply( allInten[ , which( allGroups != targetGroup), drop=F], MARGIN=1, FUN=average.FUN, na.rm=TRUE)
	if ( ! is.null( extraColumn)) {
		avgExtra1 <- apply( allExtra[ , which( allGroups == targetGroup), drop=F], MARGIN=1, FUN=average.FUN, na.rm=TRUE)
		avgExtra2 <- apply( allExtra[ , which( allGroups != targetGroup), drop=F], MARGIN=1, FUN=average.FUN, na.rm=TRUE)
	}
	if ( offset > 0) {
		avgInten1 <- avgInten1 - offset
		avgInten2 <- avgInten2 - offset
	}
	# lastly, never let the expression go negative
	avgInten1 <- pmax( avgInten1, 0)
	avgInten2 <- pmax( avgInten2, 0)

	# now we can do all pairs (from the K groups) and get a big set of fold change ranks
	cat( "\nFold Change over all pairs...")
	allRanks <- allFCs <- vector( mode="list")
	nRanks <- 0
	nPoolSkips <- 0
	for( i in 1:nFiles) {
	    if ( allGroups[i] != targetGroup) next
	    for( j in 1:nFiles) {
		# only do fold change between the 2 different groups
		if ( allGroups[i] == allGroups[j]) next
		# skip the pair if the 2 files are from different pools of samples
		if ( allPools[i] != allPools[j]) {
			nPoolSkips <- nPoolSkips + 1
			next
		}
		myFC <- allInten[ , i] / allInten[ , j]
		myFC <- log2( myFC)
		# turn fold change to ranks, we want "biggest" to be number one
		myRanks <- rank( (-myFC), ties.method="average")
		nRanks <- nRanks + 1
		allRanks[[ nRanks]] <- myRanks
		allFCs[[ nRanks]] <- myFC
	}}
	cat( "\nN_File pairs of FoldChange data: ", nRanks)
	if ( nPoolSkips > 0) cat( "\nN_File 'not same pool' pairs skipped: ", nPoolSkips)

	# with all pairs done, we can now do the Rank Product step
	allFCs <- matrix( unlist(allFCs), nrow=nGenes, ncol=nRanks)
	myFC <- apply( allFCs, MARGIN=1, FUN=mean)
	allRanks <- matrix( unlist(allRanks), nrow=nGenes, ncol=nRanks)
	avgRank <- apply( allRanks, MARGIN=1, FUN=average.FUN)

	# let's do it log based and scale by how many compares...
	myRPup <- apply( (allRanks/nGenes), MARGIN=1, FUN=average.FUN)

	# Rank Product is not symmetric for down regulation...  do that separate
	allRanksDown <- (nGenes - allRanks + 1)
	myRPdown <- apply( (allRanksDown/nGenes), MARGIN=1, FUN=average.FUN)
	avgRankDown <- apply( allRanksDown, MARGIN=1, FUN=average.FUN)

	gProds <- gene2ProductAllSpecies( allGenes)
	if ( any( gProds == "")) {
		tmp <- read.delim( fnames[1], as.is=T)
		gnams <- tmp[[ geneColumn]]
		gpros <- tmp[[ grep( "product", tolower(colnames(tmp)))]]
		needs <- which( gProds == "")
		where <- match( allGenes[needs], gnams, nomatch=0)
		gProds[ needs[ where > 0]] <- gpros[ where]
	}

	out <- data.frame( "GENE_ID"=allGenes, "PRODUCT"=gProds, "LOG2FOLD"=myFC,
			"RP_VALUE"=myRPup, "AVG_RANK"=avgRank, 
			"RP_VALUE_DOWN"=myRPdown, "AVG_RANK_DOWN"=avgRankDown, 
			"AVG_SET1"=avgInten1, "AVG_SET2"=avgInten2, stringsAsFactors=FALSE)
	if ( ! is.null(extraColumn)) {
		out$AVG_EXTRA1 <- avgExtra1
		out$AVG_EXTRA2 <- avgExtra2
	}
	ord <- order( out$RP_VALUE)
	out <- out[ ord, ]
	rownames(out) <- 1:nGenes

	# lastly, simulate to estimate an E value and false positives rate for these rank products
	cat( "\n\nE-value simulation...\n")
	bigRPset <- bigRPsetDown <- vector()
	for ( k in 1:nSimulations) {
		# each simulation repeats our K by N replicates
		allRanks <- vector( mode="list")
		allRanksDown <- vector( mode="list")
		nRanks <- 0
		# A) permute the intensities
		for( i in 1:nFiles) allInten[ , i] <- sample( allInten[ ,i], nGenes)
		# B) redo the RP calc...
		for( i in 1:nFiles) {
	    	    if ( allGroups[i] != targetGroup) next
		    myV1 <- allInten[ , i]
		    lapply( 1:nFiles, function(j) {
			if ( allGroups[i] == allGroups[j]) return()
			# skip the pair if the 2 files are from different pools of samples
			if ( allPools[i] != allPools[j]) return()
			myFC <- log2( myV1 / allInten[ , j])
			# turn fold change to ranks, we want "biggest" to be number one
			myRanks <- base::rank( (-myFC), ties.method="average")
			nRanks <<- nRanks + 1
			allRanks[[ nRanks]] <<- myRanks
			allRanksDown[[ nRanks]] <<- nGenes - myRanks + 1
			return()
		    })
		}
		allRanks <- matrix( unlist(allRanks), nrow=nGenes, ncol=nRanks)
		myRP <- apply( (allRanks/nGenes), MARGIN=1, FUN=average.FUN)
		bigRPset <- append( bigRPset, myRP)
		allRanksDown <- matrix( unlist(allRanksDown), nrow=nGenes, ncol=nRanks)
		myRPdown <- apply( (allRanksDown/nGenes), MARGIN=1, FUN=average.FUN)
		bigRPsetDown <- append( bigRPsetDown, myRPdown)
		cat( "\r",k)
	}
	cat( "\nFinalizing...")
	bigRPset <- sort( bigRPset)
	# estimate by finding location of something a 'tiny bit bigger' and then adjust
	useRP <- out$RP_VALUE * 1.00001
	nBetter <- findInterval( useRP, bigRPset)
	nBetter <- ifelse( nBetter > 0, (nBetter-1), 0)
	myEvalue <- nBetter / nSimulations
	PFP <- myEvalue / 1:nGenes
	PFP <- ifelse( PFP > 1, 1, PFP)
	out$E_VALUE <- myEvalue
	out$FP_RATE <- PFP

	# repeat for the down regulation parts
	bigRPsetDown <- sort( bigRPsetDown)
	ord <- order( out$RP_VALUE_DOWN)
	out <- out[ ord, ]
	useRP <- out$RP_VALUE_DOWN * 1.00001
	nBetter <- findInterval( useRP, bigRPsetDown)
	nBetter <- ifelse( nBetter > 0, (nBetter-1), 0)
	myEvalue <- nBetter / nSimulations
	PFP <- myEvalue / 1:nGenes
	PFP <- ifelse( PFP > 1, 1, PFP)
	out$E_VALUE_DOWN <- myEvalue
	out$FP_RATE_DOWN <- PFP

	# put back into UP order
	ord <- order( out$RP_VALUE)
	out <- out[ ord, ]

	# lastly, turn the colnames into the labels given, and re-order
	intenCol1 <- 8
	intenCol2 <- 9
	colnames(out)[intenCol1] <- targetGroup
	lab2 <- paste( "Not", targetGroup, sep=" ")
	colnames(out)[intenCol2] <- lab2
	if ( is.null( extraColumn)) {
		evalueCols <- 10:13
		out <- out[ ,c(1:7, evalueCols, intenCol1:intenCol2)]
	} else {
		extraCol1 <- 10
		extraCol2 <- 11
		evalueCols <- 12:15
		colnames(out)[extraCol1] <- paste( "Extra", targetGroup, sep=" ")
		colnames(out)[extraCol2] <- paste( "Extra Not", targetGroup, sep=" ")
		out <- out[ , c( 1:7, evalueCols, intenCol1:extraCol2)]
	}
	return( out)
}


calcRP <- function( x) {

	# let's do it log based for numerical stability and normalize by how many compares...
	return( 2 ^ (sum( log2(x)) / length(x)))
}


rankProduct <- function( rankM, nSimulations=1000, average.FUN=logmean) {

	# given a matrix of rank positions, where the rows are genes and the columns
	# are files, do the Rank Product calculation to get: average, E-value, FP rate
	nGenes <- nrow( rankM)
	nFiles <- ncol( rankM)

	avgRank <- apply( rankM, MARGIN=1, FUN=average.FUN)

	# turn the ranks in each row into their RP probabilites
	myRP <- apply( (rankM/nGenes), MARGIN=1, FUN=average.FUN)

	out <- data.frame( "RP_VALUE"=myRP, "AVG_RANK"=avgRank, stringsAsFactors=FALSE)
	rownames(out) <- rownames( rankM)

	# lastly, simulate to estimate an E value and false positives rate for these rank products
	cat( "\nE-value simulation...\n")
	bigRPset <- vector( length=(nGenes*nSimulations))
	nout <- 0
	rankMfake <- rankM
	for ( k in 1:nSimulations) {
		# A) permute each column of ranks
		for( i in 1:nFiles) rankMfake[ , i] <- sample( rankMfake[ ,i], nGenes)
		# B) redo the RP calc...
		myRP <- apply( (rankMfake/nGenes), MARGIN=1, FUN=average.FUN)
		bigRPset[(nout+1):(nout+nGenes)] <-  myRP
		nout <- nout + nGenes
		cat( "\r",k)
	}

	# to do the E-value and False Positive rate, we need to put the result in sorted order
	myorder <- order( out$RP_VALUE)
	out2 <- out[ myorder, ]
	cat( "\nFinalizing...")
	bigRPset <- sort( bigRPset)
	# estimate by finding location of something a 'tiny bit bigger' and then adjust
	useRP <- out2$RP_VALUE * 1.00001
	nBetter <- findInterval( useRP, bigRPset)
	nBetter <- ifelse( nBetter > 0, (nBetter-1), 0)
	myEvalue <- nBetter / nSimulations
	PFP <- myEvalue / 1:nGenes
	PFP <- ifelse( PFP > 1, 1, PFP)
	out2$E_VALUE <- myEvalue
	out2$FP_RATE <- PFP

	# revert this to the original given order
	redoOrder <- 1:nGenes
	redoOrder[ myorder] <- 1:nGenes

	return( out2[ redoOrder, ])
}
