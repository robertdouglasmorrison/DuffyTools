# metaRank.R -- combine several ranked files of genes

 
metaRanks <- function( fnames, fids, weightset=rep(1, length(fnames)), 
			geneColumn="GENE_ID", valueColumn="LOG2FOLD", 
			pvalueColumn="PVALUE", productColumn="PRODUCT", sep="\t",
			rank.average.FUN=sqrtmean, value.average.FUN=mean,
			keepIntergenics=FALSE, missingGenes=c("drop", "fill", "na"), 
			missingValue=0, naDropPercent=0.5, nFDRsimulations=0,
			diffExpressPvalueCorrection=FALSE, 
			direction=c("UP","DOWN"), cleanHyperlinks=FALSE, verbose=TRUE) {

	missingGenes <- match.arg( missingGenes)
	direction <- match.arg( direction)

	if (verbose) cat( "\nReading in files...")
	allDF <- vector( "mode"="list")
	allGenes <- vector()
	nfiles <- 0
	missing <- vector()
	for (i in 1:length(fnames)) {
		filename <- fnames[i]
		if ( ! file.exists( filename)) {
			cat( "\nNot found: ", filename)
			missing <- c( missing, i)
			next
		}
		tmp <- read.delim( filename, as.is=T, sep=sep)
		if ( ! geneColumn %in% colnames(tmp)) {
			cat( "\nGeneID column not found.   File=", filename)
			next
		}
		nfiles <- nfiles + 1
		if ( ! keepIntergenics) {
			isNG <- grep( "(ng)", tmp[[ geneColumn]], fixed=T)
			if ( length(isNG) > 0) tmp <- tmp[ -isNG, ]
			nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
			drops <- which( tmp[[ geneColumn]] %in% nonGenes)
			if ( length(drops) > 0) tmp <- tmp[ -isNG, ]
		}
		# optionally strip out hyperlinks, like in GeneSet PathName usages...
		if (cleanHyperlinks) {
			tmp[[ geneColumn]] <- cleanGeneSetModuleNames( tmp[[ geneColumn]], wrapParentheses=F)
		}

		# by default, all the files are assumed to be already sorted in "Up-regulated" order
		allDF[[ nfiles]] <- tmp
		if ( direction == "DOWN") {
			# the general rule is to make it be "Down-regulated", flip it over
			nr <- nrow( tmp)
			tmp <- tmp[ rev( 1:nr), ]
			# special cases need special treatment
			# 1)  Rank Product is not symmetrical.  Instead it has explicit columns for down
			if ( "RP_VALUE_DOWN" %in% colnames(tmp)) {
				tmp <- tmp[ order( tmp$RP_VALUE_DOWN), ]
			}
			allDF[[ nfiles]] <- tmp
		}
		allGenes <- base::union( allGenes, tmp[[geneColumn]])
		if (verbose) cat( "\n",filename, "\tN_Genes: ", nrow(tmp))
	}

	allG <- sort( unique( allGenes))
	NG <- length( allG)

	if ( length( missing) > 0) {
		fnames <- fnames[ -missing]
		fids <- fids[ -missing]
	}
	if (verbose) cat( "\nN_Files: ", nfiles, "   N_Genes: ", NG)
	
	# get the rank order in all datasets
	rankM <- foldM <- pvalM <- matrix( NA, nrow=NG, ncol=nfiles)
	rownames(rankM) <- rownames(foldM) <- rownames(pvalM) <- allG
	colnames(rankM) <- colnames(foldM) <- colnames(pvalM) <- fids
	allProds <- rep( "", NG)

	for( i in 1:nfiles) {
		thisDF <- allDF[[i]]
		theseGenes <- thisDF[[geneColumn]]
		where <- match( allG, theseGenes, nomatch=0)
		# so DE tools throw away genes with no difference and/or expression
		# they will be 'missing', so give them values saying 'not differentially expressed'
		isMissing <- which( where == 0)
		missingFold <- 0
		missingPval <- 1
		missingRank <-  round( (NG+nrow(thisDF))/4)  # thus the 'midpoint' of the average # of genes
		rankM[ where > 0, i] <- where[ where > 0]
		rankM[ isMissing, i] <- missingRank
		logFoldColumn <- grep( valueColumn, colnames(thisDF))
		if ( length( logFoldColumn) > 0) {
			foldM[ where > 0, i] <- thisDF[[ logFoldColumn[1]]][where]
			foldM[ isMissing, i] <- missingFold
		}
		pvalColumn <- grep( pvalueColumn, colnames(thisDF))
		if ( length( pvalColumn) > 0) {
			pvalM[ where > 0, i] <- thisDF[[ pvalColumn[1]]][where]
			pvalM[ isMissing, i] <- missingPval
		}
		hasPROD <- ( !is.na(productColumn) && (productColumn %in% colnames(thisDF)))
		if ( hasPROD) {
			allProds[ where > 0] <- ifelse( allProds[where > 0] == "", 
							thisDF[[ productColumn]][where], allProds[ where > 0])
		}
	}

	# check for missing genes
	whoNA <- vector()
	for ( i in 1:nfiles) {
		who <- which( is.na( rankM[ ,i]))
		if ( length(who) > 0) {
			whoNA <- sort( unique( c( whoNA, who)))
			if ( missingGenes == "fill") {
				rankM[ who, i] <- nrow(rankM)
				foldM[ who, i] <- missingValue
				pvalM[ who, i] <- 1.0
			}
		}
	}
	if ( missingGenes == "drop" && length(whoNA) > 0) {
		rankM <- rankM[ -whoNA, ]
		foldM <- foldM[ -whoNA, ]
		pvalM <- pvalM[ -whoNA, ]
		allG <- allG[ -whoNA]
		allProds <- allProds[ -whoNA]
		NG <- nrow(rankM)
	}

	# if 'na', and a gene row has too many na's, drop the whole thing
	if ( missingGenes == "na") {
		nna <- apply( rankM, 1, function(x) sum(is.na(x)))
		whoNA <- which( nna > (nfiles * naDropPercent))
		if ( length( whoNA) > 0) {
			rankM <- rankM[ -whoNA, ]
			foldM <- foldM[ -whoNA, ]
			pvalM <- pvalM[ -whoNA, ]
			allG <- allG[ -whoNA]
			allProds <- allProds[ -whoNA]
			NG <- nrow(rankM)
		}
	}
	if (verbose) cat( "    after 'drop' and/or 'na' filtering: ", NG)

	if (verbose) cat( "\nAveraging Ranks...")
	if ( all( weightset == 1)) {
		avgRank <- apply( rankM, MARGIN=1, FUN=rank.average.FUN, na.rm=T)
	} else {
		if ( identical( rank.average.FUN, logmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( log2(x) * weightset, na.rm=T)
						denom <- sum( weightset, na.rm=T)
						return( 2 ^ (numerator / denom))
					})
		} else if ( identical( rank.average.FUN, sqrtmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( sqrt(x) * weightset, na.rm=T)
						denom <- sum( weightset, na.rm=T)
						return( (numerator / denom) ^ 2)
					})
		} else {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( x * weightset, na.rm=T)
						denom <- sum( weightset, na.rm=T)
						return( numerator / denom)
					})
		}
	}

	if (diffExpressPvalueCorrection) {
		# the p-values of strongly down-regulated genes are (usually) very small numbers
		# if wanted, turn these into bad p-values of being up-regulated
		for ( i in 1:nfiles) {
			myfolds <- foldM[ ,i]
			mypvals <- pvalM[ ,i]
			if ( direction == "UP") {
				needsCorrect <- which( myfolds < 0 & mypvals < 0.5)
			} else {
				needsCorrect <- which( myfolds > 0 & mypvals < 0.5)
			}
			if ( length(needsCorrect)) {
				mypvals[ needsCorrect] <- 1.0 - mypvals[ needsCorrect] 
				pvalM[ ,i] <- mypvals
			}
		}
	}

	avgFold <- apply( foldM, MARGIN=1, FUN=value.average.FUN, na.rm=T)
	if (pvalueColumn != "") {
		#avgPval <- apply( pvalM, MARGIN=1, FUN=logmean, na.rm=T)
		avgPval <- apply( pvalM, MARGIN=1, FUN=p.combine)
	} else {
		avgPval <- rep.int( NA, nrow(pvalM))
	}

	if ( ! is.na( productColumn)) {
		out <- data.frame( allG, allProds, avgFold, avgPval, avgRank, rankM, stringsAsFactors=FALSE)
		colnames(out) <- c( geneColumn, productColumn, valueColumn, "AVG_PVALUE", "AVG_RANK", colnames( rankM))
	} else {
		out <- data.frame( allG, avgFold, avgPval, avgRank, rankM, stringsAsFactors=FALSE)
		colnames(out) <- c( geneColumn, valueColumn, "AVG_PVALUE", "AVG_RANK", colnames( rankM))
	}

	# do the final ordering by Average Rank
	if ( pvalueColumn != "") {
		ord <- order( out$AVG_RANK, out$AVG_PVALUE, -out[[ valueColumn]])
	} else {
		ord <- order( out$AVG_RANK, -out[[ valueColumn]])
	}
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)

	# do a simulation of random permutations of these ranks
	# do the FDR after sorting final row order
	if ( nFDRsimulations > 0) {
		simM <- rankM
		NR <- nrow(simM)
		randomAvgRank <- vector( length=NR*nFDRsimulations)
		nnow <- 0
		if (verbose) cat( "  estimating FDR..")
		for ( k in 1:nFDRsimulations) {
			for ( i in 1:nfiles) simM[ , i] <- sample( NR)
			randomNow <- apply( simM, MARGIN=1, FUN=rank.average.FUN, na.rm=T)
			randomAvgRank[ (nnow+1):(nnow+NR)] <- randomNow
			nnow <- nnow + NR
			if (verbose) cat( ".")
		}
		# with this pool of 'by-chance average ranks, we can estimate the likelihood of ours
		randomAvgRank <- sort( randomAvgRank)
		avgRank <- out$AVG_RANK
		myLocs <- findInterval( avgRank * 1.00000001, randomAvgRank)
		myLocs <- ifelse( myLocs > 0, myLocs - 1, 0)
		myEvalue <- myLocs / nFDRsimulations
		myFPrate <- myEvalue / (1:NR)
		myFPrate <- ifelse( myFPrate > 1, 1, myFPrate)
		out$FDR <- round( myFPrate, digits=4)
	}

	# correlation test...
	ccM <- matrix( NA, nrow=nfiles, ncol=nfiles)
	for( i in 1:(nfiles-1)) {
	for( j in (i+1):nfiles) {
		thisCC <- cor( rankM[ ,i], rankM[ ,j], use="complete")
		ccM[i,j] <- ccM[j,i] <- thisCC
	}}
	colnames( ccM) <- colnames(rankM)
	cc <- apply( ccM, MARGIN=2, mean, na.rm=T)
	metaRankCC <<- cc

	if (verbose) {
		cat( "\nRank Correlations:\n")
		print( sort( cc, decreasing=T))
	}

	return( out)
}


metaRank2html <- function( tbl, fileout="metaRanks.html", title="", maxRows=100, 
			valueColumn="LOG2FOLD", ...) {

	# clean up any formatting...
	if ( nrow(tbl)) {
		tbl[[ valueColumn]] <- formatC( tbl[[ valueColumn]], format="f", digits=3)
		if ("AVG_PVALUE" %in% colnames(tbl)) tbl$AVG_PVALUE <- formatC( tbl$AVG_PVALUE, format="e", digits=3)
		if ("AVG_RANK" %in% colnames(tbl)) tbl$AVG_RANK <- formatC( tbl$AVG_RANK, format="f", digits=2)
		NC <- ncol(tbl)
		if ( NC > 3) colnames(tbl)[4:NC] <- gsub( "_", " ", colnames(tbl)[4:NC], fixed=T)
	}

	title <- paste( "Meta Ranks:  &nbsp; ", title)
	table2html( tbl, fileout=fileout, title=title, maxRows=maxRows, ...)
	return()
}


metaRank.data.frames <- function( df.list, weightset=rep(1, length(df.list)), 
			geneColumn="GENE_ID", valueColumn="LOG2FOLD", 
			pvalueColumn="PVALUE", productColumn="PRODUCT",
			rank.average.FUN=sqrtmean, value.average.FUN=mean,
			missingGenes=c("drop", "fill", "na"), missingValue=0,
			naDropPercent=0.5, diffExpressPvalueCorrection=FALSE, 
			direction="UP", nFDRsimulations=0, cleanHyperlinks=FALSE, 
			verbose=TRUE) {

	missingGenes <- match.arg( missingGenes)

	allDF <- vector( "mode"="list")
	nDF <- 0
	allGenes <- vector()
	nGenesPerDF <- vector()
	for (i in 1:length(df.list)) {
		tmp <- df.list[[ i]]
		if ( is.null(tmp)) {
			cat( "\nData frame is NULL  DataFrame: ", i, "\t", names(df.list)[i])
			next
		}
		if ( nrow(tmp) < 1) {
			cat( "\nData frame is empty...  DataFrame: ", i, "\t", names(df.list)[i])
			next
		}
		if ( ! geneColumn %in% colnames(tmp)) {
			cat( "\nGeneID column not found.   DataFrame: ", i, "\t", names(df.list)[i])
			next
		}
		# optionally strip out hyperlinks, like in GeneSet PathName usages...
		if (cleanHyperlinks) {
			tmp[[ geneColumn]] <- cleanGeneSetModuleNames( tmp[[ geneColumn]], wrapParentheses=F)
		}
		nDF <- nDF + 1
		allDF[[ nDF]] <- tmp
		names(allDF)[nDF] <- names(df.list)[i]
		nGenesPerDF[nDF] <- nrow(tmp)
		allGenes <- base::union( allGenes, tmp[[geneColumn]])
		if (verbose) cat( "\n",i, names(df.list)[i], "\tN_Genes: ", nrow(tmp))
	}
	fids <- names( allDF)

	allG <- sort( unique( allGenes))
	NG <- length( allG)
	if (verbose) cat( "\nN_DataFrames: ", nDF, "   N_Genes: ", NG)
	
	# get the rank order in all datasets
	rankM <- foldM <- pvalM <- matrix( NA, nrow=NG, ncol=nDF)
	rownames(rankM) <- rownames(foldM) <- rownames(pvalM) <- allG
	colnames(rankM) <- colnames(foldM) <- colnames(pvalM) <- fids
	allProds <- rep( "", NG)

	for( i in 1:nDF) {
		thisDF <- allDF[[i]]
		hasPROD <- ( !is.na(productColumn) && (productColumn %in% colnames(thisDF)))
		theseGenes <- thisDF[[geneColumn]]
		where <- match( allG, theseGenes, nomatch=0)
		rankM[ where > 0, i] <- where[ where > 0]
		logFoldColumn <- grep( valueColumn, colnames(thisDF))
		if ( length( logFoldColumn) > 0) {
			foldM[ where > 0, i] <- thisDF[[ logFoldColumn[1]]][where]
		}
		pvalColumn <- grep( pvalueColumn, colnames(thisDF))
		if ( length( pvalColumn) > 0) {
			pvalM[ where > 0, i] <- thisDF[[ pvalColumn[1]]][where]
		}
		if (hasPROD) {
			allProds[ where > 0] <- ifelse( allProds[where > 0] == "", 
				thisDF[[ productColumn]][where], allProds[ where > 0])
		}
	}

	# check for missing genes
	whoNA <- vector()
	for ( i in 1:nDF) {
		who <- which( is.na( rankM[ ,i]))
		if ( length(who) > 0) {
			whoNA <- sort( unique( c( whoNA, who)))
			if ( missingGenes == "fill") {
				rankM[ who, i] <- nrow(rankM)
				foldM[ who, i] <- missingValue
				pvalM[ who, i] <- 1.0
			}
		}
	}
	if ( missingGenes == "drop" && length(whoNA) > 0) {
		rankM <- rankM[ -whoNA, , drop=F]
		foldM <- foldM[ -whoNA, , drop=F]
		pvalM <- pvalM[ -whoNA, , drop=F]
		allG <- allG[ -whoNA]
		allProds <- allProds[ -whoNA]
		NG <- length(allG)
	}

	# if 'na', and a gene row has too many na's, drop the whole thing
	if ( missingGenes == "na") {
		nna <- apply( rankM, 1, function(x) sum(is.na(x)))
		whoNA <- which( nna > (nDF * naDropPercent))
		if ( length( whoNA) > 0) {
			rankM <- rankM[ -whoNA, , drop=F]
			foldM <- foldM[ -whoNA, , drop=F]
			pvalM <- pvalM[ -whoNA, , drop=F]
			allG <- allG[ -whoNA]
			allProds <- allProds[ -whoNA]
			NG <- length(allG)
		}
	}
	if (verbose) cat( "    after 'drop' and/or 'na' filtering: ", NG)
	if ( NG < 1) return( data.frame())

	if (verbose) cat( "\nAveraging Ranks...")
	if ( all( weightset == 1)) {
		avgRank <- apply( rankM, MARGIN=1, FUN=rank.average.FUN, na.rm=T)
	} else {
		if ( identical( average.FUN, logmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( log2(x) * weightset)
						denom <- sum( weightset)
						return( 2 ^ (numerator / denom))
					})
		} else if ( identical( average.FUN, sqrtmean)) {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( sqrt(x) * weightset)
						denom <- sum( weightset)
						return( (numerator / denom) ^ 2)
					})
		} else {
			avgRank <- apply( rankM, MARGIN=1, FUN=function(x) {
						numerator <- sum( x * weightset)
						denom <- sum( weightset)
						return( numerator / denom)
					})
		}
	}

	if (diffExpressPvalueCorrection && direction == "UP") {
		# the p-values of strongly down-regulated genes are (usually) very small numbers
		# if wanted, turn these into bad p-values of being up-regulated
		for ( i in 1:nDF) {
			myfolds <- foldM[ ,i]
			mypvals <- pvalM[ ,i]
			needsCorrect <- which( myfolds < 0 & mypvals < 0.5)
			if ( length(needsCorrect)) {
				mypvals[ needsCorrect] <- 1.0 - mypvals[ needsCorrect] 
				pvalM[ ,i] <- mypvals
			}
			cat( "\nDebug: ", i, "  N_Corrected_Pvalues: ", length( needsCorrect))
		}
	}

	avgFold <- apply( foldM, MARGIN=1, FUN=value.average.FUN, na.rm=T)
	#avgPval <- apply( pvalM, MARGIN=1, FUN=logmean, na.rm=T)
	avgPval <- apply( pvalM, MARGIN=1, FUN=p.combine)
	if ( ! is.na(productColumn)) {
		out <- data.frame( allG, allProds, avgFold, avgRank, avgPval, rankM, stringsAsFactors=FALSE)
		colnames(out) <- c( geneColumn, productColumn, valueColumn, "AVG_RANK", "AVG_PVALUE", colnames( rankM))
	} else {
		out <- data.frame( allG, avgFold, avgRank, avgPval, rankM, stringsAsFactors=FALSE)
		colnames(out) <- c( geneColumn, valueColumn, "AVG_RANK", "AVG_PVALUE", colnames( rankM))
	}

	# do the final ordering by Average Rank
	ord <- order( out$AVG_RANK, out$AVG_PVALUE, -out[[ valueColumn]])
	out <- out[ ord, ]
	rownames( out) <- 1:nrow(out)

	# do a simulation of random permutations of these ranks
	# do the FDR after sorting final row order
	if ( nFDRsimulations > 0) {
		simM <- rankM
		NR <- nrow(simM)
		randomAvgRank <- vector( length=NR*nFDRsimulations)
		nnow <- 0
		if (verbose) cat( "  estimating FDR..")
		for ( k in 1:nFDRsimulations) {
			for ( i in 1:nDF) simM[ , i] <- sample( NR)
			randomNow <- apply( simM, MARGIN=1, FUN=rank.average.FUN, na.rm=T)
			randomAvgRank[ (nnow+1):(nnow+NR)] <- randomNow
			nnow <- nnow + NR
			if (verbose) cat( ".")
		}
		# with this pool of 'by-chance average ranks, we can estimate the likelihood of ours
		randomAvgRank <- sort( randomAvgRank)
		avgRank <- out$AVG_RANK
		myLocs <- findInterval( avgRank * 1.00000001, randomAvgRank)
		myLocs <- ifelse( myLocs > 0, myLocs - 1, 0)
		myEvalue <- myLocs / nFDRsimulations
		myFPrate <- myEvalue / (1:NR)
		myFPrate <- ifelse( myFPrate > 1, 1, myFPrate)
		out$FDR <- round( myFPrate, digits=4)
	}

	# correlation test...
	ccM <- matrix( NA, nrow=nDF, ncol=nDF)
	if ( nDF > 1) {
		for( i in 1:(nDF-1)) {
		for( j in (i+1):nDF) {
			thisCC <- cor( rankM[ ,i], rankM[ ,j], use="complete")
			ccM[i,j] <- ccM[j,i] <- thisCC
		}}
	}
	colnames( ccM) <- colnames(rankM)
	cc <- apply( ccM, MARGIN=2, mean, na.rm=T)
	metaRankCC <<- cc

	if (verbose) {
		cat( "\nRank Correlations:\n")
		print( sort( cc, decreasing=T))
	}

	return( out)
}
