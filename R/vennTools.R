# vennTools.R  -- assorted Venn Diagram tools


`vennOverlap` <- function( file1, file2, fid1=basename(file1), fid2=basename(file2), 
				gene1Column="GENE_ID", gene2Column=gene1Column, 
				value1Column="RPKM_M", value2Column=value1Column, 
				minValue1=1, minValue2=minValue1, extraColumns=NULL,
				minCount=NULL, maxCount=NULL, keepIntergenics=FALSE, 
				geneUniverse=NULL,
				label="Your label goes here...", sep1="\t", sep2=sep1, cex=2,
				col1='blue', col2='red') {

	df1 <- read.delim( file1, as.is=T, sep=sep1)
	cat( "\nRead file: ", file1, "\nN_Genes read in: ", nrow(df1))
	if ( !( all( c( gene1Column, value1Column) %in% colnames(df1)))) {
		cat( "\nMissing columns:  looked for: ", gene1Column, value1Column, "  Found: ", colnames(df1))
		return()
	}

	df2 <- read.delim( file2, as.is=T, sep=sep2)
	cat( "\nRead file: ", file2, "\nN_Genes read in: ", nrow(df2))
	if ( !( all( c( gene2Column, value2Column) %in% colnames(df2)))) {
		cat( "\nMissing columns:  looked for: ", gene2Column, value2Column, "  Found: ", colnames(df2))
		return()
	}
	
	return( vennOverlap.data.frames( df1, df2, fid1, fid2, 
				gene1Column=gene1Column, gene2Column=gene2Column, 
				value1Column=value1Column, value2Column=value2Column, 
				minValue1=minValue1, minValue2=minValue2, extraColumns=extraColumns,
				minCount=minCount, maxCount=maxCount, keepIntergenics=keepIntergenics, 
				geneUniverse=geneUniverse,
				label=label, cex=cex, col1=col1, col2=col2))
}


`vennOverlap.data.frames` <- function( df1, df2, fid1, fid2, 
				gene1Column="GENE_ID", gene2Column=gene1Column, 
				value1Column="RPKM_M", value2Column=value1Column, 
				minValue1=1, minValue2=minValue1, extraColumns=NULL,
				minCount=NULL, maxCount=NULL, keepIntergenics=FALSE, 
				geneUniverse=NULL, 
				label="Your label goes here...", cex=2, col1='blue', col2='red') {

	checkX11()

	# local function
	`gatherWantedSubset` <- function( df, geneColumn, valueColumn, minValue=1, minCount=NULL, maxCount=NULL) {

		if ( ! geneColumn %in% colnames(df)) {
			cat( "\nMissing GENE_ID column:  looked for: ", geneColumn, "\nFound: ", colnames(df))
			stop()
		}
		if ( ! is.null( valueColumn)) {
			if ( ! valueColumn %in% colnames(df)) {
				cat( "\nMissing VALUE column:  looked for: ", valueColumn, "\nFound: ", colnames(df))
				stop()
			}
			# order by the value column
			ord <- order( df[[valueColumn]], decreasing=T)
			df <- df[ ord, ]
		}

		# only the real genes?
		if ( ! keepIntergenics) {
			drops <- grep( "(ng)", df[[geneColumn]], fixed=T)
			if ( length( drops)) {
				cat( "\n  Dropped intergenics: ", length(drops))
				df <- df[ -drops, ]
			}
		}

		# main selection is by value
		if ( is.null( valueColumn)) {
			want <- 1:nrow(df)
		} else {
			want <- which( df[[valueColumn]] >= minValue)
		}
		if ( ! is.null( maxCount)) {
			if ( length(want) > maxCount) want <- want[ 1:maxCount]
		}
		if ( ! is.null( minCount)) {
			if ( length(want) < minCount) want <- union( want, 1:minCount)
		}
		return( df[ want, ])
	}

	smlDF1 <- gatherWantedSubset( df1, gene1Column, value1Column, minValue1, minCount=minCount, maxCount=maxCount)
	smlDF2 <- gatherWantedSubset( df2, gene2Column, value2Column, minValue2, minCount=minCount, maxCount=maxCount)
	genes1 <- smlDF1[[gene1Column]]
	genes2 <- smlDF2[[gene2Column]]
	if ( is.null( geneUniverse)) {
		geneUniverse <- union( df1[[gene1Column]], df2[[gene2Column]])
	} else {
		geneUniverse <- as.GeneUniverse( geneUniverse)
	}
	N_TotalGenes <- length( geneUniverse)

	# combine and resolve
	either <- union( genes1, genes2)
	both <- intersect( genes1, genes2)
	Neither <- length(either)
	Nboth <- length(both)
	N1 <- length(genes1)
	N2 <- length(genes2)
	N1plus2 <- N1 + N2
	Pct1 <- Nboth / N1
	Pct2 <- Nboth / N2
	PctBoth <- Nboth / Neither
	only1 <- setdiff( genes1, both)
	only2 <- setdiff( genes2, both)
	Nonly1 <- length( only1)
	Nonly2 <- length( only2)
	cat( "\nSet 1 = ", fid1, " \tSet 2 = ", fid2)
	cat( "\nN_Genes in union:      ",  Neither)
	cat( "\nN_Genes in common:     ",  Nboth)
	cat( "\nN_Genes only in set 1: ",  Nonly1)
	cat( "\nN_Genes only in set 2: ",  Nonly2)
	cat( "\nN_Genes in Universe:   ",  N_TotalGenes)
	cat( "\n")

	# create a joint meta ranked table of the gene intersection
	if ( length(both)) {
		dfList <- vector( mode="list")
		dfList[[1]] <- smlDF1
		dfList[[2]] <- smlDF2
		names(dfList) <- c( fid1, fid2)
		# do the meta ranks with all genes in both...
		metaAns <- metaRank.data.frames( dfList, geneColumn=gene1Column, valueColumn=value1Column,
						pvalueColumn="AVG_PVALUE", missingGenes="fill", verbose=F)
		# make sure the short version of gene names are present
		shortID <- shortGeneName( metaAns[[gene1Column]], keep=1)
		if ( any( shortID != metaAns[[gene1Column]])) {
			metaAns <- cbind( "GENE_NAME"=shortID, metaAns, stringsAsFactors=F)
		}

		# and then keep just the intersection, and both edge groups too...
		who <- which( metaAns[[gene1Column]] %in% both)
		intersectDF <- metaAns[ who, ]
		who <- which( metaAns[[gene1Column]] %in% only1)
		only1DF <- metaAns[ who, ]
		who <- which( metaAns[[gene1Column]] %in% only2)
		only2DF <- metaAns[ who, ]
	} else {
		intersectDF <- data.frame()
		only1DF <- smlDF1
		only2DF <- smlDF2
	}
	# we can be more precise than the Meta Rank tool, as to where these genes were in the other data
	where2 <- match( only1DF[[gene1Column]], df2[[gene2Column]], nomatch=0)
	only1DF[[ncol(only1DF)]][where2 > 0] <- where2[ where2 > 0]
	where1 <- match( only2DF[[gene2Column]], df1[[gene1Column]], nomatch=0)
	only2DF[[ncol(only2DF)-1]][where1 > 0] <- where1[ where1 > 0]

	# allow keeping extra columns from the input that we want to carry forward
	# do the merge in 'reverse" 2, then 1, order so they end up in the 1,2 order...
	if ( ! is.null(extraColumns)) {
		extras <- c( extraColumns, value1Column, value2Column)
	} else {
		extras <- c( value1Column, value2Column)
	}
	extras <- rev( unique( extras))

	# local function to do the merging
	mergeExtraColumn <- function( df, extrasDF, column, geneColumn="GENE_ID", extraGeneColumn=geneColumn, 
					suffix="Extra") {

		where <- column
		v <- extrasDF[[where]]
		extrasName <- colnames(extrasDF)[where]
		# get the gene IDs from the extras
		extraIDs <- extrasDF[[extraGeneColumn]]
		# find where those genes are in the result we merge into
		whG <- match( extraIDs, df[[geneColumn]], nomatch=0)
		# create the new column and fill it
		outV <- rep.int( NA, nrow(df))
		outV[whG] <- v[ whG > 0]
		# choose where to insert this new column
		colPtr <- max( match( c("GENE_ID","PRODUCT","GENE_NAME"), colnames(df), nomatch=0))
		if ( ! colPtr) colPtr <- nrow(df)
		# append after this site, and then add on the rest
		leftSide <- df[ , 1:colPtr, drop=F]
		newDF <- cbind( leftSide, outV, stringsAsFactors=F)
		colnames(newDF)[colPtr+1] <- paste( extrasName, suffix, sep="_")
		if (colPtr < ncol(df)) {
			rightSide <- df[ , (colPtr+1):ncol(df), drop=F]
			newDF <- cbind( newDF, rightSide, stringsAsFactors=F)
		}
		return( newDF)
	}

	for ( extraCol in extras) {
		where1 <- match( extraCol, colnames(df1), nomatch=0)
		where2 <- match( extraCol, colnames(df2), nomatch=0)
		if (where2 > 0 && length(both)) {
			intersectDF <- mergeExtraColumn( intersectDF, df2, column=where2, geneColumn=gene2Column, suffix=fid2)
			only1DF <- mergeExtraColumn( only1DF, df2, column=where2, geneColumn=gene2Column, suffix=fid2)
			only2DF <- mergeExtraColumn( only2DF, df2, column=where2, geneColumn=gene2Column, suffix=fid2)
		}
		if (where1 > 0 && length(both)) {
			intersectDF <- mergeExtraColumn( intersectDF, df1, column=where1, geneColumn=gene1Column, suffix=fid1)
			only1DF <- mergeExtraColumn( only1DF, df1, column=where1, geneColumn=gene1Column, suffix=fid1)
			only2DF <- mergeExtraColumn( only2DF, df1, column=where1, geneColumn=gene1Column, suffix=fid1)
		}
	}

	# create some descriptive text
	subText1 <- ""
	if ( ! is.null( value1Column)) subText1 <- paste( "Gene Selection Criteria:   Value of '", value1Column, "' >= ", minValue1, sep="")
	if ( ! is.null( minCount)) subText1 <- paste( subText1, "   Min_Genes = ", minCount, sep="")
	if ( ! is.null( maxCount)) subText1 <- paste( subText1, "   Max_Genes = ", maxCount, sep="")
	subTextGenome <- paste( "Enrichment:   Likelihood of  ", Nboth, "  in common:   P = ", 
				formatC(pvalGenome, format='e', digits=2),
				"    (Expected = ", formatC( nExpectGenome, format='f', digits=2), ")", sep="")

	# at this point we have everything we need to calc and plot the venns
	vennAns <- vennOverlap.counts( N1=N1, N2=N2, Nboth=Nboth, Ntotal=N_TotalGenes, label=label, 
			subText1=subText1, subText2=subTextGenome, cex=cex, col1=col1, col2=col2,
			fid1=fid1, fid2=fid2)

	return( list( "intersection"=intersectDF, "only1"=only1DF, "only2"=only2DF, 
			"overlap"=Nboth, "PercentOverlap"=PctBoth, "Percent1"=Pct1, "Percent2"=Pct2, 
			vennAns))
}


`vennOverlap.counts` <- function( N1, N2, Nboth, Ntotal, label="", subText1="", subText2="", col1='blue', col2='red', 
				cex=2, fid1="", fid2="") {

	enrichAnsGenome <- enrichment( nMatch=Nboth, nYourSet=N1, nTotal=Ntotal, nTargetSubset=N2)
	nExpectGenome <- enrichAnsGenome$nExpected
	if ( Nboth >= nExpectGenome) {
		pvalGenome <- enrichAnsGenome$P_atLeast
	} else {
		pvalGenome <- enrichAnsGenome$P_atMost
	}

	# let's calc how to draw these 2 circles
	r1 <- sqrt( N1 / pi)
	r2 <- sqrt( N2 / pi)
	#r1 <- r2 <- min( r1, r2)

	# start by assuming no overlap, and then bring them closer based on size of overlap
	x1 <- -r1
	x2 <- r2
	y <- 0
	powerFactor <- 1
	# when the pools are very different in size, we need to push them together a bit more
	ratio <- (min( N1,N2) / max( N1,N2)) ^ 0.05   #  0.33
	powerFactor <- 1 / (2 - ratio)
	pctOverlap1 <- Nboth / N1
	shift1 <- pctOverlap1 ^ powerFactor
	x1 <- x1 + if (shift1 > 0) (sqrt(shift1)*r1) else -(r1*0.1)
	pctOverlap2 <- Nboth / N2
	shift2 <- pctOverlap2 ^ powerFactor
	x2 <- x2 - if (shift2 > 0) (sqrt(shift2)*r2) else -(r2*0.1)

	# draw them
	require( plotrix)
	xlim <- max(c(r1,r2)*2) * c( -1.05,1.05)
	ylim <- max(r1,r2) * c( -1.05,1.05)
	mainText <- paste( "Venn Overlap:    ", label)
	plot( 1,1, type='n', main=paste( "Venn Overlap:    ", label), xaxt='n', yaxt='n', xlab=NA, ylab=NA,
			xlim=xlim, ylim=ylim, frame.plot=FALSE)
	if (subText1 != "") mtext( subText1, side=1, line=1, font=2, cex=1.1)
	if (subText2 != "") mtext( subText2, side=1, line=2, font=2, cex=1.1)

	draw.circle( x1, y, r1, nv=200, border=col1, lwd=4)
	draw.circle( x2, y, r2, nv=200, border=col2, lwd=4)

	# show the gene counts
	Nonly1 <- N1 - Nboth
	Nonly2 <- N2 - Nboth
	Neither <- Nonly1 + Nonly2 + Nboth
	text( 0, y, Nboth, col=1, font=2, cex=cex)
	text( ((x1-r1) + (x2-r2))/2, y, Nonly1, col=col1, font=2, cex=cex)
	text( ((x1+r1) + (x2+r2))/2, y, Nonly2, col=col2, font=2, cex=cex)
	asPcts <- gsub( " ", "", as.percent( c(Nboth, Nonly1, Nonly2), big.value=Neither, digits=0))
	yPct <- y - (max(r1,r2)*0.13)
	text( 0, yPct, asPcts[1], col=1, font=2, cex=cex*0.5)
	text( ((x1-r1) + (x2-r2))/2, yPct, asPcts[2], col=col1, font=2, cex=cex*0.5)
	text( ((x1+r1) + (x2+r2))/2, yPct, asPcts[3], col=col2, font=2, cex=cex*0.5)

	text( x1-r1*0.5, r1*1.005, fid1, col=col1, font=2, cex=cex*0.5)
	text( x2+r2*0.5, r2*1.005, fid2, col=col2, font=2, cex=cex*0.5)

	dev.flush()
	return( list( "expected"=nExpectGenome, "p.value"=pvalGenome))
}


`vennOverlap.list` <- function( listOfGeneVectors, geneUniverse=NULL) {

	# given one list, with arbitrary number of vectors, try to do a fully generic Venn Overlap
	N <- length( listOfGeneVectors)
	if ( N < 2) return(NULL)
	setName <- names( listOfGeneVectors)

	# step over all possible counts of sets of genes
	# record the group compare string and count as we go
	outCompText <- outOverlap <- outFirstOnly <- outSecondOnly <- outPval <- vector()
	whoOverlap <- vector()
	nOut <- 0

	# set up to do enrichment analysis on the overlaps, of what genes could be seen
	if ( is.null( geneUniverse)) {
		gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)
		geneUniverse <- unique( shortGeneName( gmap$GENE_ID, keep=1))
	}
	nTotal <- length( geneUniverse)

	# create all the permutations, generically.  Use a call to expand.grid()
	oneSeqText <- paste( 1, N, sep=":")
	nSeqText <- paste( "list(", paste( rep( oneSeqText, times=N), collapse=", "), ")", sep="")
	permM <- expand.grid( eval( parse( text=nSeqText)))

	# this is a data.frame of N columns, with every possible permutation of 1 to N
	for ( set1size in 1:(N-1)) {
		# grab the first K columns to be set 1
		set1M <- as.matrix( permM[ , 1:set1size, drop=F])
		for ( set2size in 1:(N-set1size)) {
			# grab the next K columns to be set 2
			set2M <- as.matrix( permM[ , (1:set2size)+set1size, drop=F])
			# now visit every row of the permutations
			for ( k in 1:nrow(permM)) {
				is1 <- sort( unique( set1M[ k, ]))
				n1 <- length(is1)
				is2 <- sort( unique( set2M[ k, ]))
				n2 <- length(is2)
				# if any sets are in both lists, skip this compare
				if ( length( intersect( is1, is2)) > 0) next
				text1 <- paste( "(", paste( setName[is1], collapse="+"), ")", sep="")
				text2 <- paste( "(", paste( setName[is2], collapse="+"), ")", sep="")
				# make the smaller set be first
				if ( length(is1) < length(is2)) {
					thisCompText <- paste( text1, text2, sep=".vs.")
				} else if ( length(is1) > length(is2)) {
					thisCompText <- paste( text2, text1, sep=".vs.")
				} else {
					# when both sets same size, have the 'earlier' one go first
					if ( is1[1] < is2[1]) {
						thisCompText <- paste( text1, text2, sep=".vs.")
					} else {
						thisCompText <- paste( text2, text1, sep=".vs.")
					}
				}
				# if we have this already, skip it
				if ( thisCompText %in% outCompText) next
				# we have a new compare, so get those genes.  
				#Take the intersection of all genes in each set of sets
				genes1 <- listOfGeneVectors[[ is1[1]]]
				if ( n1 > 1) for ( j in 2:n1) genes1 <- intersect( genes1, listOfGeneVectors[[ is1[j]]])
				genes2 <- listOfGeneVectors[[ is2[1]]]
				if ( n2 > 1) for ( j in 2:n2) genes2 <- intersect( genes2, listOfGeneVectors[[ is2[j]]])
				# we have a new compare to store results for
				nOut <- nOut + 1
				outCompText[nOut] <- thisCompText
				outOverlap[nOut] <- nInBoth <- length( genesInBoth <- intersect( genes1, genes2))
				outFirstOnly[nOut] <- length( setdiff( genes1, genes2))
				outSecondOnly[nOut] <- length( setdiff( genes2, genes1))
				whoOverlap[nOut] <- if ( nInBoth) paste( sort( genesInBoth), collapse="; ") else ""
				# calc the enrichment of the overlap, and extract the appropriate P-value
				enrichAns <- enrichment( nMatch=nInBoth, nYourSet=length(genes1), nTotal=nTotal, nTargetSubset=length(genes2))
				outPval[nOut] <- myPval <- if ( enrichAns$nExpected < nInBoth) enrichAns$P_atLeast_N else enrichAns$P_atMost_N
				# when either of the two groups were themselves a overlap of previous groups, the P-values need to be combined
				if ( any( c(n1,n2) > 1)) {
					if ( n1 > 1) {
						# undo the 'combining' text, and find that earlier overlap comparison
						grpTextString <- sub( "+", ").vs.(", text1, fixed=T)
						wh <- match( grpTextString, outCompText, nomatch=0)
						if (wh > 0) myPval <- myPval * outPval[wh]
					}
					if ( n2 > 1) {
						grpTextString <- sub( "+", ").vs.(", text2, fixed=T)
						wh <- match( grpTextString, outCompText, nomatch=0)
						if (wh > 0) myPval <- myPval * outPval[wh]
					}
					outPval[nOut] <- myPval
				}
			}
		}
	}

	out <- data.frame( "Groups.Compared"=outCompText, "Overlap"=outOverlap, "Group1.Only"=outFirstOnly, 
				"Group2.Only"=outSecondOnly, "P.Value"=outPval, "Who.In.Overlap"=whoOverlap, 
				stringsAsFactors=F)
	rownames(out) <- 1:nOut

	# when given more than 2 groups, we need one final intersection of all groups
	if ( N > 2) {
		genesAll <- listOfGeneVectors[[ 1]]
		for ( j in 2:N) genesAll <- intersect( genesAll, listOfGeneVectors[[ j]])
		nOverlapAll <- length( genesAll)
		whoAll <- if ( nOverlapAll) paste( sort( genesAll), collapse="; ") else ""
		allCompText <- paste( "All(", paste( setName, collapse="+"), ")", sep="")
		out2 <- data.frame( "Groups.Compared"=allCompText, "Overlap"=nOverlapAll, "Group1.Only"=NA, 
				"Group2.Only"=NA, "P.Value"=NA, "Who.In.Overlap"=whoAll, stringsAsFactors=F)
		
		out <- rbind( out, out2)
		rownames(out) <- 1:nrow(out)
	}

	return( out)
}
