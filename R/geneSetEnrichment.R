# geneSetEnrichment.R

# take a subset of genes (like from a DE analysis, and a list of 'GeneSets'
# like a pathway groupings or GO groupings, and find the enrichment by these groups


geneSetEnrichment <- function( genes, geneSets=defaultGeneSets(),
				upOnly=TRUE, minEnrich=2.0, maxPvalue=0.05, wt.enrich=1, wt.pvalue=2,
				geneUniverse=NULL, reportGenes=FALSE, dropPathNameLinks=FALSE, verbose=T) {
	
	# the 'geneSets' can be a vector of strings of geneSet datasets..., 
	# or the usual list of vectors of genes
	cat( "\nGathering GeneSets..")
	GSanswer <- gatherGeneSets( geneSets, descriptor="geneSets", mode="combined")
	geneSets <- GSanswer$geneSets[[1]]



	getPathwaySetData <- function( pathwaySets, geneUniverse=NULL) {

		# given the list of all gene sets
		setNames <- names( pathwaySets)
		cnts <- sapply( pathwaySets, length)
		allGenes <- unlist( pathwaySets)
		allNames <- rep( setNames, times=cnts)
		pathTable <- data.frame( "GeneID"=allGenes, "PathwayID"= allNames, stringsAsFactors=F)
		rownames(pathTable) <- 1:nrow(pathTable)

		# we may be given just a subset of valid genes in the 'universe'
		if ( ! is.null( geneUniverse)) {
			#shrink the universe of all known genes
			keepers <- which( pathTable$GeneID %in% geneUniverse)
			pathTable <- pathTable[ keepers, ]
		}

		nGenesInTable <- length( unique( pathTable$GeneID))
		pathID.table <- table( pathTable$PathwayID)
		pathGenePcts <- pathID.table * 100 / nGenesInTable
		pathGeneCnts <- pathID.table
	
		cat( "\nFound ", length(pathID.table), " Gene Sets to evaluate")
		out <- list( "nRowsTotal"=nrow( pathTable), "nUniqueGenes"=nGenesInTable, 
				"PathwayTable"=pathTable, "idTable"=pathID.table, 
				"pctByGene"=pathGenePcts, "countByGene"=pathGeneCnts)
		return( out)
	}


	getGeneSetData <- function( genes, pathTable) {
	
		mySet <- subset( pathTable, GeneID %in% genes)
		if ( nrow(mySet) < 1) {
			cat( "\nWarning:  None of the given genes are found in any GeneSet")
			cat( "\nVerify the current species matches your genes.\n")
			out <- list( "nRowsTotal"=0, "nUniqueGenes"=0, 
					"PathwayTable"=data.frame(), "idTable"=vector(), 
					"pctByGene"=vector(), "countByGene"=vector())
			return( out)
		}
		rownames(mySet) <- 1:nrow(mySet)
	
		nGenesInSet <- length( unique( mySet$GeneID))
		pathID.table <- table( mySet$PathwayID)
		pathGenePcts <- pathID.table * 100 / nGenesInSet
		pathGeneCnts <- pathID.table
	
		out <- list( "nRowsTotal"=nrow( mySet), "nUniqueGenes"=nGenesInSet, 
				"PathwayTable"=mySet, "idTable"=pathID.table, 
				"pctByGene"=pathGenePcts, "countByGene"=pathGeneCnts)
		return( out)
	
	}
	
	
	findEnrich <- function( gSet, mData, gData, minEnrich=2.0, maxPvalue=0.1,
					upOnly=TRUE, wt.fold=1, wt.pvalue=1) {
	
		# given a set of genes, metabolic table count, and the gene set counts
		namesG <- names( gData$countByGene)
		namesM <- names( mData$countByGene)
	
		# lets look at all, not just those present
		# both <- intersect( namesG, namesM)
		both <- union( namesG, namesM)
		if ( length( both) < 1) return( data.frame())
		whereG <- match( both, namesG, nomatch=0)
		whereM <- match( both, namesM, nomatch=0)
	
		# gather the gene counts and percentages for each pathway
		Mcnt <- ifelse( whereM > 0, mData$countByGene[ whereM], 0)
		Mpct <- ifelse( whereM > 0, mData$pctByGene[ whereM], 0)
		Gcnt <- ifelse( whereG > 0, gData$countByGene[ whereG], 0)
		Gpct <- ifelse( whereG > 0, gData$pctByGene[ whereG], 0)
	
		# enrichment is how much more in our set than the entire table
		enrich <- Gpct / Mpct
	
		# the probalities are for each pathway...  calculate them
		pvals <- vector( length=length(enrich))
		cat( "\nCalculating enrichment P-values..\n")
		for ( i in 1:length(pvals)) {
			x <- Gcnt[i]
			m <- Mcnt[i]
			n <- mData$nUniqueGenes - m
			k <- gData$nUniqueGenes
			lowTailWanted <- (enrich[i] < 1.0)
			# get the entire prob. dist. and sum up the half we want...
			allPs <- dhyper( 0:k, m, n, k)
			if ( lowTailWanted) {
				# the probability of 0 to X genes
				pvals[i] <- sum( allPs[1:(x+1)])
			} else {
				# the probability of X up to K genes
				pvals[i] <- sum( allPs[(x+1):length(allPs)])
			}
			if (verbose && i %% 10 == 0) cat( "\r", i, sub( "<a.+","",both[i]), enrich[i], pvals[i])
		}
		out <- data.frame( both, enrich, pvals, Mcnt, Mpct, Gcnt, Gpct, stringsAsFactors=FALSE)
		colnames( out) <- c( "PathName", "Enrichment", "P_value", "N_Total", "Pct_Total", "N_Given", "Pct_Given")
	
		# having a gene count of zero forces a zero enrichment, which is hard to rank
		# use a small fudge factor to prevent exact zeros
		enrichDE <- enrich
		isZero <- which( Gcnt == 0)
		if ( length( isZero)) {
			smallCnt <- 1 / gData$nUniqueGenes
			Gpct[ isZero] <- smallCnt / gData$nUniqueGenes
			enrichDE[isZero] <- Gpct[isZero] / Mpct[isZero]
		}

		# sort by enrichment and P-value, where 1.0 means not enriched...
		#ord <- diffExpressRankOrder( out$Enrichment, out$P_value, wt.fold=wt.enrich, wt.pvalue=wt.pvalue,
		ord <- diffExpressRankOrder( enrichDE, out$P_value, wt.fold=wt.enrich, wt.pvalue=wt.pvalue,
					notDE.value=1.0)
		out <- out[ ord, ]
		rownames( out) <- 1:nrow( out)
	
		testValue <- out$Enrichment
		testP <- out$P_value

		who <- which( testValue >= minEnrich & testP <= maxPvalue)
		if ( ! upOnly) {
			who2 <- which( testValue <= (1/minEnrich) & testP <= maxPvalue)
			who <- unique( c( who, who2))
		}
		out <- out[ who, ]
		if ( nrow(out) > 0) rownames( out) <- 1:nrow( out)
		return( out)
	}

					
	whichGenesEnriched <- function( genes, pathTable, enrichDF, wt.fold=1, wt.pvalue=1,
					mode=c('vector','data.frame')) {

		mode <- match.arg( mode)
		foundGenes <- foundPaths <- foundEnrich <- foundPvals <- vector()
		foundNgiven <- foundNtotal <- vector()
		nout <- 0

				
		if ( mode == "data.frame") {
			if ( ! nrow(enrichDF)) return( data.frame())
			cat( "\nGathering genes in enriched pathways..\n")
			for (k in 1:nrow(enrichDF)) {
				myPath <- enrichDF$PathName[k]
				sml <- subset( pathTable, PathwayID == myPath)
				myGenes <- intersect( genes, sml$GeneID)
				N <- length(myGenes)
				if ( N) {
					now <- (nout+1) : (nout+N)
					foundGenes[now] <- myGenes
					foundPaths[now] <- myPath
					foundEnrich[now] <- enrichDF$Enrichment[k]
					foundPvals[now] <- enrichDF$P_value[k]
					foundNgiven[now] <- enrichDF$N_Given[k]
					foundNtotal[now] <- enrichDF$N_Total[k]
					nout <- max( now)
				}
				if ( k %% 10 == 0) cat( "\r", k, sub( "<a.+","",myPath), nout)
			}
			
			out <- data.frame( "GeneID"=foundGenes, "Product"=gene2Product(foundGenes),
					"PathName"=foundPaths, 
					"N_Given"=foundNgiven, "N_Total"=foundNtotal,
					"Enrichment"=foundEnrich, "P_value"=foundPvals, 
					stringsAsFactors=FALSE)
			#ord <- order( out$GeneID, -out$Enrichment, out$PathName)
			#ord <- diffExpressRankOrder( out$Enrichment, out$P_value, wt.fold=wt.enrich, 
			#			wt.pvalue=wt.pvalue, notDE.value=1.0)
			#out <- out[ ord, ]
			if ( nrow(out)) rownames(out) <- 1:nrow(out)

			# don't keep the hyperlinks
			out$PathName <- trimGeneSetNameLink( out$PathName)

			return( out)
		} else {
			if ( ! nrow(enrichDF)) return( vector())
			cat( "\nGathering genes in enriched pathways..\n")
			for (k in 1:nrow(enrichDF)) {
				myPath <- enrichDF$PathName[k]
				sml <- subset( pathTable, PathwayID == myPath)
				myGenes <- intersect( genes, sml$GeneID)
				N <- length(myGenes)
				now <- (nout+1)
				if ( N) {
					foundGenes[now] <- paste( sort( myGenes), collapse=", ")
				} else {
					foundGenes[now] <- ""
				}
				nout <- now
				if ( k %% 10 == 0) cat( "\r", k, sub( "<a.+","",myPath), nout)
			}
			return( foundGenes)
		}
	}

	
	# the pathway data is always the short gene names, never full IDs
	genes <- shortGeneName( genes, keep=1)
					
	setData <- getPathwaySetData( geneSets, geneUniverse=geneUniverse)
	geneData <- getGeneSetData( genes, setData$PathwayTable)
	enrichDF <- findEnrich( genes, setData, geneData, upOnly=upOnly, 
					minEnrich=minEnrich, maxPvalue=maxPvalue,
					wt.fold=wt.enrich, wt.pvalue=wt.pvalue)
	
	# with the enrichment result in hand, see what genes are the ones that hit
	if ( is.logical(reportGenes) && reportGenes == FALSE) {
		if (dropPathNameLinks) {
			enrichDF$PathName <- trimGeneSetNameLink( enrichDF$PathName)
		}
		return( enrichDF)
	}
	if ( is.logical(reportGenes) && reportGenes == TRUE) {
		enrichGenes <- whichGenesEnriched( genes, setData$PathwayTable, enrichDF, 
					wt.fold=wt.enrich, wt.pvalue=wt.pvalue, mode='data.frame')
		if (dropPathNameLinks) {
			enrichDF$PathName <- trimGeneSetNameLink( enrichDF$PathName)
			enrichGenes$PathName <- trimGeneSetNameLink( enrichGenes$PathName)
		}
		out <- list( "pathways"=enrichDF, "genes"=enrichGenes)
		return( out)
	}
	if ( is.character(reportGenes)) {
		enrichGenes <- whichGenesEnriched( genes, setData$PathwayTable, enrichDF, 
					wt.fold=wt.enrich, wt.pvalue=wt.pvalue, mode=reportGenes)
		if ( is.data.frame( enrichGenes)) {
			if (dropPathNameLinks) {
				enrichDF$PathName <- trimGeneSetNameLink( enrichDF$PathName)
				enrichGenes$PathName <- trimGeneSetNameLink( enrichGenes$PathName)
			}
			out <- list( "pathways"=enrichDF, "genes"=enrichGenes)
			return( out)
		}

		if (dropPathNameLinks) {
			enrichDF$PathName <- trimGeneSetNameLink( enrichDF$PathName)
		}
		out <- data.frame( enrichDF, "GeneList"=enrichGenes, stringsAsFactors=F)
		return( out)
	}
}


enrichment2html <- function( tbl, file="enrichment.html", title="Pathway Enrichment results", maxRows=500) {

	if ( ! nrow(tbl)) {
		cat( "\nEmpty Enrichment table. ", "\n Nothing to write to file: ", file)
		return()
	}

	# add a column to be more obvious
	tbl$Over_Under <- ifelse( tbl$Enrichment > 1.0, "Over", "Under")

	# make the names a bit easier to wrap around
	linkPart <- sub( "(.+)(<a.+)", "\\2", tbl$PathName)
	linkPart <- ifelse( linkPart == tbl$PathName, "", linkPart)
	namePart <- sub( "(.+)(<a.+)", "\\1", tbl$PathName)
	namePart <- gsub( "_", " ", namePart, fixed=T)
	# the dots in the module source are needed...
	#namePart <- gsub( ".", " ", namePart, fixed=T)
	tbl$PathName <- paste( namePart, linkPart)

	# do minor reformating
	#c( "PathName", "N_Total", "Pct_Total", "N_Given", "Pct_Given", "Enrichment", "P_value")
	colnames(tbl)[1] <- "Pathway Name"
	tbl$Pct_Total <- as.percent( tbl$Pct_Total, big.value=100, digits=2, sep="")
	tbl$Pct_Given <- as.percent( tbl$Pct_Given, big.value=100, digits=2, sep="")
	tbl$Enrichment <- formatC( tbl$Enrichment, format="f", digits=3)
	tbl$P_value <- formatC( tbl$P_value, format="e", digits=2)
	colnames(tbl) <- gsub( "_", " ", colnames(tbl), fixed=T)
	colnames(tbl) <- gsub( ".", " ", colnames(tbl), fixed=T)

	table2html( tbl, file=file, title=title, linkColumnNames=NULL, maxRows=maxRows)
	return()
}


geneSetMetaEnrichment <- function( genes, Ngenes=seq( 100, length(genes), by=100),
				geneSets=defaultGeneSets(),
				upOnly=TRUE, minEnrich=1.05, maxPvalue=0.5,
				wt.enrich=1, wt.pvalue=2, verbose=T) {
	
	# the Enrichment tool with multiple subsets and then do the meta ranking
	ansList <- vector( mode="list")
	nAns <- 0

	for ( n in Ngenes) {
		if ( n > length(genes)) next

		cat( "\n\nDoing Enrichment of top N = ", n)
		ans <- geneSetEnrichment( genes[1:n], geneSets=geneSets, upOnly=T, minEnrich=minEnrich, 
				maxPvalue=maxPvalue, wt.enrich=wt.enrich, wt.pvalue=wt.pvalue, 
				reportGenes=FALSE, geneUniverse=NULL, verbose=verbose)

		# let's augment the pathway name with how many genes are in it...
		ans$PathName <- paste( ans$PathName,"  (N=", ans$N_Total,")", sep="")

		nAns <- nAns + 1
		ansList[[nAns]] <- ans
		names(ansList)[nAns] <- paste( "N_", n, sep="")
	}
	cat( "\nDone.\n")

	ans <- metaRank.data.frames( ansList, geneColumn="PathName", valueColumn="Enrichment", pvalueColumn="P_value",
				missingGenes="fill", missingValue=1)

	# tidy up...
	colnames(ans)[1] <- "PathName"
	ans <- ans[ , -grep( "PRODUCT", colnames(ans))]

	return( ans)
}


metaEnrichment2html <- function( tbl, file="metaEnrichment.html", title="Pathway Meta Enrichment Results",
				maxRows=500) {

	# make the names a bit easier to wrap around
	linkPart <- sub( "(.+)(<a.+)", "\\2", tbl$PathName)
	linkPart <- ifelse( linkPart == tbl$PathName, "", linkPart)
	namePart <- sub( "(.+)(<a.+)", "\\1", tbl$PathName)
	namePart <- gsub( "_", " ", namePart, fixed=T)
	#namePart <- gsub( ".", " ", namePart, fixed=T)
	tbl$PathName <- paste( namePart, linkPart)

	# do minor reformating
	colnames(tbl)[1] <- "Pathway Name"
	tbl$Enrichment <- formatC( tbl$Enrichment, format="f", digits=3)
	tbl$AVG_RANK <- formatC( tbl$AVG_RANK, format="f", digits=2)
	colnames(tbl)[3] <- "Avg Rank"
	tbl$AVG_PVALUE <- formatC( tbl$AVG_PVALUE, format="e", digits=2)
	colnames(tbl)[4] <- "Avg Pvalue"
	colnames(tbl) <- gsub( "_", " ", colnames(tbl), fixed=T)
	colnames(tbl) <- gsub( ".", " ", colnames(tbl), fixed=T)

	table2html( tbl, file=file, title=title, linkColumnNames=NULL, maxRows=maxRows)
	return()
}

