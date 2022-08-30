# geneSetTools.R -- various tools for manipulating 'GeneSets'


`defaultGeneSets` <- function( speciesID=getCurrentSpecies()) {

	# expand this to look at the cell type defintions, to add custom names to this list
	out <- vector()
	if ( speciesID %in% MAMMAL_SPECIES) out <-  MAMMAL_GENE_SETS
	if ( speciesID %in% PARASITE_SPECIES) out <- PARASITE_GENE_SETS
	if ( speciesID %in% BACTERIA_SPECIES) out <- BACTERIA_GENE_SETS

	# look at the current cell type reference name
	reference <- getCellTypeReference()
	if ( ! is.null( reference)) {
		out2 <- paste( reference, c( "AllGeneSets", "TopGeneSets"), sep=".")
		out <- c( out, out2)
	}

	return( out)
}


`gatherGeneSets` <- function( geneSets, descriptor="CombinedGeneSets", mode=c("combined", "separate")) {

	mode <- match.arg( mode)

	bigSetOfGeneSets <- vector( "mode"="list")
	bigSetOfDescriptors <- vector()
	nbig <- 0

	if ( is.character( geneSets)) {
		# we can be given a set of data object names, that hold an object named 'allGeneSets'
		for( nam in geneSets) {
			file <- paste( getCurrentSpeciesFilePrefix(), nam, sep=".")
			allGeneSets <- NULL
			data( list=file, envir=environment())
			if ( !is.null(allGeneSets)) {
				nbig <- nbig + 1
				bigSetOfDescriptors[nbig] <- nam
				bigSetOfGeneSets[[nbig]] <- allGeneSets
			}
		}
	} else {
		# old way is one explicit list
		nbig <- 1
		bigSetOfDescriptors[1] <- descriptor
		bigSetOfGeneSets[[1]] <- geneSets
	}

	# combine the separate geneSets if we want to...
	if ( mode == "combined" && nbig > 1) {
		newBig <- vector( mode="list")
		for ( i in 1:nbig) {
			thisName <- bigSetOfDescriptors[i]
			thisSet <- bigSetOfGeneSets[[i]]
			myNames <- names( thisSet)
			myNewNames <- paste( thisName, myNames, sep=": ")
			names(thisSet) <- myNewNames
			newBig <- c( newBig, thisSet)
		}
		nbig <- 1
		bigSetOfDescriptors[1] <- descriptor
		bigSetOfGeneSets[[1]] <- newBig
		length(bigSetOfDescriptors) <- length(bigSetOfGeneSets) <- 1
	}

	out <- list( "n"=nbig, "descriptor"=bigSetOfDescriptors, "geneSets"=bigSetOfGeneSets)
	return( out)
}


`cleanGeneSetName` <- function( text) {

	# make the names a bit easier to wrap around in a html page
	# but preserve any HTML link
	linkPart <- sub( "(.+)(<a.+)", "\\2", text)
	linkPart <- ifelse( linkPart == text, "", linkPart)
	namePart <- sub( "(.+)(<a.+)", "\\1", text)
	namePart <- gsub( "_", " ", namePart, fixed=T)
	newtext <- paste( namePart, linkPart)
	newtext <- sub( " +$", "", newtext)
	return( newtext)
}


`cleanGeneSetModuleNames` <- function( nams, wrapParentheses=TRUE) {

	# make the names a bit shorter, and remove any HTML link
	out <- nams

	# 1. strip out any hyperlink anchors
	out <- sub( " ?<a.+/a>", "", out)
	# 1.5 strip out any hypertext 'blank space'
	out <- sub( " ?\\&nbsp;", "", out)

	# 2.  any () get a second line and may get clipped
	if (wrapParentheses) {
		hasParen <- grep( "(", out, fixed=T)
		fronts <- sub( "\\(.+", "", out[hasParen])
		backs <- sub( "(.+)(\\(.+)", "\\2", out[hasParen])
		doClip <- which( nchar( backs) > 40)
		backs[ doClip] <- paste( substr( backs[doClip], 1, 37), "...)", sep="")
		out[ hasParen] <- paste( fronts, backs, sep="\n")
	}
	out <- sub( " +$", "", out)
	out
}


# make a table of top gene set names from a set of genes
`geneSetHits` <- function( genes, geneSets=defaultGeneSets()) {
	
	# the 'geneSets' can be a vector of strings of geneSet datasets..., 
	# or the usual list of vectors of genes
	cat( "\nGathering GeneSets..")
	GSanswer <- gatherGeneSets( geneSets, descriptor="geneSets", mode="combined")
	geneSets <- GSanswer$geneSets[[1]]
	cat( "  Done.\n")

	# we don't want any HTML links for this step...
	names(geneSets) <- trimGeneSetNameLink( names( geneSets))


	getPathwaySetData <- function( pathwaySets) {

		# given the list of all gene sets
		setNames <- names( pathwaySets)
		cnts <- sapply( pathwaySets, length)
		allGenes <- unlist( pathwaySets)
		allNames <- rep( setNames, times=cnts)
		pathTable <- data.frame( "GeneID"=allGenes, "PathwayID"= allNames, stringsAsFactors=F)
		rownames(pathTable) <- 1:nrow(pathTable)

		nGenesInTable <- length( unique.default( pathTable$GeneID))
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
	
		nGenesInSet <- length( unique.default( mySet$GeneID))
		pathID.table <- table( mySet$PathwayID)
		pathGenePcts <- pathID.table * 100 / nGenesInSet
		pathGeneCnts <- pathID.table
	
		out <- list( "nRowsTotal"=nrow( mySet), "nUniqueGenes"=nGenesInSet, 
				"PathwayTable"=mySet, "idTable"=pathID.table, 
				"pctByGene"=pathGenePcts, "countByGene"=pathGeneCnts)
		return( out)
	
	}
	
	
	setData <- getPathwaySetData( geneSets)
	geneData <- getGeneSetData( genes, setData$PathwayTable)

	return( geneData)
}


`geneSetHitsToPie` <- function( hits, min.genes=1, min.pct=1.5, max.show=20, label="", allGeneSets=NULL, 
				txt.padX=0, txt.padY=0,
				sort.order=c("alphabetical","increasing.value","decreasing.value")) {

	# take the output of 'geneSetHits()' and render as a Pie chart
	countData <- hits$countByGene
	nams <- names(countData)
	pcts <- countData * 100 / hits$nUniqueGenes

	# if the 'allGeneSets' full set of pathways gets passed in, use that for coloring
	# so all subsets retain the parent coloring
	if ( ! is.null( allGeneSets)) {
		allColors <- rainbow( length(allGeneSets), end=0.86)
		setNames <- trimGeneSetNameLink( names( allGeneSets))
		where <- match( nams, setNames, nomatch=0)
		if ( any( where == 0)) {
			cat( "\nError: some gene set names not found: ", names(countData)[where == 0])
		}
		mycolors <- allColors[where]
	} else {
		mycolors <- rainbow( length(countData), end=0.86)
	}

	# two possible reasons to not show all wedges, too few genes and too small a percentage
	keep1 <- which( countData >= min.genes)
	keep2 <- which( pcts >= min.pct)
	keep <- intersect( keep1, keep2)
	countData <- countData[ keep]
	nams <- nams[ keep]
	pcts <- pcts[ keep]
	mycolors <- mycolors[ keep]
	N <- length(countData)

	# arrange the ordering of the wedges
	sort.order <- match.arg( sort.order)
	if ( sort.order == "alphabetical") ord <- order( nams)
	if ( sort.order == "increasing.value") ord <- order( countData, decreasing=F)
	if ( sort.order == "decreasing.value") ord <- order( countData, decreasing=T)
	countData <- countData[ ord]
	nams <- nams[ ord]
	pcts <- pcts[ ord]
	mycolors <- mycolors[ ord]
	if ( N > max.show) {
		N <- length(countData) <- length(nams) <- length(pcts) <- length(mycolors) <- max.show
	}

	# append the percentage to the names...
	nams <- paste( nams, " (", round(pcts), "%)", sep="")
	names(countData) <- nams

	# push the labels a bit to minimize overlap, at the top (1,N) and at the bottom
	txt.padX <- rep( txt.padX, length.out=N)
	txt.padY <- rep( txt.padY, length.out=N)
	txt.padX[ 1] <- 0.18
	txt.padY[ 1] <- 0.07
	txt.padY[ N] <- 0.07
	txt.padX[ N] <- 0.18
	txt.padX[ 2] <- 0.16
	txt.padX[ N-1] <- 0.16
	whoSmaller <- if (countData[1] < countData[N]) 1 else N
	txt.padY[ whoSmaller] <- 0.16
	whoBottom <- which.min( abs( (cumsum(countData)-countData/2) - sum(countData)/2))
	txt.padY[ whoBottom] <- 0.12
	txt.padY[ whoBottom-1] <- 0.03
	txt.padY[ whoBottom+1] <- 0.03
	txt.padX[ whoBottom-1] <- 0.10
	txt.padX[ whoBottom+1] <- 0.10

	geneSetPie( countData, col=mycolors, radius=0.7, init.angle=90, edges=300, border='white', 
			main=label, txt.padX=txt.padX, txt.padY=txt.padY)

}


`geneSetPie` <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
			init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
			col = NULL, border = NULL, lty = NULL, main = NULL, txt.padX = 0, txt.padY = 0) 
{
	if (!is.numeric(x) || any(is.na(x) | x < 0)) 
	    stop("'x' values must be positive.")
	if (is.null(labels)) 
	    labels <- as.character(seq_along(x))
	else labels <- as.graphicsAnnot(labels)

	x <- c(0, cumsum(x)/sum(x))
	dx <- diff(x)
	nx <- length(dx)
	plot.new()
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1)
	if (pin[1L] > pin[2L]) 
	    xlim <- (pin[1L]/pin[2L]) * xlim
	else ylim <- (pin[2L]/pin[1L]) * ylim
	dev.hold()
	on.exit(dev.flush())
	plot.window(xlim, ylim, "", asp = 1)
	if (is.null(col)) 
	    col <- if (is.null(density)) 
	        c("white", "lightblue", "mistyrose", "lightcyan", 
	            "lavender", "cornsilk")
	    else par("fg")
	col <- rep(col, length.out = nx)
	border <- rep(border, length.out = nx)
	angle <- rep(angle, length.out = nx)
	#lty <- rep(lty, length.out = nx)
	#density <- rep(density, length.out = nx)
	twopi <- if (clockwise) 
	    -2 * pi
	else 2 * pi
	t2xy <- function(t) {
	    t2p <- twopi * t + init.angle * pi/180
	    list(x = radius * cos(t2p), y = radius * sin(t2p))
	}
	for (i in 1L:nx) {
	    n <- max(2, floor(edges * dx[i]))
	    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
	    polygon(c(P$x, 0), c(P$y, 0), density = density, angle = angle[i], 
	        border = border[i], col = col[i], lty = lty)
	    P <- t2xy(mean(x[i + 0:1]))
	    lab <- as.character(labels[i])
	    if (!is.na(lab) && nzchar(lab)) {
	        lines( (c(1, 1.15+txt.padX[i]) * P$x), (c(1, 1.25+txt.padY[i]) * P$y))
	        text( (1.2+txt.padX[i]) * P$x, (1.35+txt.padY[i]) * P$y, labels[i], xpd = TRUE, 
                adj = if (P$x < -0.15) 1 else if (P$x > 0.15) 0 else 0.5, font=2)
	    }
	}
	title( main=main)
	title( sub="Percent of Genes per Group", font.sub=2, line=0.75)
	invisible(NULL)
}


`geneSetGrep` <- function( pattern, geneSets=defaultGeneSets(), fixed=FALSE, geneProducts=FALSE) {
	
	# the 'geneSets' can be a vector of strings of geneSet datasets..., 
	# or the usual list of vectors of genes
	cat( "\nGathering GeneSets..")
	GSanswer <- gatherGeneSets( geneSets, descriptor="geneSets", mode="combined")
	geneSets <- GSanswer$geneSets[[1]]
	cat( "  Done.\n")

	# we don't want any HTML links for this step...
	geneSetNames <- names(geneSets) <- trimGeneSetNameLink( names( geneSets))

	hits <- grep( pattern, geneSetNames, fixed=fixed)
	if ( ! length( hits)) {
		# try case insensitive
		hits <- grep( tolower(pattern), tolower(geneSetNames), fixed=fixed)
	}
	if ( ! length( hits)) return( NULL)

	# default is to return the list of pathways and the genes as a vector
	if ( ! geneProducts) return( geneSets[hits])

	# optionally, return a list of small data frames that include the gene product
	annoSets <- geneSets[ hits]
	for ( i in 1:length(annoSets)) {
		setIn <- annoSets[[i]]
		setOut <- data.frame( "GENE_ID"=setIn, "PRODUCT"=gene2Product(setIn), stringsAsFactors=F)
		annoSets[[i]] <- setOut
	}
	return( annoSets)
}


`getCellTypeGeneSetAssociation` <- function( reference=getCellTypeReference(), optionsFile="Options.txt", 
					speciesID=getCurrentSpecies()) {

	geneSetCellTypes <- NULL

	# we are now using the generic cell type tools to know which resource to load
	CellTypeSetup( reference=reference)
	reference <- getCellTypeReference()

	# see if this is what is already loaded
	if ( exists( "GeneSetAssociation", envir=CellTypeEnv)) {
		geneSetCellTypes <- CellTypeEnv[[ "GeneSetAssociation"]]
	} else {
		geneSetFile <- paste( reference, "GeneSetAssociation", sep=".")
		localFile <- paste( geneSetFile, "rda", sep=".")
		if ( file.exists( localFile)) {
			cat( "\n  loading gene set association from local file: ", localFile)
			load( localFile, envir=environment())
		} else {
			data( list=geneSetFile, envir=environment())
		}
		# save a copy for future use
		if ( ! is.null( geneSetCellTypes)) CellTypeEnv[[ "GeneSetAssociation"]] <- geneSetCellTypes
	}
	return( geneSetCellTypes)
}


`geneSetCellType` <- function( geneSetNames, max.types=1, reference=getCellTypeReference()) {

	out <- rep.int( "", N <- length(geneSetNames))
	if ( ! N) return(out)

	# we are now using the generic cell type tools to know which resource to load
	geneSetCellTypes <- getCellTypeGeneSetAssociation( reference=reference)
	if ( is.null( geneSetCellTypes)) return(out)

	# prep what was passed in:  1) clean the names, and then strip any module HTML links
	namesIn <- cleanGeneSetName( geneSetNames)
	namesIn <- cleanGeneSetModuleNames( namesIn, wrapParentheses=F)
	knownGeneSets <- geneSetCellTypes$GeneSetName
	# next try is to strip off any gene set prefix from the names
	knownGeneSets2 <- sub( "^[A-Za-z0-9.]+: ", "", knownGeneSets)
	#knownGeneSetsB <- cleanGeneSetName( knownGeneSets)

	# we can now have 2+ cell types, so perhaps change how we ask
	if ( max.types == 1) {
		# first try is direct perfect match
		where <- match( namesIn, knownGeneSets, nomatch=0)
		out[ where > 0] <- geneSetCellTypes$CellType[ where]
		stillBlank <- which( out == "")
		if ( length( stillBlank)) {
			where2 <- match( namesIn[stillBlank], knownGeneSets2, nomatch=0)
			out[ stillBlank[ where2 > 0]] <- geneSetCellTypes$CellType[ where2]
		}
	} else {
		
		# find 1+ hits for each gene set
		smlCT <- subset( geneSetCellTypes, GeneSetName %in% namesIn)
		smlGenesetFac <- factor( smlCT$GeneSetName)
		whereOut <- match( levels(smlGenesetFac), namesIn)
		k <- 0
		if ( nrow(smlCT)) {
			tapply( 1:nrow(smlCT), smlGenesetFac, function(x) {
				if ( length(x) > max.types) x <- x[ 1:max.types]
				cellTypeStr <- paste( smlCT$CellType[x], ":", smlCT$PctExpression[x], "%", sep="", collapse="; ")
				k <<- k + 1
				out[ whereOut[k]] <<- cellTypeStr
				return(NULL)
			})
		}
	}
	out
}


`geneSetOverlap` <- function( gs1, gs2) {

	# do set overlap calculation (Jaccard Index and hypergeometric P-value)

	# allow the genesets to be a character string name
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.character(gs1)) {
		allGeneSets <- NULL
		f1 <- paste( prefix, gs1[1], sep=".")
		data( list=f1, envir=environment())
		if ( is.null( allGeneSets)) stop( paste( "GeneSet not found: ", f1))
		gs1 <- allGeneSets
	} else if ( ! is.list(gs1)) {
		stop( "GeneSet must be a list or a character string")
	}
	if ( is.character(gs2)) {
		allGeneSets <- NULL
		f2 <- paste( prefix, gs2[1], sep=".")
		data( list=f2, envir=environment())
		if ( is.null( allGeneSets)) stop( paste( "GeneSet not found: ", f2))
		gs2 <- allGeneSets
	} else if ( ! is.list(gs2)) {
		stop( "GeneSet must be a list or a character string")
	}
		

	# set up to visit all pairs of gene sets
	names1 <- cleanGeneSetModuleNames( names(gs1), wrap=F)
	names2 <- cleanGeneSetModuleNames( names(gs2), wrap=F)
	N1 <- length( gs1)
	N2 <- length( gs2)
	outID1 <- outID2 <- outN1 <- outN2 <- outNover <- outJaccard <- outExpect <- outPval <- vector()
	nout <- 0

	# set the universe of all possible gene names
	all1 <- unique.default( unlist( gs1))
	all2 <- unique.default( unlist( gs2))
	nAllGenes <- length( union( all1, all2))

	for ( i in 1:N1) {
		genes1 <- gs1[[i]]
		len1 <- length( genes1)
		if ( ! len1) next
		for ( j in 1:N2) {
			genes2 <- gs2[[j]]
			len2 <- length( genes2)
			if ( ! len2) next

			nInter <- length( intersect( genes1, genes2))
			nUnion <- len1 + len2 - nInter
			ji <- nInter / nUnion

			ans <- enrichment( nMatch=nInter, nYourSet=len1, nTotal=nAllGenes, nTargetSubset=len2)
			nExp <- ans$nExpected
			pv <- if (nExp >= nInter) ans$P_atMost_N else ans$P_atLeast_N

			nout <- nout + 1
			outID1[nout] <- names1[i]
			outID2[nout] <- names2[j]
			outN1[nout] <- len1
			outN2[nout] <- len2
			outNover[nout] <- nInter
			outJaccard[nout] <- ji
			outExpect[nout] <- nExp
			outPval[nout] <- pv
		}
	}

	out <- data.frame( "GeneSet1"=outID1, "GeneSet2"=outID2, "Jaccard"=round( outJaccard, digits=4), 
			"Pvalue"=outPval, "N_Overlap"=outNover, "N_Expect"=round( outExpect,digits=2),
			"N_Set1"=outN1, "N_Set2"=outN2, stringsAsFactors=F)
	ord <- order( out$Jaccard, decreasing=T)
	out <- out[ ord, ]
	keep <- which( out$N_Overlap > 0 & out$Pvalue <= 0.05)
	out <- out[ keep, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)

	return( out)
}


`geneSetBestMatch` <- function( genes, geneSets=defaultGeneSets(), nBest=5, speciesID=getCurrentSpecies()) {

	# do set overlap calculation (Jaccard Index and hypergeometric P-value)
	curID <- getCurrentSpecies()
	if ( speciesID != curID) {
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( curID))
	}

	# for a given set of gene names, and a GeneSet name(s)
	if ( ! is.list( geneSets)) {
		tmp <- gatherGeneSets( geneSets)
		geneSets <- tmp$geneSets[[1]]
	}

	genes1 <- shortGeneName( genes, keep=1)
	genes1 <- unique.default( genes1)
	len1 <- length( genes1)

	# set up to visit all the gene sets
	names2 <- cleanGeneSetModuleNames( names(geneSets), wrap=F)
	N2 <- length( geneSets)
	outID2 <- outN1 <- outN2 <- outNover <- outJaccard <- outEnrich <- outExpect <- outPval <- vector()
	nout <- 0

	# set the universe of all possible gene names
	all1 <- unique.default( genes1)
	all2 <- unique.default( unlist( geneSets))
	nAllGenes <- length( union( all1, all2))

	for ( j in 1:N2) {
		genes2 <- geneSets[[j]]
		len2 <- length( genes2)
		if ( ! len2) next

		nInter <- length( intersect( genes1, genes2))
		nUnion <- len1 + len2 - nInter
		ji <- nInter / nUnion

		ans <- enrichment( nMatch=nInter, nYourSet=len1, nTotal=nAllGenes, nTargetSubset=len2)
		nExp <- ans$nExpected
		pv <- if (nExp >= nInter) ans$P_atMost_N else ans$P_atLeast_N

		nout <- nout + 1
		outID2[nout] <- names2[j]
		outN1[nout] <- len1
		outN2[nout] <- len2
		outNover[nout] <- nInter
		outJaccard[nout] <- ji
		outExpect[nout] <- nExp
		outPval[nout] <- pv
		outEnrich[nout] <- ans$Enrichment
	}

	out <- data.frame( "BestGeneSet"=outID2, "Jaccard"=round( outJaccard, digits=5), 
			"Pvalue"=formatC( outPval, format="e", digits=3), 
			"N_Given"=outN1, "N_InSet"=outN2, "N_Overlap"=outNover, "N_Expect"=round( outExpect,digits=3), 
			"Enrichment"=round( outEnrich, digits=3), stringsAsFactors=F)
	#ord <- order( out$Jaccard, decreasing=T)
	ord <- order( out$Jaccard, -(as.numeric(out$Pvalue)), out$Enrichment, out$N_Overlap, decreasing=T)
	out <- out[ ord, ]
	keep <- which( out$N_Overlap > 0)
	if ( length( keep) > nBest) length(keep) <- nBest
	out <- out[ keep, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)

	return( out)
}


`geneSetFitDataFrame` <- function( tbl, geneSets=defaultGeneSets(), speciesID=getCurrentSpecies(),
				geneColumn="GENE_ID", foldColumn="LOG2FOLD", min.fold=0.2, 
				pvalueColumn="AVG_PVALUE", max.pvalue=0.05, max.geneSets=100) {

	# do set overlap calculation (Jaccard Index and hypergeometric P-value)
	curID <- getCurrentSpecies()
	if ( speciesID != curID) {
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( curID))
	}

	# for a given data frame of Diff Expression, see what gene sets best explain the UP genes
	if ( ! geneColumn %in% colnames(tbl)) {
		cat( "\n'GeneID' column not found.  Tried: ", geneColumn)
		cat( "\nColumn names seen: ", colnames(tbl))
		return( NULL)
	}
	if ( ! foldColumn %in% colnames(tbl)) {
		cat( "\n'LogFold' column not found.  Tried: ", foldColumn)
		cat( "\nColumn names seen: ", colnames(tbl))
		return( NULL)
	}
	if ( ! pvalueColumn %in% colnames(tbl)) {
		cat( "\n'P-value' column not found.  Tried: ", pvalueColumn)
		cat( "\nColumn names seen: ", colnames(tbl))
		return( NULL)
	}

	cat( "\nTrimming to significant UP genes...")
	cat( "\nFold change cutoff:    ", min.fold)
	cat( "\nP-value cutoff:        ", max.pvalue)
	keepFold <- which( as.numeric( tbl[[foldColumn]]) >= min.fold)
	keepPval <- which( as.numeric( tbl[[pvalueColumn]]) <= max.pvalue)
	keep <- sort( intersect( keepFold, keepPval))
	tbl <- tbl[ keep, ]
	obsGenes <- shortGeneName( tbl[[ geneColumn]], keep=1)
	obsFolds <- tbl[[foldColumn]]
	# undo the log2 transform, so we have "how many times over expected" each gene was seen
	obsFolds <- 2 ^ obsFolds

	NG <- length( obsGenes)
	cat( "\nN_Genes call Sig & UP: ", NG)

	# get all the gene sets
	cat( "\nGathering GeneSets and evaluating overlap..")
	if ( ! is.list( geneSets)) {
		tmp <- gatherGeneSets( geneSets)
		geneSets <- tmp$geneSets[[1]]
	}
	gsNames <- names(geneSets) <- cleanGeneSetModuleNames( names(geneSets), wrap=F)
	# find the gene sets that have the best overlap with these genes
	gsAns <- geneSetBestMatch( obsGenes, geneSets=geneSets, speciesID=speciesID, nBest=max.geneSets)
	#keep those gene sets that may be useful
	keep <- which( gsAns$Jaccard > 0.001 & as.numeric(gsAns$Pvalue) < 0.2 & gsAns$N_Overlap >= 2)
	NGS <- length( keep)
	if ( NGS < 2) {
		cat( "\nNot enough gene sets overlap with selected genes.")
		cat( "\nLoosen 'min.fold' and/or 'max.pvalue'")
		return(NULL)
	}
	gsAns <- gsAns[ keep, ]
	rownames(gsAns) <- 1:NGS
	# and map from these good overlappers back to the full set
	whereGS <- match( gsAns$GeneSet, names(geneSets))
	cat( "  Done.\nN_GeneSets with usable overlap: ", NGS)

	# build a matrix of 'gene detected' calls for every gene in every useful gene set
	# zeros for genes not in each gene set, +1 for genes in the gene set
	geneCalls <- matrix( 0, nrow=NG, ncol=NGS)
	rownames(geneCalls) <- obsGenes
	colnames(geneCalls) <- make.names( 1:NGS)
	for ( i in 1:NGS) {
		thisSet <- geneSets[[ whereGS[i] ]]
		wh <- match( obsGenes, thisSet, nomatch=0)
		geneCalls[ wh > 0, i] <- 1
	}

	# OK, ready to do a linear model of gene sets to best describe the fold change we see
	cat( "\nCalling 'stepAIC()' to find the best linear model..")
	require( MASS)
	formula <- reformulate( termlabels=colnames(geneCalls), response="obsFolds", intercept=FALSE)
	lmAns1 <- lm( formula, data=as.data.frame(geneCalls))
	lmAns2 <- stepAIC( lmAns1, trace=0)
	cat( "  Done.")

	# extract the coefficients
	cat( "\nReformat and package up results..")
	sumAns <- summary.lm( lmAns2)
	finalAns <- coef( sumAns)
	# clean them a bit
	finalAns[,1] <- round( finalAns[,1], digits=4)
	finalAns[,2] <- round( finalAns[,2], digits=3)
	finalAns[,3] <- round( finalAns[,3], digits=3)
	colnames(finalAns) <- c( "LM_Estim", "LM_StErr", "LM_Tval", "LM_Pval")
	# we only need/want the Estimate & P-value
	finalAns <- finalAns[ ,c(1,4)]
	# put the gene set names back in place
	newnames <- rownames(finalAns)
	wh <- match( newnames, colnames(geneCalls), nomatch=0)
	wh2 <- match( newnames, colnames(geneCalls), nomatch=NA)
	newnames[ wh > 0] <- gsAns$GeneSet[wh]
	rownames(finalAns) <- newnames
	# also grab some of the Jaccard details
	otherBits <- gsAns[ wh2, c(2,4:7)]
	colnames(otherBits) <- c( "Jaccard", "Enrich", "Overlap", "Expect", "SetSize")

	# prepare the final ordering
	# 1)  line 1 is the model origin
	isIntercept <- grep( "^\\(Intercept", rownames(finalAns))
	if ( length(isIntercept)) {
		finalAns <- finalAns[ -isIntercept, ]
		otherBits <- otherBits[ -isIntercept, ]
	}
	# 2)  we only want the positive coefficients, those that say being present added
	isUP <- which( finalAns[ ,1] > 0.0)
	if ( ! length( isUP)) {
		cat( "\nModel found no UP gene sets.")
		return(NULL)
	}
	finalAns <- finalAns[ isUP, , drop=FALSE]
	otherBits <- otherBits[ isUP, , drop=FALSE]
	# 3)  we only want the terms that were significant
	#isSIG <- which( finalAns[ ,2] <= 0.06)
	#if ( ! length( isSIG)) {
	#	cat( "\nModel found no significant gene sets.")
	#	return(NULL)
	#}
	#finalAns <- finalAns[ isSIG, ]
	#otherBits <- otherBits[ isSIG, ]
	# 3)  use our standard fold + P-value to rank these
	ord <- diffExpressRankOrder( finalAns[,1], finalAns[,2], wt.fold=1, wt.pvalue=1)
	finalAns <- finalAns[ ord, ]
	otherBits <- otherBits[ ord, ]

	out <- data.frame( "GeneSet"=rownames(finalAns), finalAns, otherBits, stringsAsFactors=F)
	rownames(out) <- 1:nrow(out)
	cat( "  Done.\n")

	return( out)
}


`evaluateGeneSets` <- function( tbl, geneSets=defaultGeneSets(), speciesID=getCurrentSpecies(),
				geneColumn="GENE_ID", foldColumn="LOG2FOLD", max.pvalue.keep=0.05) {

	# evaluate the "average" fold change for every gene set, given a set of gene fold change data
	curID <- getCurrentSpecies()
	if ( speciesID != curID) {
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( curID))
	}

	# for a given data frame of Diff Expression, see what gene sets best explain the fold changes we see
	if ( ! geneColumn %in% colnames(tbl)) {
		cat( "\n'GeneID' column not found.  Tried: ", geneColumn)
		cat( "\nColumn names seen: ", colnames(tbl))
		return( NULL)
	}
	if ( ! foldColumn %in% colnames(tbl)) {
		cat( "\n'LogFold' column not found.  Tried: ", foldColumn)
		cat( "\nColumn names seen: ", colnames(tbl))
		return( NULL)
	}

	obsGenes <- shortGeneName( tbl[[ geneColumn]], keep=1)
	obsFolds <- tbl[[foldColumn]]
	NG <- length( obsGenes)
	cat( "\nN_Genes in table: ", NG)

	# get all the gene sets
	cat( "\nGathering GeneSets and evaluating overlap..")
	if ( ! is.list( geneSets)) {
		tmp <- gatherGeneSets( geneSets)
		geneSets <- tmp$geneSets[[1]]
	}
	gsNames <- names(geneSets) <- cleanGeneSetModuleNames( names(geneSets), wrap=F)
	NGS <- length( geneSets)
	cat( "\nN_GeneSets: ", NGS, "\n")

	# find the gene sets that have the best overlap with these genes
	nOverlap <- avgFold <- pval <- jaccard <- vector( length=NGS)
	for ( i in 1:NGS) {
		thisSet <- geneSets[[i]]
		who <- which( obsGenes %in% thisSet)
		N <- length(who)
		if ( ! N) {
			nOverlap[i] <- 0
			avgFold[i] <- 0
			jaccard[i] <- 0
			pval[i] <- 1
			next
		}
		nOverlap[i] <- N
		avgFold[i] <- mean( obsFolds[who], na.rm=T)
		pval[i] <- if (N>1) t.test( obsFolds[who], mu=0)$p.value else 1
		mySet <- obsGenes[who]
		jaccard[i] <- length(intersect(mySet,thisSet)) / length(union(mySet,thisSet))
		if ( i %% 100 == 0) cat( "\r", i, avgFold[i], pval[i], "  ", substr( gsNames[i], 1, 60))
	}

	out <- data.frame( "GeneSet"=gsNames, "Avg_Log2Fold"=round(avgFold, digits=3), "P_Value"=pval, 
				"N_Genes_Total"=sapply(geneSets,length), "N_Genes_Given"=nOverlap, 
				"Jaccard"=round( jaccard, digits=3), stringsAsFactors=F)

	drops <- which( out$P_Value > max.pvalue.keep)
	if ( length( drops)) out <- out[ -drops, ]
	ord <- diffExpressRankOrder( out$Avg_Log2Fold, out$P_Value)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)
	cat( "  Done.\n")

	return( out)
}

