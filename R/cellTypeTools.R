# cellTypeTools.R -- tools for working with our immune cell subset datasets, as per Gene Sets, etc.
#			we will/have combinded human and mouse data

# circa 2020:   cloning all the LifeCycle data analysis & deconstruction tools to
#				apply them to mammalian immune cell type datasets.  See down below.

# 2022:  Allowing 2+ cell types per call.  Forces many small changes...


# load gene to cell type association data, by looking up the named data object in the Options.txt file
`getCellTypeGeneAssociation` <- function( reference=getCellTypeReference(), optionsFile="Options.txt", 
					speciesID=getCurrentSpecies()) {

	# build the name of the R data object to be loaded
	if ( is.na( reference)) {
		reference <- defaultCellTypeReference( optionsFile=optionsFile, speciesID=speciesID)
	}

	# see if this is what is already loaded
	if ( exists( "Reference", envir=CellTypeEnv)) {
		curReference <- CellTypeEnv[[ "Reference"]]
		if ( ! is.null( curReference) && curReference == reference) {
			if ( exists( "GeneAssociation", envir=CellTypeEnv)) return( CellTypeEnv[[ "GeneAssociation"]])
		}
	}
	
	prefix <- getOtherSpeciesFilePrefix( speciesID)
	referenceName <- paste( prefix, reference, "GeneAssociation", sep=".")

	# try to load that
	geneCellTypes <- NULL

	# allow the file to exist in the current folder, while default is the installed data object
	localFile <- file.path( ".", paste( referenceName, "rda", sep="."))
	if ( file.exists( localFile)) {
		cat( "\n  loading gene association from local file: ", localFile)
		load( localFile, envir=environment())
	} else {
		data( list=referenceName, package="DuffyTools", envir=environment())
	}
	if ( is.null(geneCellTypes)) {
		cat( "\nWarning:  tried to load a 'DuffyTools/data/' installed object file called: ", referenceName)
		cat( "\n          expected it to contain a data object named:  'geneCellTypes' \n")
		return( NULL)
	}

	neededColumns <- c( "GENE_ID", "CellType", "PctExpression")
	if ( ! all( neededColumns %in% colnames(geneCellTypes))) {
		cat( "\nWarning: some required gene cell type columns not found.  Expected: ", neededColumns, "\n")
		return(NULL)
	}
	
	# stash a copy for fast retrieval
	CellTypeEnv[[ "GeneAssociation" ]] <- geneCellTypes

	return( geneCellTypes)
}


# map genes to cell types, using the current cellTypeReference dataset.  Allow returning more than one name.
`gene2CellType` <- function( genes, max.types=1, speciesID=getCurrentSpecies()) {

	out <- rep.int( "", N <- length(genes))

	genesIn <- shortGeneName( genes, keep=1)
	geneCellTypes <- getCellTypeGeneAssociation( reference=getCellTypeReference(), speciesID=speciesID)
	if ( is.null(geneCellTypes)) {
		cat( "\nWarning:  No 'GeneCellTypes' object loaded for species: ", speciesID)
		cat( "\n  Tried to load 'GeneAssociation' data object for: ", getCellTypeReference())
		return( out)
	}

	if ( max.types == 1) {
		where <- match( genesIn, geneCellTypes$GENE_ID, nomatch=0)
		out[ where > 0] <- geneCellTypes$CellType[ where]
	} else {
		# need to find 2+ cell type entries per gene
		# so crop the data to contain just the genes we want, for speed.
		# and map to where in the output list each gene ends up
		smlCT <- subset( geneCellTypes, GENE_ID %in% genesIn)
		if (nrow(smlCT)) {
			smlGeneFac <- factor(smlCT$GENE_ID)
			whereOut <- match( levels(smlGeneFac), genesIn)
			k <- 0
			tapply( 1:nrow(smlCT), smlGeneFac, function(x) {
				if ( length(x) > max.types) x <- x[ 1:max.types]
				cellTypeStr <- paste( smlCT$CellType[x], ":", smlCT$PctExpression[x], "%", sep="", collapse="; ")
				k <<- k + 1
				out[ whereOut[k]] <<- cellTypeStr
				return(NULL)
				})
		}
		# Note: a gene could appear twice or more in the output.  Make sure we fill all locations.
		# because of the way match works, first location always has the full answer
		dupsIn <- which( duplicated( genesIn))
		if ( length( dupsIn)) {
			for ( k in dupsIn) {
				allWho <- which( genesIn == genesIn[k])
				out[ allWho] <- out[ allWho[1]]
			}
		}
	}
	if ( all( out == "")) warning( paste( "No genes mapped to cell types.  Verify current SpeciesID"))

	out
}


`cellTypeEnrichment` <- function( cellTypes, mode=c("genes", "geneSets"), 
				upOnly=TRUE, minEnrich=1.25, maxPvalue=0.01, 
				wt.enrich=1, wt.pvalue=2, speciesID=getCurrentSpecies(), 
				geneUniverse=NULL, correct=TRUE, verbose=T) {
	
	# grab the universe of cell type data
	mode <- match.arg( mode)
	if ( mode == "geneSets") {
		geneSetCellTypes <- getCellTypeGeneSetAssociation()
		if ( is.null(geneSetCellTypes)) {
			cat( "\nWarning:  No 'GeneSetCellTypes' object loaded for species: ", speciesID)
			cat( "\n  Tried to load 'GeneAssociation' data object for: ", getCellTypeReference())
			return(NULL)
		}
		cellTypeUniverse <- geneSetCellTypes$CellType
	}
	if ( mode == "genes") {
		# speciesID is needed, since orthologging may be required
		if (is.null(speciesID)) stop( "'speciesID' must not be NULL for gene cell types.")
		prefix <- getOtherSpeciesFilePrefix( speciesID)
		geneCellTypes <- getCellTypeGeneAssociation( reference=getCellTypeReference(), speciesID=speciesID)
		if ( is.null(geneCellTypes)) {
			cat( "\nWarning:  No 'GeneCellTypes' object loaded for species: ", speciesID)
			cat( "\n  Tried to load 'GeneAssociation' data object for: ", getCellTypeReference())
			return(NULL)
		}

		# this table often has more than one cell type per gene, to show the diversity
		# limit it to the first/best cell type choice for each gene
		# Note:  with allowing 2+ cell types per gene, we need to NOT remove duplicates, as it distorts 
		# the true number of genes per cell type
		# dups <- which( duplicated( geneCellTypes$GENE_ID))
		# if ( length(dups)) geneCellTypes <- geneCellTypes[ -dups, ]

		# allow the caller to shrink the universe, as from arrays with less than all genes present
		if ( ! is.null( geneUniverse)) {
			geneUniverse <- shortGeneName( geneUniverse, keep=1)
			keep <- which( geneCellTypes$GENE_ID %in% geneUniverse)
			if (verbose) cat( "\nReducing Gene Universe from ", nrow(geneCellTypes), " to ", length(keep))
			geneCellTypes <- geneCellTypes[ keep, ]
		}

		# now we know how big the universe of genes is
		cellTypeUniverse <- geneCellTypes$CellType
	}

	# any missing cell types in the input should not influence our counts or our enrichment calls
	# remove them
	drops <- c( which( cellTypes == ""), which( is.na( cellTypes)))
	if ( length( drops)) cellTypes <- cellTypes[ -drops]

	# almost ready to calculate enrichment.
	# in the past, these cell type strings were always just a single cell type.  Now they can instead
	# be the relative percentages of 2+ cell types.  Always of the form:   CellType1:XX%; CellType2:YY%; etc
	maxTermsPerCellType <- 1
	isNewFormat <- any( grepl( "; ", cellTypes))
	if ( isNewFormat) {
		# first we need to know how many cell types any element could contain, so count the terms
		# to know the upper limit
		cellTermLists <- strsplit( cellTypes, split="; ")
		nTermsEach <- sapply( cellTermLists, length)
		maxTermsPerCellType <- max( nTermsEach, na.rm=T)
		# now make it one common pool of terms, and then proportionalize them based on their relative %
		cellTerms <- unlist( cellTermLists)
		cellTermCellNames <- sub( "\\:[0-9]+%", "", cellTerms)
		cellTermPcts <- sub( "(.+\\:)([0-9]+)(%)", "\\2", cellTerms)
		# use the percentages to prorate all the cell terms, and then scale back to how many cell types we were given
		cellTbl <- table( rep( cellTermCellNames, times=as.numeric(cellTermPcts)))
		cellTbl <- round( cellTbl * length(cellTerms) / sum(cellTbl))
		cellTbl <- cellTbl[ cellTbl > 0]
		# given this scaled by proportion list, recreate what the full flat list would have been
		cellTypes <- rep( names(cellTbl), times=as.numeric(cellTbl))
		givenTbl <- table( cellTypes)
	} else {
		# old format is just a single cell type per gene
		givenTbl <- table( cellTypes)
	}
	# now we must do the same to the cell universe side.  These used to have exactly one entry each,
	# but the new format allows 2+ entries each.  Trim the universe to have no more than that many each too.
	# By definition, each element of the universe is sorted by decreasing percentage, so use the first K each.
	if ( mode == "genes") {
		use <- unlist( tapply( 1:nrow(geneCellTypes), factor(geneCellTypes$GENE_ID), function(x) x[1:maxTermsPerCellType]))
		geneCellTypes <- geneCellTypes[ use, ]
		cellTypeUniverse <- geneCellTypes$CellType
	} else {
		use <- unlist( tapply( 1:nrow(geneSetCellTypes), factor(geneSetCellTypes$GeneSetName), function(x) x[1:maxTermsPerCellType]))
		geneSetCellTypes <- geneSetCellTypes[ use, ]
		cellTypeUniverse <- geneSetCellTypes$CellType
	}

	# now fully ready to calculate enrichment.
	givenPcts <- givenTbl * 100 / sum(givenTbl)
	givenNames <- names( givenTbl)
	allTbl <- table( cellTypeUniverse)
	allPcts <- allTbl * 100 / sum(allTbl)
	allNames <- names( allTbl)
	allExpect <- allPcts * length(cellTypes) / 100
	
	# we will only consider the names in the Universe
	both <- allNames
	# gather the counts and percentages for both worlds
	whereG <- match( both, givenNames, nomatch=NA)
	Gcnt <- ifelse( is.na(whereG), 0, givenTbl[whereG])
	Gpct <- ifelse( is.na(whereG), 0, givenPcts[whereG])
	whereA <- match( both, allNames, nomatch=0)
	Acnt <- ifelse( is.na(whereA), 0, allTbl[whereA])
	Apct <- ifelse( is.na(whereA), 0, allPcts[whereA])
	ExpCnt <- ifelse( is.na(whereA), 0, allExpect[whereA])

	# enrichment is how much more in our set than the entire table
	enrich <- Gpct / Apct
	
	# the probabilities are for each cell type,  calculate them
	N_All <- length( cellTypeUniverse)
	N_Given <- length( cellTypes)
	pvals <- vector( length=length(both))
	if (verbose) cat( "\nCalculating enrichment P-values..\n")
	for ( i in 1:length(pvals)) {
		x <- Gcnt[i]
		m <- Acnt[i]
		n <- N_All - m
		k <- N_Given
		lowTailWanted <- (enrich[i] < 1.0)
		# get the entire prob. dist. and sum up the half we want...
		allPs <- dhyper( 0:k, m, n, k)
		if ( lowTailWanted) {
			# the probability of 0 to X genes
			pvals[i] <- sum( allPs[1:(x+1)], na.rm=T)
		} else {
			# the probability of X up to K genes
			pvals[i] <- sum( allPs[(x+1):length(allPs)], na.rm=T)
		}
		if (verbose) cat( "\n", i, both[i], enrich[i], pvals[i])
	}
	pvals[ is.na(pvals)] <- 1
	if (verbose) cat( "\n")

	# do correction for multiple comparisons?
	if (correct) pvals <- p.adjust( pvals, method="BH")

	enrich <- round( enrich, digits=3)
	Apct <- round( Apct, digits=2)
	Gpct <- round( Gpct, digits=2)
	ExpCnt <- round( ExpCnt, digits=1)
	signif <- rep.int( "", length(pvals))
	signif[ pvals < 0.1] <- "."
	signif[ pvals < 0.05] <- "*"
	signif[ pvals < 0.01] <- "**"
	signif[ pvals < 0.001] <- "***"
	pvals <- as.numeric( formatC( pvals, format="e", digits=2))

	out <- data.frame( both, enrich, pvals, Acnt, Apct, Gcnt, Gpct, ExpCnt, stringsAsFactors=FALSE)
	colnames( out) <- c( "CellType", "Enrichment", "P_Value", "N_Total", "Pct_Total", "N_Given", 
				"Pct_Given", "N_Expect")
	out$OverUnder <- ifelse( enrich >= 1, "Over", "Under")
	out$OverUnder[ enrich >= 0.99 & enrich <= 1.01] <- ""
	out$Signif <- signif
	
	# having a gene count of zero forces a zero enrichment, which is hard to rank
	# use a small fudge factor to prevent exact zeros
	enrichDE <- enrich
	isZero <- which( Gcnt == 0)
	if ( length( isZero)) {
		smallCnt <- 1 / N_Given
		Gpct[ isZero] <- smallCnt / N_Given
		enrichDE[isZero] <- Gpct[isZero] / Apct[isZero]
	}

	# sort by enrichment and P-value, where 1.0 means not enriched...
	ord <- diffExpressRankOrder( enrichDE, out$P_Value, wt.fold=wt.enrich, wt.pvalue=wt.pvalue,
					notDE.value=1.0)
	out <- out[ ord, ]
	rownames( out) <- 1:nrow( out)
	
	# now we have all data, send back just what was asked for
	testValue <- out$Enrichment
	testP <- out$P_Value
	who <- which( testValue >= minEnrich & testP <= maxPvalue)
	if ( ! upOnly) {
		who2 <- which( testValue <= (1/minEnrich) & testP <= maxPvalue)
		who <- sort( unique( c( who, who2)))
	}
	out <- out[ who, ]
	if ( nrow(out) > 0) rownames( out) <- 1:nrow( out)
	return( out)
}


`cellTypePies` <- function( fileSet, radius=0.9, 
			pie.main=sub( "\\..+", "", basename(fileSet)), ...) {

	# given a set of filenames of "CellTypeEnrichment.csv" files, turn them into pies
	pie.main <- rep( pie.main, length.out=length(fileSet))
	fileFound <- file.exists( fileSet)
	N <- sum( fileFound)
	fileSet <- fileSet[ fileFound]
	pie.main <- pie.main[ fileFound]

	saveMAI <- par( "mai")
	par( mai=c(0.2, 0.1,0.4,0.1))
	par( mfrow=c(1,1))
	if (N > 1) par( mfrow=c(1,2))
	if (N > 2) par( mfrow=c(2,2))
	if (N > 4) par( mfrow=c(2,3))
	if (N > 6) par( mfrow=c(3,3))
	if (N > 9) par( mfrow=c(3,4))


	for ( i in 1:N) {
		tbl <- read.csv( fileSet[i], as.is=T)
		myCells <- tbl$CellType
		myCounts <- tbl$N_Given
		if ( is.null( myCells)) next
		if ( is.null( myCounts)) next
		ord <- order( myCells)
		myCells <- myCells[ord]
		myCounts <- myCounts[ord]
		NC <- length(myCells)
		# 'Whole Blood' is last...
		myColors <- rainbow( NC, end=0.85)
		myColors[NC] <- 'brown'

		pie( myCounts, labels=myCells, clockwise=T, col=myColors, main=pie.main[i], radius=radius, ...)
	}
	par( mfrow=c(1,1))
	par( mai=saveMAI)
}


`cellTypeExpression` <- function( m, mode=c("sqrtmean","log2"), min.expression=0, verbose=TRUE) {

	# given a matrix of gene expression, reduce it to a matrix of cell type expression
	if ( is.null( rownames(m))) stop( "Matrix must have genes as rownames.")

	# drop genes with very low expression?
	if ( min.expression > 0) {
		if (verbose) cat( "\nDropping gene rows with expression < ", min.expression)
		bigV <- apply( m, 1, max, na.rm=T)
		drops <- which( bigV < min.expression)
		if ( length( drops)) {
			m <- m[ -drops, ]
			if (verbose) cat( "\nN_Genes dropped: ", length(drops))
		}
	}
	
	gids <- shortGeneName( rownames(m), keep=1)
	cids <- gene2CellType( gids)

	# now that cell type calls can be 2+ types, may need to accumulate differently
	isNewFormat <- any( grepl( "; ", cids))
	if ( isNewFormat) {
		cellTerms <- strsplit( cids, split="; ")
		cellsPerGene <- sapply( cellTerms, length)
		cellTermsV <- unlist( cellTerms)
		cellTermCellNames <- sub( "\\:[0-9]+%", "", cellTermsV)
		cellTermPcts <- as.numeric( sub( "(.+\\:)([0-9]+)(%)", "\\2", cellTermsV))
		cellTermGenes <- rep( gids, times=cellsPerGene)
		ctGenePtrs <- match( cellTermGenes, gids)
		NCTterms <- length(cellTermCellNames)
		cellTypes <- sort( unique( cellTermCellNames))
	} else {
		cellTypes <- sort( unique( cids))
	}

	# create the storage to hold the result
	NCT <- length( cellTypes)
	if ( NCT < 2) {
		cat( "\nUnexpected low number of cell types: ", NCT)
		cat( "\nDebug:  Gene IDs rownames(m): ", head( gids))
		cat( "\nDebug:  Cell Type calls:      ", head( cids))
		return( NULL)
	}
	out <- matrix( NA, nrow=NCT, ncol=ncol(m))
	colnames(out) <- colnames(m)
	rownames(out) <- gsub( "\\.\\.+", ".", make.names( cellTypes))

	# now fill the new matrix with the average expression by cell type
	# Old cell types of exactly 1 type per gene are easy.  New style needs prorating
	if ( isNewFormat) {
		cellFac <- factor( cellTypeCellNames)
		mode <- match.arg( mode)
		for ( i in 1:ncol(m)) {
			geneV <- m[ , i]
			if ( mode == "sqrtmean") {
				cellV <- tapply( 1:NCTterms, cellFac, function(x) {
						myGeneV <- geneV[ ctGenePtrs[x]]
						myPcts <- cellTermPcts[x]
						# use the 1-100 pcts as relative votes for each gene
						myV <- rep( myGeneV, times=myPcts)
						sqrtmean( myV, na.rm=T)})
			} else {
				cellV <- tapply( 1:NCTterms, cellFac, function(x) {
						myGeneV <- geneV[ ctGenePtrs[x]]
						myPcts <- cellTermPcts[x]
						# use the 1-100 pcts as relative votes for each gene
						myV <- rep( myGeneV, times=myPcts)
						logmean( myV+1, na.rm=T)})
			}
			out[ , i] <- cellV
		}
	} else {
		cellFac <- factor( cids)
		mode <- match.arg( mode)
		for ( i in 1:ncol(m)) {
			geneV <- m[ , i]
			if ( mode == "sqrtmean") {
				cellV <- tapply( geneV, cellFac, sqrtmean, na.rm=T)
			} else {
				cellV <- tapply( geneV, cellFac, function(x) logmean( x+1, na.rm=T))
			}
			out[ , i] <- cellV
		}
	}

	drops <- which( rownames(out) %in% c( "", " ", "X"))
	if ( length( drops)) out <- out[ -drops, ]

	return( out)
}


`cellTypeHeatmap` <- function( m, nColors=50, min.expression=0, verbose=FALSE, ...) {

	# given a matrix of gene expression, make a heatmap that reduces it to a matrix of 
	# cell types by samples

	# step 1:  reduce the matrix of genes to a matrix of cell type expression values
	cellM <- cellTypeExpression( m, min.expression=min.expression, verbose=verbose)
	cellIDs <- rownames(cellM)

	heatColors <- heatMapColors( nColors, palette="red-white-blue", rampExponent=1, plotRamp=F)

	# step 2:  turn that into heatmap data, and show it
	# convert  expression to M-values
	cellM2 <- expressionMatrixToMvalue( cellM)

	# and pass it to the heatmap tool
	ans <- heatmap( cellM2, heatColors=heatColors, ...)
	return( invisible( ans))
}


`cellTypeGeneHeatmap` <- function( m, nColors=50, nGenesPerCellType=c("25","50","100","250","500","1000"), 
					excludeCellTypePattern=NULL, includeCellTypePattern=NULL, 
					speciesID=getCurrentSpecies(), verbose=T, ...) {

	# given a matrix of gene expression, keep only the top N genes from the selected cell types,
	# and pass that subset on to the heatmap tool.
	if ( is.null( rownames(m))) stop( "Matrix must have genes as rownames.")
	rownames(m) <- shortGeneName( rownames(m), keep=1)

	# step 1:  get the cell type data, and down-select the gene sets to use
	CellTypeSetup()
	reference <- getCellTypeReference()
	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	dataSetName <- paste( getCurrentSpeciesFilePrefix(), reference, "AllGeneSets", sep=".")
	allGeneSets <- NULL
	data( list=dataSetName, envir=environment())
	if ( is.null( allGeneSets)) {
		cat( "\nFailed to load needed cell type dataset: ", dataSetName, "\nUnable to plot...")
		return(NULL)
	}
	geneSetNames <- names(allGeneSets)

	# keep just the right sized sets
	nGenesPerCellType <- match.arg( nGenesPerCellType)
	keep <- grep( paste( "top", as.character(nGenesPerCellType), "genes", sep=" "), geneSetNames)
	geneSetNames <- geneSetNames[ keep]
	# then the gene set name include/exclude patterns
	if ( ! is.null( excludeCellTypePattern)) {
		drops <- grep( excludeCellTypePattern, geneSetNames)
		if ( length(drops)) geneSetNames <- geneSetNames[ -drops]
	}
	if ( ! is.null( includeCellTypePattern)) {
		keep <- grep( includeCellTypePattern, geneSetNames)
		geneSetNames <- geneSetNames[ keep]
	}
	NGS <- length(geneSetNames)
	if (verbose) cat( "\nFinal set of Cell Type Gene Sets: ", NGS, "\n", geneSetNames)
	if ( ! NGS) {
		cat( "\nNo Cell Type subsets selected..")
		return(NULL)
	}

	# find the genes for each of these gene sets
	geneSetGenes <- vector( mode="list", length=NGS)
	for ( i in 1:NGS) {
		whereGS <- match( geneSetNames[i], names(allGeneSets))
		myGenes <- allGeneSets[[ whereGS]]
		whereM <- match( myGenes, rownames(m), nomatch=0)
		geneSetGenes[[i]] <- whereM[ whereM > 0]
		# make sure we see genes...
		if ( all( whereM == 0)) cat( "\nWarning: no genes found for: ", geneSetNames[i], "  Check current species.")
	}
	nGenesPerSet <- sapply( geneSetGenes, length)

	# step 2:  with the down seletion done, ready to make the matrix of genes we want
	allGenePtrs <- unlist( geneSetGenes)
	mGeneSetGenes <- m[ allGenePtrs, ]

	# step 3:  get ready to send it to the heatmap tool
	heatColors <- heatMapColors( nColors, palette="red-white-blue", rampExponent=1, plotRamp=F)

	# convert gene expression to M-values
	M2 <- expressionMatrixToMvalue( mGeneSetGenes, verbose=verbose)

	# create some side bar coloring
	cellColors <- rainbow( NGS, end=0.85)
	rowColors <- rep( cellColors, times=nGenesPerSet)
	#cellColors <- rainbow( length(allGeneSets), end=0.85)
	#rowColors <- rep( cellColors[ match(geneSetNames, names(allGeneSets))], times=nGenesPerSet)
	# try to put each cell type name in the center of that band
	rowNames <- rep.int( NA, nrow(M2))
	useGeneSetNames <- sub( ": top.+", "", geneSetNames)
	yBig <- cumsum( nGenesPerSet)
	yHalf <- round( nGenesPerSet/2)
	rowNames[ yBig - yHalf] <- useGeneSetNames
	# and put them in alphabetical order, top to bottom
	M2 <- M2[ rev( 1:nrow(M2)), ]
	rowColors <- rev( rowColors)
	rowNames <- rev( rowNames)

	# and pass it to the heatmap tool
	ans <- heatmap( M2, heatColors=heatColors, Rowv=NA, RowSideColors=rowColors, 
			labRow=rowNames, cexRow=(0.2+1/log10(max(NGS,3))), verbose=verbose, ...)
	return( invisible( ans))
}


# porting of Life Cycle tools to be used as Cell Type tools below here...

`CellTypeSetup` <- function( reference=getCellTypeReference(), optionsFile="Options.txt", reload=FALSE) {

	# get the name of the R data object all ready loaded or to be loaded
	speciesID <- getCurrentSpecies()
	if ( is.na( reference)) {
		reference <- defaultCellTypeReference( optionsFile=optionsFile, speciesID=speciesID)
	}

	# compare this reference name against what is loaded already
	isReady <- ( exists( "VectorSpace", envir=CellTypeEnv) && !reload)
	isRightSpecies <- isRightReference <- FALSE
	if ( isReady) {
		curSpecies <- get( "Species", envir=CellTypeEnv)
		if( ! is.null( curSpecies)) isRightSpecies <- ( curSpecies == getCurrentSpecies())
		curReference <- get( "Reference", envir=CellTypeEnv)
		if( ! is.null( curReference)) isRightReference <- ( curReference == reference)
	}
	if ( !isReady || !isRightSpecies || !isRightReference) {
		# load the new target matrix
		loadCellTypeMatrix( reference, unitVectorMode="absolute")
		# since it changed, remove any items that relate to previous
		if ( exists( "GeneAssociation", envir=CellTypeEnv)) rm( "GeneAssociation", envir=CellTypeEnv) 
		if ( exists( "GeneSetAssociation", envir=CellTypeEnv)) rm( "GeneSetAssociation", envir=CellTypeEnv) 
	}
	return()
}


`loadCellTypeMatrix` <- function( reference, unitVectorMode=c("absolute","relative","none","squared","cubed"), 
				custom.file=NULL, custom.colors=NULL, preNormalize=TRUE, postNormalize=TRUE, 
				doRMA=FALSE, min.value=0, min.spread=2.0, verbose=FALSE) {

	unitVectorMode <- match.arg( unitVectorMode)
	ans <- NULL

	cat( "\nSetting up CellType dataset:  ", reference)

	# get the data from the package
	if ( tolower(reference) != "custom") {
	
		# turn the name into a DuffyTools data object
		prefix <- getCurrentSpeciesFilePrefix()
		referenceName <- paste( prefix, reference, "TargetMatrix", sep=".")
		targetM <- targetColors <- NULL

		# allow the file to exist in the current folder, while default is the installed data object
		localFile <- file.path( ".", paste( referenceName, "rda", sep="."))
		if ( file.exists( localFile)) {
			cat( "\n  loading target matrix from local file: ", localFile)
			load( localFile, envir=environment())
		} else {
			data( list=referenceName, package="DuffyTools", envir=environment())
		}
		if ( is.null( targetM)) {
			cat( "\nFailed to load Cell Type target matrix data object: ", referenceName)
			cat( "\n  expected a data object containing 'targetM'")
			return(NULL)
		}
		if ( ! exists( "targetColors")) targetColors <- rainbow( ncol(targetM), end=0.8)

		ans <- buildCellTypeVectorSpace( file=NULL, tbl=targetM, 
			unitVectorMode=unitVectorMode, min.value=min.value, min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=FALSE, verbose=verbose)
		ans2 <- buildCellTypeVectorSpace( file=NULL, tbl=targetM, 
			unitVectorMode="none", min.value=min.value, min.spread=NULL, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=FALSE, verbose=FALSE)
	}

	if ( tolower(reference) == "custom") {

		# pre read the file to make the dimensions construct
		if ( ! file.exists( custom.file)) stop( paste( "Custom Cell Type Target Matrix file not found: ", custom.file))
		tmp <- read.delim(custom.file, as.is=T)
		# should verify it has rownames that are the geneIDs....

		cat( "\nBuilding custom Cell Type dataset:   Dimensions:\n")
		print( colnames(tmp))
		cat( "\nFirst genes: \n")
		print( head( tmp))

		ans <- buildCellTypeVectorSpace( tbl=tmp, 
			unitVectorMode=unitVectorMode, min.value=min.value, min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=verbose)
		ans2 <- buildCellTypeVectorSpace( tbl=tmp,
			unitVectorMode="none", min.value=min.value, min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=doRMA, verbose=verbose)
			
		targetColors <- custom.colors
	}

	# no matter what, force these 2 views to have identical genes
	finalGenes <- intersect( ans$GENE_ID, ans2$GENE_ID)
	ans <- ans[ match( finalGenes, ans$GENE_ID), ]
	ans2 <- ans2[ match( finalGenes, ans2$GENE_ID), ]

	CellTypeEnv[[ "VectorSpace" ]] <- ans
	CellTypeEnv[[ "IntensitySpace" ]] <- ans2
	CellTypeEnv[[ "Species" ]] <- getCurrentSpecies()
	CellTypeEnv[[ "Reference" ]] <- reference
	if ( ! is.null( ans)) {
		CellTypeEnv[[ "N_STAGES" ]] <- ncol(ans) - 2   # geneID, Product come first
		CellTypeEnv[[ "STAGE_NAMES" ]] <- colnames(ans)[3:ncol(ans)]
		if ( ! is.null( targetColors)) CellTypeEnv[[ "STAGE_COLORS"]] <- targetColors
	}
	cat( "\nDone.\n")
}


`reshapeCellTypeMatrix` <- function( f, geneColumn="GENE_ID", intensityColumn="RPKM_M", verbose=TRUE) {

	# allow the reference cell type matrix to alter its distribution to match a given transcriptome
	
	# read in the given transcriptome, that we will use as the new distribution
	if ( is.character(f)) {
		tbl <- read.delim(f, as.is=T)
		cat( "\nReshaping CellType data to mimic: ", basename(f))
	} else if ( is.data.frame(f)) {
		tbl <- f
	}
	if ( ! all( c(geneColumn, intensityColumn) %in% colnames(tbl))) stop( "Required columns not found in transcriptome file")
	givenGenes <- shortGeneName( tbl[[ geneColumn]], keep=1)
	givenInten <- as.numeric( tbl[[ intensityColumn]])
	
	# get the current cell type matrix
	cellDF <- CellTypeEnv[[ "IntensitySpace" ]]
	if ( is.null(cellDF)) stop( "run 'CellTypeSetup()' first, before reshaping the Cell Type Matrix")
	cellM <- as.matrix( cellDF[ , 3:ncol(cellDF)])	# GeneID, Product, then the cell types...
	cellGenes <- shortGeneName( cellDF$GENE_ID, keep=1)
	
	# about to use REI (quantile normalization) to do the reshaping
	# only use the given genes that are in the matrix already
	use <- match( cellGenes, givenGenes, nomatch=NA)
	cat( "\n  Doing Rank Equivalent Intensity normalization..")
	newInten <- rankEquivIntensity( cellM, targetIntensity=givenInten[use], blendMode="targetOnly")
	
	# we cannot let the reshaping turn 'ON' genes that were completely off originally.  Check and prevent
	newInten[ cellM <= 0] <- 0
	
	# now restore the global normalization magnitude that we started with
	globalSum <- median( apply( cellM, 2, sum, na.rm=T))
	for ( i in 1:ncol(newInten)) newInten[ , i] <- newInten[ , i] * globalSum / sum( newInten[ , i], na.rm=T)
	
	# with these modified expression profiles done, now we can redo the unit vectors
	newVect <- newInten
	for( i in 1:nrow( newInten)) {
		oneV <- newInten[ i, ]
		newVect[i, ] <- oneV / sum(oneV, na.rm=T)
	}
	
	# put these reshaped matrices back in place
	ans <- cbind( cellDF[ ,1:2], newVect, stringsAsFactors=F)
	ans2 <- cbind( cellDF[ ,1:2], round(newInten,digits=5), stringsAsFactors=F)
	CellTypeEnv[[ "VectorSpace" ]] <- ans
	CellTypeEnv[[ "IntensitySpace" ]] <- ans2
	return(NULL)
}


# always do the setup check, so we force the cell type target matrix into existence
`getCellTypeMatrix` <- function( mode=c("IntensitySpace", "VectorSpace")) {

	CellTypeSetup()
	mode <- match.arg( mode)

	tbl <- CellTypeEnv[[ mode]]
	if ( is.null(tbl)) {
		cat( "\n  Warning: no cell type matrix available for current species and cell type reference name.")
		return(NULL)
	}
	genes <- tbl$GENE_ID

	# get just the values, not the IDs and Products
	m <- as.matrix( tbl[ ,3:ncol(tbl)])
	rownames(m) <- genes
	return(m)
}


# change this function to be a 'query only' mode.  Just return what is already loaded
`getCellTypeColors` <- function() {

	if ( ! exists( "VectorSpace", envir=CellTypeEnv)) return(NULL)
	colors <- CellTypeEnv[[ "STAGE_COLORS"]]
	names(colors) <- CellTypeEnv[[ "STAGE_NAMES"]]
	if ( is.null( colors)) {
		tbl <- CellTypeEnv[[ "VectorSpace"]]
		nCellTypes <- ncol(tbl) - 2 	# GeneID, Product,  then data....
		colors <- rainbow( nCellTypes, end=0.8)
	}
	return( colors)
}


# make this function be a 'query only' mode.  Use it to ask for the currently defined reference
`getCellTypeReference` <- function() {

	if ( ! exists( "VectorSpace", envir=CellTypeEnv)) return(NA)
	refSpecies <- CellTypeEnv[[ "Species" ]] 
	if ( refSpecies != getCurrentSpecies()) return(NA)
	reference <- CellTypeEnv[[ "Reference"]]
	return( reference)
}


# get the default reference to use for the given species
`defaultCellTypeReference` <- function( optionsFile="Options.txt", speciesID=getCurrentSpecies()) {

	defaultReference <- NA
	if (speciesID %in% PARASITE_SPECIES) defaultReference <- "Parasite.LifeCycle"
	if (speciesID %in% MAMMAL_SPECIES) defaultReference <- "27.Blood.CellTypes"
	if (speciesID %in% BACTERIA_SPECIES) defaultReference <- "TB.Phenotypes"
	if ( file.exists( optionsFile)) {
		reference <- getOptionValue( optionsFile, arg="cellTypeReference", speciesID=speciesID, 
					notfound=defaultReference, verbose=FALSE)
	} else {
		reference <- defaultReference
	}
	return( reference)
}


`buildCellTypeVectorSpace` <- function( file="", tbl=NULL,
					unitVectorMode=c( "absolute","relative","none","squared","cubed"), 
					min.value=0, min.spread=2, preNormalize=FALSE, 
					postNormalize=TRUE, doRMA=TRUE, verbose=FALSE) {

	unitVectorMode <- match.arg( unitVectorMode)

	# get the original reference set of gene data
	if ( is.null( tbl)) {
		tbl <- read.delim( file, as.is=TRUE)
	}

	# gene IDs are the rownames, and dimension names are the column names
	geneIDs <- shortGeneName( rownames(tbl), keep=1)
	cellIDs <- colnames(tbl)

	# no nned for any orthologging, each matrix is organism specific
	thisSpecies <- getCurrentSpecies()

	# extract just the wanted columns, impose a minimun intensity, and perhaps normalize the intensity columns
	mIn <- as.matrix( tbl)
	colnames( mIn) <- cellIDs
	mIn[ is.na(mIn)] <- min.value
	mIn[ mIn < min.value] <- min.value

	# with microarrays, RMA is more appropriate, otherwise just plain quantile normalize
	if ( preNormalize) {
		if ( doRMA) {
			mIn <- duffyRMA( mIn, verbose=verbose)
		} else {
			mIn <- duffyMagnitudeNormalize( mIn, verbose=verbose)
		}
	}

	# there may be duplicate rows because of gene re-annotation or orthollogging...
	# combine those by mean average
	uGenes <- unique.default( geneIDs)
	nUnique <- length( uGenes)
	if ( nUnique < nrow(mIn)) {
		mNew <- matrix( 0, nrow=nUnique, ncol=ncol(mIn))
		colnames(mNew) <- colnames(mIn)
		for ( i in 1:nUnique) {
			who <- which( geneIDs == uGenes[i])
			if ( length(who) == 1) {
				mNew[ i, ] <- mIn[ who, ]
			} else {
				mNew[ i, ] <- apply( mIn[ who, ], MARGIN=2, FUN=sqrtmean)
			}
		}
	} else {
		mNew <- mIn
	}

	# with microarrays, RMA is more appropriate, otherwise just plain quantile normalize
	if ( postNormalize) {
		if ( doRMA) {
			mNew <- duffyRMA( mNew, verbose=verbose)
		} else {
			mNew <- duffyMagnitudeNormalize( mNew, verbose=verbose)
		}
	}
	
	# now normalize each row to sum to unit vector (length=1)
	if ( unitVectorMode != "none") {
	    for( i in 1:nrow( mNew)) {
		oneV <- mNew[ i, ]
		# subtract the min so the proportions have at least one 0%
		if ( unitVectorMode == "relative") oneV <- oneV - min(oneV)
		oneV <- oneV / sum(oneV)
		mNew[i, ] <- oneV
	    }
	    # explore allowing a change of relative shape, by applying a power correction
	    if ( unitVectorMode %in% c( "squared","cubed")) {
	    	for( i in 1:nrow( mNew)) {
			oneV <- mNew[ i, ]
			oneV <- if ( unitVectorMode == "squared") oneV^2 else oneV^3
			oneV <- oneV / sum(oneV)
			mNew[i, ] <- oneV
	        }
	    }
	}

	if ( ! is.null( min.spread)) {
		min.spread <- as.numeric( min.spread)
		mins <- apply( mNew, MARGIN=1, min, na.rm=T)
		maxs <- apply( mNew, MARGIN=1, max, na.rm=T)
		spreads <- maxs / mins
		drops <- which( spreads < min.spread)
		if ( length(drops) > 0) {
			mNew <- mNew[ -drops, ]
			uGenes <- uGenes[ -drops]
			cat( "\nDropping Genes with too little intensity change between dimensions: ", length(drops))
		}
	}

	ans <- data.frame( "GENE_ID"=uGenes, "PRODUCT"=gene2Product( uGenes), mNew,
			stringsAsFactors=FALSE)
	rownames(ans) <- 1:nrow(ans)
	return( ans)
}
 

`calcCellTypeProfile` <- function( genes, inten) {

	CellTypeSetup()

	# get the cell type data...  there is 'geneID, product', and then the intensities...
	vectorSpace <- CellTypeEnv[[ "VectorSpace"]]
	unitVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build up a tally of how much intensity goes to each cell type
	# do it for all genes, and the make stage histogram from the marker genes
	genes <- shortGeneName( genes, keep=1)

	# the Vector Space table is now in the current species, no need to ortholog here!
	where <- base::match( genes, vectorSpace$GENE_ID, nomatch=0)

	# build the matrix of stage fractions for every gene we were given.
	# switching to only use genes that are defined... others will not contribute
	myFractions <- matrix( 0, nrow=length(genes), ncol=N_STAGES)
	for( i in 1:N_STAGES) {
		myFractions[ where > 0, i] <- unitVectors[ where, i]
	}

	# account for background intensity to put microarray and RNA, etc., on equal footing .. now done elsewhere.!
	useInten <- inten

	# spread these intensitys over those stages
	myVectors <- myIntens <- matrix( 0, nrow=length(genes), ncol=N_STAGES)
	rownames(myVectors) <- rownames(myIntens) <- genes
	colnames(myVectors) <- colnames(myIntens) <- STAGE_NAMES
	for( i in 1:N_STAGES) {
		myVectors[ , i] <- useInten * myFractions[ , i]
		myIntens[ , i] <- inten * myFractions[ , i]
	}

	# now do a summation by stage, and express as percentages...
	allSums <- apply( myVectors, MARGIN=2, FUN=sum, na.rm=T)
	bigSum <- sum( allSums)
	# prevent divide by zero
	ans <- allSums * 100 / max( bigSum, 1, na.rm=T)

	out <- list( "Profile"=ans, "IntensityVectors"=myIntens)
	return( out)
}


`calcCellTypeProfileFromFile` <- function( f, geneColumn="GENE_ID", intensityColumn="RPKM_M", sep="\t") {

	CellTypeSetup()

	# open that file and find the needed columns
	tmp <- read.delim( f, as.is=T, sep=sep)
	if ( nrow( tmp) < 1) return(NULL)

	gset <- tmp[[ geneColumn]]
	if ( is.null(gset)) {
		cat( "calcCellTypeProfile:  gene name column not found.\nFile: ",f,
			"\nTried: ", geneColumn)
		return( NULL)
	}

	inten <- tmp[[ intensityColumn]]
	if ( is.null(inten)) {
		cat( "calcCellTypeProfile:  gene intensity column not found.\nFile: ",f,
			"\nTried: ", intensityColumn)
		return( NULL)
	}

	return( calcCellTypeProfile( gset, inten))
}


`plotCellTypeProfileFromFileSet` <- function( fnames, fids, fcolors=NULL, geneColumn="GENE_ID", 
		intensityColumn="RPKM_M", yMax=NULL, legend.cex=0.8, max.labels=20, 
		mask.low.pct=NULL, mask.low.diff=NULL, label="your label goes here...", sep="\t") {

	CellTypeSetup()

	# make sure we can read those files
	filesOK <- file.exists( fnames)
	if ( !all( filesOK)) {
		cat( "\nCellTypeProfile:  Some transcript files not found:\n")
		print( fnames[ !filesOK])
		return(NULL)
	}

	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build the storage
	nFiles <- length( fnames)
	m <- matrix( nrow=nFiles, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids

	# load each file in turn
	cat( "\nLoading:  ")
	for( i in 1:nFiles) {
		cat( " ", basename(fnames[i]))
		ans <- calcCellTypeProfileFromFile( fnames[i], geneColumn=geneColumn,
				intensityColumn=intensityColumn, sep=sep)
		m[ i, ] <- ans$Profile
	}
	cat( "\n")

	# plot it
	plotCellTypeProfiles(m, col=fcolors, label=label, yMax=yMax, 
				legend.cex=legend.cex, max.labels=max.labels, mask.low.pct=mask.low.pct,
				mask.low.diff=mask.low.diff)

	return( t(m))
}


`plotCellTypeProfileFromMatrix` <- function( geneSet, intenMatrix, fids=colnames(intenMatrix), 
			fcolors=NULL, yMax=NULL, legend.cex=0.8, max.labels=20, 
			mask.low.pct=NULL, mask.low.diff=NULL, label="your label goes here...") {

	CellTypeSetup()

	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build the storage
	nColumns <- ncol( intenMatrix)
	m <- matrix( nrow=nColumns, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids

	# load each dataset in turn
	for( i in 1:nColumns) {
		ans <- calcCellTypeProfile( geneSet, intenMatrix[ ,i])
		m[ i, ] <- ans$Profile
	}

	# plot it
	plotCellTypeProfiles(m, col=fcolors, label=label, yMax=yMax, 
				legend.cex=legend.cex, max.labels=max.labels, mask.low.pct=mask.low.pct,
				mask.low.diff=mask.low.diff)

	return( t(m))
}


`plotCellTypeProfiles` <- function( m, col=NULL, yMax=NULL, label="", 
				legend.cex=0.8, max.labels=20, mask.low.pct=NULL, mask.low.diff=NULL) {

	N <- nrow(m)
	NC <- ncol(m)

	if ( all( is.na(m))) {
		cat( "\nWarning:  no non-zero gene intensities !!")
		cat( "\nPerhaps expression data does not match current species...")
		return(NULL)
	}
	
	# allow masking of low percentage cell types to better use the plot region
	if ( ! is.null( mask.low.pct)) {
		cellMaxes <- apply( m, 2, max, na.rm=T)
		drops <- which( cellMaxes < as.numeric( mask.low.pct))
		if ( length( drops)) {
		 	m <- m[ , -drops, drop=F]
		 	NC <- ncol(m)
		 	if (NC < 2) {
		 		cat( "\nMasked away too many cell types..")
		 		return( NULL)
		 	}
		}
	}
	if ( ! is.null( mask.low.diff)) {
		cellDiffs <- apply( m, 2, function(x) diff( range( x, na.rm=T)))
		drops <- which( cellDiffs < as.numeric( mask.low.diff))
		if ( length( drops)) {
		 	m <- m[ , -drops, drop=F]
		 	NC <- ncol(m)
		 	if (NC < 2) {
		 		cat( "\nMasked away too many cell types..")
		 		return( NULL)
		 	}
		}
	}

	par( "mai"=c( 1.5, 0.95, 0.85, 0.2))
	las <- 3
	border <- par( "fg")
	if ( is.null( col)) {
		col=gray( seq( 0.2, 1.0, length.out=N))
	} else {
		if ( N*NC >= 200) border <- NA
	}

	barSpace <- c( 0, N/4)

	if ( is.null( yMax)) yMax <- max(m, na.rm=T) * 1.15
	mainText <- paste( "Cell Type Expression Profile Plot:\n", label)

	mp <- barplot(m, beside=T, col=col, border=border, main=mainText, 
		ylab="Percent of Total Gene Intensity", 
		space=barSpace, las=las, font.lab=2, font.axis=2, cex.lab=1, cex.axis=1, cex.names=0.8,
		xaxs="i", ylim=c(0,yMax), xlim=c( -1, (N*1.2)*(ncol(m)*1.2)+2))
	
	# limit the legend to a reasonable number
	who <- 1:N
	isSubset <- FALSE
	if ( N > 12) legend.cex <- legend.cex * 0.95
	if ( N > 18) legend.cex <- legend.cex * 0.95
	if ( N > max.labels) {
		who <- seq.int( 1, N, length.out=max.labels)
		isSubset <- TRUE
		legend.cex <- legend.cex * 0.95
	}
	ans <- legend( "topright", rownames(m)[who], fill=col[who], cex=legend.cex, bg="white")
	if ( isSubset) mtext( "not all labeled", side=3, adj=1, cex=legend.cex*0.95)
	dev.flush()
	return( invisible( mp))
}


`plotCellTypeProfileUnitVectors` <- function( gSet, col=1, lwd=1, legend=NA, plot=TRUE, yMax=1,
				legend.cex=1, label="", annotate=FALSE) {

	CellTypeSetup()

	# get the life cycle data...  there is 'geneID, product', and then the intensities...
	vectorSpace <- CellTypeEnv[[ "VectorSpace"]]
	unitVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	if (plot) {

		par( "mai"=c( 1.5, 0.95, 0.85, 0.4))

		plot( 1,1, type="n", main=paste( "Cell Type Profile:   Gene Unit Vectors\n",label),
			xlim=c(0.5,N_STAGES+0.5), ylim=c(0,yMax), xaxt="n", xlab=NA,
			ylab="Gene 'Projection' per Cell Type", las=3, font.axis=2, font.lab=2)
		axis( 1, at=1:N_STAGES, label=STAGE_NAMES, las=3, font=2, cex.axis=0.8)
	}

	# draw the lines a bit nicer...as step steps...
	colUse <- rep( col, length.out=length(gSet))
	lwdUse <- rep( lwd, length.out=length(gSet))

	gSet <- shortGeneName( gSet, keep=1)
	where <- base::match( gSet, vectorSpace$GENE_ID, nomatch=NA)
	if ( any( is.na(where))) cat( "\nSome Genes not found: ", gSet[ is.na(where)])
	who <- which( !is.na(where))

	for ( k in who) {
		y <- as.numeric( vectorSpace[ where[k], 3:ncol(vectorSpace)])
		ans <- drawCellTypeProfileDensityLine( y, col=colUse[k], lwd=lwdUse[k])
		xLocation <- ans$x
		if ( annotate) {
			useY <- which.max( y)
			text( xLocation[useY], y[useY], gSet[k], col=colUse[k], pos=3, font=2, cex=1)
		}
	}
	if ( ! is.na( legend) && length(who)) {
		legend( x=legend, legend=gSet[who], col=colUse[who], lwd=3, cex=legend.cex)
	}
	dev.flush()
	
	# send back useful info
	units <- vectorSpace[ where, ]
	out <- list( "unitVectors"=units, "x"=xLocation)
	return( out)
}


`plotCellTypeProfileIntensity` <- function( gSet, col=1, lwd=1, pt.cex=1, legend=NA, plot=TRUE, 
				legend.cex=1, label="", minYmin=1, minYmax=NULL, threshold=NULL, annotate=FALSE) {

	CellTypeSetup()

	# get the life cycle data...  there is 'geneID, product', and then the intensities...
	intenSpace <- CellTypeEnv[[ "IntensitySpace"]]
	intenVectors <- as.matrix( intenSpace[ , 3:ncol(intenSpace)])
	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	gSet <- shortGeneName( gSet, keep=1)
	where <- base::match( gSet, intenSpace$GENE_ID, nomatch=0)
	if (sum(where > 0) < 1) {
		cat( "\nNo Matching Genes found..")
		return()
	}

	intenSpace <- intenSpace[ where, , drop=F]
	intenVectors <- intenVectors[ where, , drop=F]
	gSet <- gSet[ where > 0]
	NG <- length(gSet)

	bigValue <- max( intenVectors, na.rm=T)
	if ( ! is.null( minYmax)) bigValue <-max( bigValue, minYmax, na.rm=T)
	if ( ! is.null( threshold)) bigValue <- max( bigValue, threshold, na.rm=T)
	yMax <- bigValue * 1.2
	smallValue <- min( intenVectors, na.rm=T)
	if ( ! is.null( threshold)) minYmin <- min( minYmin, threshold)
	yMin <- min( smallValue, minYmin)
	if ( yMin < 0.01) yMin <- 0.01

	if (plot) {

		par( "mai"=c( 1.5, 0.95, 0.85, 0.4))
		mainText <- paste( "Cell Type Profile:    Gene Intensity Vectors\n",label)
		if ( NG == 1) {
			mainText <- paste( "Cell Type Profile:    Gene Intensity:    ", gSet)
			mainText <- paste( mainText, "\n", gene2Product( gSet))
		}

		plot( 1,1, type="n", main=mainText, xlim=c(0.5,N_STAGES+0.5), ylim=c(yMin,yMax), 
			log="", xaxt="n", xlab=NA, ylab="Gene Intensity per Stage  (RPKM)", las=3, font.axis=2, font.lab=2)
		axis( 1, at=1:N_STAGES, label=STAGE_NAMES, las=3, font=2)
	}

	# draw the lines a bit nicer...as step steps...
	colUse <- rep( col, length.out=length(gSet))
	lwdUse <- rep( lwd, length.out=length(gSet))

	if ( ! is.null( threshold)) {
		lines( c(-2,20), rep.int(threshold,2), lwd=1, lty=3, col="brown")
		text( c(0.3,0.3), c( threshold,threshold), c( "Expression", "Threshold"), col="brown", pos=c(3,1))
	}

	# when only one gene, color by stage instead
	if ( NG == 1) {
		colUse <- CellTypeEnv[[ "STAGE_COLORS"]]
		y <- intenVectors[ 1, ]
		y <- pmax( y, yMin)
		ans <- drawCellTypeProfileIntensityLine( y, col=colUse, lwd=lwdUse[1], pt.cex=pt.cex, col2=1, useLog=F)
		xLocation <- ans$x
		if ( annotate) {
			useY <- which.max( y)
			text( xLocation[useY], y[useY], gSet[1], col=colUse[useY], pos=3, font=2, cex=1)
		}
	} else {
		for ( k in 1:length(gSet)) {
			y <- intenVectors[ k, ]
			y <- pmax( y, yMin)
			ans <- drawCellTypeProfileIntensityLine( y, col=colUse[k], lwd=lwdUse[k], pt.cex=pt.cex, useLog=F)
			xLocation <- ans$x
			if ( annotate) {
				useY <- which.max( y)
				text( xLocation[useY], y[useY], gSet[k], col=colUse[k], pos=3, font=2, cex=1)
			}
		}
	}
	if ( ! is.na( legend)) {
		if ( NG == 1) colUse <- 1
		legend( x=legend, legend=gSet, col=colUse, lwd=3, cex=legend.cex)
	}
	dev.flush()
	
	# send back useful info
	out <- list( "intensityVectors"=intenVectors, "x"=xLocation)
	return( out)
}


`drawCellTypeProfileDensityLine` <- function(y, col=1, lwd=1) {

	yUse <- c( 0, rep( y, each=2), 0)
	xmids <- seq( 0, length(y), by=1) + 0.5
	xUse <- base::sort( c( xmids-0.05, xmids+0.05))
	lines( xUse, yUse, type="l", col=col, lwd=lwd) 
	return( invisible( list( "x"=(xmids+0.5))))
}


`drawCellTypeProfileIntensityLine` <- function(y, at=1:length(y), col=1, lwd=1, pt.cex=1, col2=col, useLog=FALSE) {

	N <- length(y)
	if (useLog) {
		splineAns <- spline( at, log2(y+1), n=N*5)
		xSpline <- splineAns$x
		ySpline <- 2 ^ (splineAns$y) - 1
		# spline in log space can look extreme...
		minY <- min(y, na.rm=T) *.1
		ySpline[ ySpline < minY] <- minY
	} else {
		splineAns <- spline( at, y, n=N*5)
		xSpline <- splineAns$x
		ySpline <- splineAns$y
	}
	lines( xSpline, ySpline, col=col2, lty=1, lwd=lwd)
	points( at, y, pch=21, bg=col, cex=pt.cex)
	return( invisible( list( "x"=at)))
}


# Cell Type Fitting -   try to model a samples cell type profile as a linear combination of
# 						cell type dimensions

`fitCellTypeProfileFromFile` <- function( f, sid="Observed", col="orchid1", 
					geneColumn="GENE_ID", intensityColumn="RPKM_M", 
					sep="\t", max.iterations=100, rate=1, tolerance=0.01,
					makePlots=c("all","final","none"), plot.path=".",
					algorithm=c("steepDescent", "nls", "GenSA"), 
					geneUniverse=NULL, verbose=TRUE, ...) {

	CellTypeSetup()

	# open that file and find the needed columns
	tmp <- read.delim( f, as.is=T, sep=sep)
	if ( nrow( tmp) < 1) return(NULL)

	gset <- tmp[[ geneColumn]]
	if ( is.null(gset)) {
		cat( "calcCellTypeProfile:  gene name column not found.\nFile: ",f,
			"\nTried: ", geneColumn)
		return( NULL)
	}
	inten <- tmp[[ intensityColumn]]
	if ( is.null(inten)) {
		cat( "calcCellTypeProfile:  gene intensity column not found.\nFile: ",f,
			"\nTried: ", intensityColumn)
		return( NULL)
	}

	makePlots <- match.arg( makePlots)
	algorithm <- match.arg( algorithm)
	ans <-  fitCellTypeProfile( gset, inten, sid=sid, col=col, max.iterations=max.iterations, rate=rate, 
				tolerance=tolerance, makePlots=makePlots, plot.path=plot.path, algorithm=algorithm, 
				geneUniverse=geneUniverse, verbose=verbose, ...)
	if (makePlots != "none") {
		# new plot printing wrapper lets us not append the device type suffix
		plotFile <- paste( sid, getCurrentSpeciesFilePrefix(), getCellTypeReference(), "FitProportions", algorithm, sep=".")
		plotFile <- file.path( plot.path, plotFile)
		printPlot( plotFile)
	}
	
	return( ans)
}


`fitCellTypeProfileFromFileSet` <- function( fnames, fids, fcolors=1:length(fids), geneColumn="GENE_ID", 
						intensityColumn="RPKM_M", sep="\t",
						max.iterations=100, rate=1, tolerance=0.01, 
						makePlots=c("all","final","none"), plot.path=".",
						algorithm=c("steepDescent", "nls", "GenSA"), 
						geneUniverse=NULL, verbose=TRUE, ...) {
								
	CellTypeSetup()

	# make sure we can read those files
	filesOK <- file.exists( fnames)
	if ( !all( filesOK)) {
		cat( "\nCellTypeProfile:  Some transcript files not found:\n")
		print( fnames[ !filesOK])
		return(NULL)
	}

	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build the storage to hold the results
	nFiles <- length( fnames)
	m <- matrix( nrow=nFiles, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids
	rmsd <- vector( length=nFiles)
	names(rmsd) <- fids

	# load each file in turn
	makePlots <- match.arg( makePlots)
	algorithm <- match.arg( algorithm)
	cat( "\nFitting N_files: ", nFiles)
	for( i in 1:nFiles) {
		cat( "\n", basename(fnames[i]))
		ans <- fitCellTypeProfileFromFile( fnames[i], sid=fids[i], col=fcolors[i], geneColumn=geneColumn,
					intensityColumn=intensityColumn, sep=sep,
					max.iterations=max.iterations, rate=rate, 
					tolerance=tolerance, makePlots=makePlots, 
					plot.path=plot.path, algorithm=algorithm, 
					geneUniverse=geneUniverse, verbose=verbose, ...)
		m[ i, ] <- ans$CellProportions
		rmsd[i] <- ans$RMSD
		if (makePlots != "none") {
			# new plot printing wrapper lets us not append the device type suffix
			plotFile <- paste( fids[i], getCurrentSpeciesFilePrefix(), getCellTypeReference(), "FitProportions", algorithm, sep=".")
			plotFile <- file.path( plot.path, plotFile)
			printPlot( plotFile)
		}
	}
	cat( "\nDone.\n")

	out <- list( "CellProportions"=t(m), "RMSD"=rmsd)
	return( out)
}


`fitCellTypeProfileFromMatrix` <- function( m, fcolors=1:ncol(m),  max.iterations=100, rate=1, tolerance=0.01, 
						makePlots=c("all","final","none"), plot.path=".", 
						algorithm=c("steepDescent", "nls", "GenSA"), 
						geneUniverse=NULL, verbose=TRUE, ...) {

	CellTypeSetup()

	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build the storage to hold the results
	nSamples <- ncol(m)
	fids <- colnames(m)
	genes <- shortGeneName( rownames(m), keep=1)
	mOut <- matrix( nrow=nSamples, ncol=N_STAGES)
	colnames(mOut) <- STAGE_NAMES
	rownames(mOut) <- fids
	rmsd <- vector( length=nSamples)
	names(rmsd) <- fids

	# load each in turn
	makePlots <- match.arg( makePlots)
	algorithm <- match.arg( algorithm)
	cat( "\nFitting N_Samples: ", nSamples)
	for( i in 1:nSamples) {
		cat( "\n", fids[i])
		ans <- fitCellTypeProfile( genes=genes, inten=m[,i], sid=fids[i], col=fcolors[i], 
					max.iterations=max.iterations, rate=rate, 
					tolerance=tolerance, makePlots=makePlots, 
					plot.path=plot.path, algorithm=algorithm, 
					geneUniverse=geneUniverse, verbose=verbose, ...)
		mOut[ i, ] <- ans$CellProportions
		rmsd[i] <- ans$RMSD
		if (makePlots != "none") {
			# new plot printing wrapper lets us not append the device type suffix
			plotFile <- paste( fids[i], getCurrentSpeciesFilePrefix(), getCellTypeReference(), "FitProportions", algorithm, sep=".")
			plotFile <- file.path( plot.path, plotFile)
			printPlot( plotFile)
		}
	}
	cat( "\nDone.\n")

	out <- list( "CellProportions"=t(mOut), "RMSD"=rmsd)
	return( out)
}


# the function to do the fit to one set of gene intensity
`fitCellTypeProfile` <- function( genes, inten, sid="Observed", col="orchid1", modelCol='brown',
					max.iterations=100, rate=1, tolerance=0.01, fit.starts=NULL,
					makePlots=c("all","final","none"), plot.path=".", sleep=0.01, 
					algorithm=c("steepDescent", "nls", "GenSA"), 
					geneUniverse=NULL, verbose=TRUE, ...) {

	# grab the Cell Type data we will need:  the gene intensity in all cell types
	CellTypeSetup()
	intensitySpace <- CellTypeEnv[[ "IntensitySpace"]]
	intensityVectors <- intensitySpace[ , 3:ncol(intensitySpace)]
	intensityMatrix <- as.matrix( intensityVectors)
	vectorSpace <- CellTypeEnv[[ "VectorSpace"]]
	vectorVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	vectorMatrix <- as.matrix( vectorVectors)
	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# start with the calculated profile for this transcriptome
	genes <- shortGeneName( genes, keep=1)

	# allow being passed in a gene universe of a subset of genes to use in the fit
	if ( ! is.null( geneUniverse)) {
		keep <- which( genes %in% as.character( geneUniverse))
		if ( length( keep) < length(geneUniverse)/2) cat( "\nWarnings:  Trimming to given Gene Universe removed too many genes..")
		if ( length( keep) < N_STAGES*2) {
			stop( "\nError:  Too few genes to run cell type fit modelling.")
		}
		genes <- genes[ keep]
		inten <- inten[ keep]
	}
	NG <- length( genes)

	# make the stage profile of the given data
	obsAns <- calcCellTypeProfile( genes, inten)
	obsProfile <- obsAns$Profile
	whereGene <- match( genes, intensitySpace$GENE_ID, nomatch=0)
	useGene <- which( whereGene > 0)
	makePlots <- match.arg( makePlots)
	
	# build up an initial model of how much intensity goes to each cell type
	if ( is.null( fit.starts)) {
		fit.starts <- rep.int( 1/N_STAGES, N_STAGES)
	} else {
		if ( length(fit.starts) != N_STAGES) stop( "Error:  'fit.starts' must be length of all cell types.")
	}

	# build a small 2x matrix to hold observed and model
	cellM <- matrix( 0, nrow=N_STAGES, ncol=2)
	colnames(cellM) <- c( sid, "ModelFit")
	rownames(cellM) <- STAGE_NAMES
	cellM[ ,1] <- obsProfile

	# make some other storage for the various fit low level functions to use as 'global' variables
	# build the matrix of stage fractions for every gene we were given.
	geneCellFractions <- matrix( 0, nrow=NG, ncol=N_STAGES)
	for( i in 1:N_STAGES) {
		geneCellFractions[ useGene, i] <- vectorMatrix[ whereGene, i]
	}
	nIter <- 0
	lastDrawnPcts <- rep.int( 0, N_STAGES)
	

	# allow more than one fit method
	steepDescentFitCellTypeModel <- function() {
		# Use the starting percentages and build a small 2x matrix to hold observed and model
		model.pcts <- fit.starts
		
		# ready to iteratively compare the model to the observed. 
		# track the best and some metrics for seeing stalling
		if (verbose) cat( "\n")
		prevRMSD <- 10000
		best.model <- model.pcts
		best.rmsd <- prevRMSD
		nTimesStuck <- 0
	
		for ( i in 1:max.iterations) {
			
			# make the latest model with the current ratios, making sure the percentages sum to exactly one
			model.pcts <- model.pcts / sum( model.pcts)
			modelInten <- rep.int( 0, NG)
			for ( j in 1:N_STAGES) {
				thisV <- intensityVectors[,j] * model.pcts[j]
				modelInten[ useGene] <- modelInten[ useGene] + thisV[whereGene]
			}
			
			# calculate the profiles of this model
			modelAns <- calcCellTypeProfile( genes, modelInten)
			modelProfile <- modelAns$Profile
			cellM[ ,2] <<- modelProfile
	
			# assess the current deviation
			deltas <- cellM[,1] - cellM[,2]
			rmsd <- round( sqrt( mean( deltas^2)), digits=4)
			if (verbose) cat( "\rIter: ", i, "   RMSD: ", rmsd)
			
			# plot it to show progress
			if (makePlots == "all") {
				labelText <- paste( "'SteepDescent' Model Fit to:  ", sid, "\nIteration: ", i, "    RMS_Deviation: ", round(rmsd,dig=4))
				mp <- plotCellTypeProfiles( t(cellM), col=c( col, modelCol), label=labelText, ...)
	
				# try to show what we have now?...
				cellPercents <- round( model.pcts * 100, digits=1)
				xShow <- apply( mp, 2, mean)
				yShow <- apply( cellM, 1, max) + (max(cellM)*0.03)
				toShow <- which( cellPercents > 0)
				if (length(toShow)) text( xShow[toShow], yShow[toShow], paste(cellPercents[toShow],"%",sep=""), cex=0.7, col=1)
				dev.flush()
				if (sleep > 0) Sys.sleep( sleep)
			}
			
			if ( rmsd <= tolerance) {
				cat( "\nConverged..")
				break
			}
			deltaRMSD <- prevRMSD - rmsd
			if ( abs(deltaRMSD) < 1e-5) {
				nTimesStuck <- nTimesStuck + 1
			} else {
				nTimesStuck <- 0
			}
			if ( nTimesStuck > 5) {
				cat( "\nStalled..")
				break
			}
			# we did not bail out, see if we are better
			if (rmsd < best.rmsd) {
				best.model <- model.pcts
				best.rmsd <- rmsd
				nTimesStuck <- 0
			}
			prevRMSD <- rmsd
			
			# now use the deltas to push the model percentages
			# turn the 0-100 profiles back to 0-1 fractions, with the rate term
			# try 1:  using a linear correction
			dPcts <- deltas * rate / 100
			model.pcts <- model.pcts + dPcts
			# since the percents all need to sum to 1, try adding the tiny negative to all other cell types
			# so the net model stays about centered
			otherDpcts <- rep.int( 0, N_STAGES)
			for ( k in 1:N_STAGES) {
				myTinyPcts <- rep.int( dPcts[k]/N_STAGES, N_STAGES)
				myTinyPcts[k] <- 0   # but nothing for myself
				otherDpcts <- otherDpcts + myTinyPcts
			}
			model.pcts <- model.pcts - otherDpcts
			# prevent negative contributions
			# try 1 is just to zero them out...
			model.pcts[ model.pcts < 0] <- 0
			
			# and go around again
		}
	
		# fell out of the loop
		if ( i == max.iterations) cat( "\nMax.iterations..")
	
		# regardless of how we ended, use the best we saw
		model.pcts <- best.model
		rmsd <- best.rmsd
		# when done force it to exactly 100%
		model.pcts <- model.pcts / sum(model.pcts,na.rm=T)

		if (makePlots == "final") {
			labelText <- paste( "'SteepDescent' Model Fit to:  ", sid, "\nIteration: ", i, "    RMS_Deviation: ", round(rmsd,dig=4))
			mp <- plotCellTypeProfiles( t(cellM), col=c( col, modelCol), label=labelText, ...)
			# try to show what we have now?...
			cellPercents <- round( model.pcts * 100, digits=1)
			xShow <- apply( mp, 2, mean)
			yShow <- apply( cellM, 1, max) + (max(cellM)*0.03)
			toShow <- which( cellPercents > 0)
			if (length(toShow)) text( xShow[toShow], yShow[toShow], paste(cellPercents[toShow],"%",sep=""), cex=0.7, col=1)
			dev.flush()
		}
		out <- list( "model.pcts"=model.pcts, "rmsd"=rmsd, "iterations"=i)
		return(out)
	}


	nlsFitCellTypeModel <- function() {

		# call the Nonlinear Least Squares (NLS)
		controlList <- nls.control( maxiter=max.iterations, minFactor=1/512, warnOnly=TRUE, tol=1e-3)
		starts <- list( "modelPcts"=fit.starts)
		# allow a small tolerance, so true zero is a valid answer, and force zero limits later
		lowerBound <- rep.int( -0.0000001, N_STAGES)
		# allow any one dimension to exceed 100% and the scale it back later
		upperBound <- rep.int(  2.0, N_STAGES)
		nIter <<- 0
		fitAns <- try( nls( obsProfile ~ nlsModelCellTypeProfile( modelPcts), start=starts,
				control=controlList, algorithm="port", lower=lowerBound, 
				upper=upperBound))

		if ( class( fitAns) == "try-error") {
			cat( "\nFitting of Cell Type Proportions failed...")
			return(NULL)
		} 

		model.pcts <- coef( fitAns)
		model.pcts[ model.pcts < 0] <- 0
		names(model.pcts) <- STAGE_NAMES
		resids <- residuals(fitAns)
		rmsd <- sqrt( mean( resids^2))

		# when done force it to exactly 100%
		model.pcts <- model.pcts / sum(model.pcts,na.rm=T)

		out <- list( "model.pcts"=model.pcts, "rmsd"=rmsd, "Iterations"=nIter)

		if (makePlots == "final") {
			labelText <- paste( "'NLS' Model Fit to:  ", sid, "\nIteration: ", nIter, "    RMS_Deviation: ", round(rmsd,dig=4))
			mp <- plotCellTypeProfiles( t(cellM), col=c( col, modelCol), label=labelText, ...)
			# try to show what we have now?...
			cellPercents <- round( model.pcts * 100, digits=1)
			xShow <- apply( mp, 2, mean)
			yShow <- apply( cellM, 1, max) + (max(cellM)*0.03)
			toShow <- which( cellPercents > 0)
			if (length(toShow)) text( xShow[toShow], yShow[toShow], paste(cellPercents[toShow],"%",sep=""), cex=0.7, col=1)
			dev.flush()
		}

		return( out)
	}


	nlsModelCellTypeProfile <- function( modelPcts) {

		# low level function called by NLS to give the current cell type profile of a set of pcts
		model.pcts <- modelPcts
		model.pcts[ model.pcts < 0] <- 0
		modelInten <- rep.int( 0, NG)
		for ( j in 1:N_STAGES) {
			thisV <- intensityMatrix[,j] * model.pcts[j]
			modelInten[ useGene] <- modelInten[ useGene] + thisV[whereGene]
		}
			
		# calculate the profiles of this model using faster code
		modelProfile <- nlsCalcCellTypeProfile( modelInten)

		# force the NLS to try to keep the percentages near 1.0, but don't let that penalty
		# affect what we show in the plots
		totalPct <- sum(modelPcts)
		modelProfileShow <- modelProfile
		if ( totalPct > 1.0) {
			modelProfileShow <- modelProfile / totalPct
			modelProfile <- modelProfile * totalPct * totalPct
		}

		cellM[ ,2] <<- modelProfileShow

		# plot it to show progress
		if (makePlots == "all" && any( abs(model.pcts - lastDrawnPcts) > 0.001)) {
			nIter <<- nIter + 1
			lastDrawnPcts <<- model.pcts
			# remember the internal values are 0..1 while the external are 0..100
			deltas <- cellM[,1] - cellM[,2]
			rmsd <- round( sqrt( mean( deltas^2)), digits=4)
			labelText <- paste( "'NLS' Model Fit to:  ", sid, "\nIteration: ", nIter, "    RMS_Deviation: ", round(rmsd,dig=4))
			mp <- plotCellTypeProfiles( t(cellM), col=c( col, modelCol), label=labelText, ...)

			# try to show what we have now?...
			cellPercents <- round( model.pcts * 100, digits=1)
			xShow <- apply( mp, 2, mean)
			yShow <- apply( cellM, 1, max) + (max(cellM)*0.03)
			toShow <- which( cellPercents > 0)
			if (length(toShow)) text( xShow[toShow], yShow[toShow], paste(cellPercents[toShow],"%",sep=""), cex=0.7, col=1)
			dev.flush()
			if (sleep > 0) Sys.sleep( sleep)
		}
			
		return( modelProfile)
	}


	# faster local version for NLS repeated calls
	nlsCalcCellTypeProfile <- function( inten) {

		# spread these intensitys over those stages
		myVectors <- matrix( 0, nrow=NG, ncol=N_STAGES)
		colnames(myVectors) <- STAGE_NAMES
		for( i in 1:N_STAGES) {
			myVectors[ , i] <- inten * geneCellFractions[ , i]
		}

		# now do a summation by stage, and express as percentages...
		allSums <- apply( myVectors, MARGIN=2, FUN=sum, na.rm=T)
		bigSum <- sum( allSums)
		ans <- allSums * 100 / bigSum
		return( ans)
	}
	

	GenSA.FitCellTypeModel <- function() {

		# call the Gen Simulated Annealing
		starts <- fit.starts
		lowerBound <- rep.int( -0.0000001, N_STAGES)   # allow a small tolerance, so true zero is a valid answer
		upperBound <- rep.int(  1.0000001, N_STAGES)
		nIter <<- 0
		fitAns <- try( do.FitCellTypeModel.GenSA( obsProfile, start=starts,
				lower=lowerBound, upper=upperBound))
		if ( class( fitAns) == "try-error") {
			cat( "\nFitting of Cell Type Proportions failed...")
			return(NULL)
		} 

		# when done force it to exactly 100%
		model.pcts <-  fitAns$BestFit / sum( fitAns$BestFit)		
		rmsd <- fitAns$RMS.Deviation
		out <- list( "model.pcts"=model.pcts, "rmsd"=rmsd, "Iterations"=nIter)

		if (makePlots == "final") {
			labelText <- paste( "'GenSA' Model Fit to:  ", sid, "\nIteration: ", nIter, "    RMS_Deviation: ", round(rmsd,dig=4))
			mp <- plotCellTypeProfiles( t(cellM), col=c( col, modelCol), label=labelText, ...)
			# try to show what we have now?...
			cellPercents <- round( model.pcts * 100, digits=1)
			xShow <- apply( mp, 2, mean)
			yShow <- apply( cellM, 1, max) + (max(cellM)*0.03)
			toShow <- which( cellPercents > 0)
			if (length(toShow)) text( xShow[toShow], yShow[toShow], paste(cellPercents[toShow],"%",sep=""), cex=0.7, col=1)
			dev.flush()
		}
		return( out)
	}


	`do.FitCellTypeModel.GenSA` <- function( obsProfile, start, lower, upper, seed=NULL) {

		# wrapper to implement cell type modelling by Generalize Simulated Annealing

		# GenSA as implemented has a hard coded seed for RNG.  We want random behavior each time thru
		# GEnSA docs suggest using a negative value
		my.seed <- -(as.integer( Sys.time()))

		# the penalty function that GenSA will minimize
		genSA.profile.model <- function( wts) {

			# make the model of gene intensity these weights would generate
			modelInten <- rep.int( 0, NG)
			for ( j in 1:N_STAGES) {
				thisV <- intensityVectors[,j] * wts[j]
				modelInten[ useGene] <- modelInten[ useGene] + thisV[whereGene]
			}
			# calculate the cell type profiles of this model
			modelProfile <- nlsCalcCellTypeProfile( modelInten)
			return( modelProfile)
		}

		# the penalty function that GenSA will minimize
		genSA.profile.residual <- function( wts, obsProfile) {

			modelProfile <- genSA.profile.model( wts)
			
			# use residual sum of squares as the penalty
			resid2 <- (obsProfile - modelProfile) ^ 2
			rss <- sum( resid2, na.rm=T)
			
			# we want the weights to sum to 1.0
			ttlWt <- sum( wts)
			if ( ttlWt > 1.0) rss <- rss * ttlWt * ttlWt			
			return( rss)
		}

		# set up to call GenSA
		# say 'good enough' if we explain 99.99% of the profile
		stopValue <- 0.001
		control.list <- list( "maxit"=5000, "threshold.stop"=stopValue, 
				"temperature"=6000, "smooth"=FALSE, "max.call"=10000000,
				"max.time"=100, "trace.mat"=TRUE, "seed"=my.seed)
				
		# starting profile is equal wt for all cell types
		wts <- start

		# the GenSA tool seems unable to push any parameter all the way to zero, as if it can't 
		# equal the lower bound, instead must stay above it
		# so let's let the lower bounds go a bit bit negative, and clean up later
		lower[ lower == 0] <- -0.01
		wts[ wts <= lower] <- 0.001

		ans <- GenSA( par=wts, lower=lower, upper=upper, fn=genSA.profile.residual, 
			control=control.list, obsProfile=obsProfile)

		# extract the answers
		fractions <- ans$par
		# Clean up:  don't let any final calls be below zero percent, and scale back to exactly 1.0
		fractions[ fractions < 0] <- 0
		fractions <- fractions / sum(fractions)
		names( fractions) <- STAGE_NAMES
		modelProfile <- genSA.profile.model( fractions)
		cellM[ ,2] <<- modelProfile
		resids <- obsProfile - modelProfile
		rms <- sqrt( mean( resids^2))
		# let's do a few other stats...
		cor.ans <- cor.test( obsProfile, modelProfile)
		r2.pearson <- cor.ans$estimate ^ 2
		pv <- cor.ans$p.value
		# also do the more general coefficient of determination
		meanI <- mean( obsProfile, na.rm=T)
		SStotal <- sum( (obsProfile - meanI) ^ 2)
		SSresid <- sum( resids ^ 2)
		r2.cod <- 1.0 - (SSresid / SStotal)

		out <- list( "BestFit"=fractions, "Observed"=obsProfile, "Model"=modelProfile, "Residuals"=resids, 
			"RMS.Deviation"=rms, "R2.CoD"=r2.cod, "R2.Pearson"=r2.pearson, "Pvalue"=pv)
		return( out)
	}
	# end of local fit functions
	

	# call the one fit function
	algorithm <- match.arg( algorithm)
	if ( algorithm == "steepDescent") {
		ans <- steepDescentFitCellTypeModel()
	} else if ( algorithm == "nls") {
		ans <- nlsFitCellTypeModel()
	} else {
		require( GenSA)
		ans <- GenSA.FitCellTypeModel()
	}
	model.pcts <- ans$model.pcts
	rmsd <- ans$rmsd
	iterations <- ans$iterations

	# once we drop out, assemble the results
	cellPercents <- round( model.pcts * 100, digits=4)
	names(cellPercents) <- STAGE_NAMES

	# rather than send back percentages that must sum to 100, let's instead make this 
	# be how many units of each cell type, such that the total intensity would best match
	# the total intensity of the given input.
	observedSum <- sum( inten, na.rm=T)
	modelSums <- apply( intensityMatrix[ whereGene, ], 2, sum, na.rm=T)
	scaleFac <- observedSum / mean( modelSums)
	cellProportions <- cellPercents * scaleFac

	out <- list( "Iterations"=iterations, "RMSD"=rmsd, "CellProportions"=cellProportions)
	return( out)
}


`plotCellTypeVolcano` <- function( file, geneColumn="GENE_ID", foldColumn="LOG2FOLD", pvalueColumn="AVG_PVALUE", 
			keepIntergenics=FALSE, label="Plot", cex=1, sep="\t", cut.fold=1, cut.pvalue=0.05, pch=21,
			marker.genes=NULL, marker.col=1, marker.cex=1, marker.labels=TRUE, marker.pch=21, 
			marker.pos=NULL, min.intensity=0, intensityColumn="RPKM_1", 
			left.label=NULL, right.label=NULL, forceYmax=NULL, ...) {

	if ( is.character(file)) {
		tmp <- read.delim( file, as.is=T, sep=sep)
		cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tmp))
	} else if (is.data.frame(file)) {
		tmp <- file
	} else {
		stop( "Argument 'file' must be a character string or a data frame.")
	}

	if ( !( all( c( geneColumn, foldColumn, pvalueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, foldColumn, pvalueColumn, 
				"\n  \tFound: ", colnames(tmp))
		return()
	}

	# extract the parts we want
	genes <- tmp[[ geneColumn]]
	genes <- shortGeneName( genes, keep=1)
	fold <- as.numeric( tmp[[ foldColumn]])
	pval <- as.numeric( tmp[[ pvalueColumn]])
	celltype <- if ( "CellType" %in% colnames(tmp)) tmp$CellType else gene2CellType(genes)

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		# gmap <- getCurrentGeneMap()
		# dropableGenes <- subset( gmap, SEQ_ID %in% c( "Pf3D7_PFC10_API_IRAB", "Pf3D7_M76611"))$GENE_ID
		# drops2 <- which( genes %in% dropableGenes)
		# drops <- base::sort( base::union( drops, drops2))
	}
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		fold <- fold[ -drops]
		pval <- pval[ -drops]
		celltype <- celltype[ -drops]
		cat( "\nAfter dropping non-genes: ", length(genes))
	}

	# allow the removal of very low expression genes
	if ( min.intensity > 0) {
		inten <- tmp[[ intensityColumn]]
		drops <- which( inten < min.intensity)
		if ( length(drops)) {
			genes <- genes[ -drops]
			fold <- fold[ -drops]
			pval <- pval[ -drops]
			celltype <- celltype[ -drops]
			cat( "\nAfter dropping low intensity: ", length(genes))
		}
	}

	# get the cell type colors to use/show.  The color map is short names, but the data passed in
	# may have the longer human readable names
	cellColors <- getCellTypeColors()
	cellNames <- names(cellColors)
	geneCellColor <- cellColors[ match( celltype, cellNames)]

	# do the volcano in line now
	# we are plotting -log10(pval) on Y axis, Fold on X axis
	clip.fold <- 10
	clip.pvalue <- 1e-10
	# prevent absurd numbers from dominating
	fold[ fold > clip.fold] <- clip.fold
	fold[ fold < -clip.fold] <- -clip.fold
	pval[ pval < clip.pvalue] <- clip.pvalue

	# add extra room on X for the labels, and the cell types too
	bigX <- max( 1, abs( fold), na.rm=T)
	myRangeX <- range( c( -1, 1, fold), na.rm=T) * 1.15
	myRangeX[1] <- myRangeX[1] - diff(myRangeX)*0.1

	y <- -( log10( pval))
	myRangeY <- c( 0, max( 1, y, na.rm=F)*1.05)
	if ( ! is.null( forceYmax)) myRangeY[2] <- as.numeric(forceYmax)

	# plot so the importent genes get drawn last
	ord <- order(y)
	plot ( fold[ord], y[ord], type="p", main=label, xlab="Log2 Fold Change",
		ylab="-Log10 P", xlim=myRangeX, ylim=myRangeY, 
		pch=pch, col=geneCellColor[ord], bg=geneCellColor[ord], cex=cex, font.axis=2, font.lab=2, ...)

	# label selected genes, assumes the genes 'up in set1' are in the first half...
	# add labels when the P-value is great enough
	if ( ! is.null( marker.genes)) {
		doLabel <- whereMarker
		myPos <- ifelse( fold[whereMarker] > 0, 4, 2)
	} else {
		doLabel <- which( pval < cut.pvalue & abs(fold) > cut.fold)
		myPos <- ifelse( fold[doLabel] > 0, 4, 2)
	}
		
	if ( is.null( marker.genes) && length( doLabel)) {
		if ( length(doLabel) > 2) {
			require( plotrix)
			myPos <- thigmophobe( fold[doLabel], y[doLabel])
			# but force these to take right sides
			myPos[ fold[doLabel] < 0 & myPos == 4] <- 2
			myPos[ fold[doLabel] > 0 & myPos == 2] <- 4
		}
		text( fold[doLabel], y[doLabel], genes[doLabel], pos=myPos, cex=marker.cex)
	}

	if ( length( marker.genes) > 0) {
		if ( marker.genes[1] == "identify") {
			identify( fold, y, shortGeneName(genes, keep=1), col=marker.col, cex=marker.cex)
		} else {
			who <- whereMarker
			marker.genes <- marker.genes[ who > 0]
			if ( length( marker.col) > 1) marker.col <- marker.col[ who > 0]
			who <- who[ who > 0]
			if ( length(who) > 0) {
				# put name left for down reg, right for up
				if ( is.null( marker.pos)) {
					pos <- ifelse( fold[who] < 0, 2, 4)
				} else {
					pos <- marker.pos
				}
				points( fold[who], y[who], col=marker.col, bg=marker.col, pch=marker.pch, cex=marker.cex)
				if (marker.labels) text( fold[who], y[who], genes[who], pos=pos, col=col[3], cex=marker.cex)
			}
		}
	}

	# optional labels to remind which group is which
	if ( !is.null(left.label)) text( myRangeX[1]*.75, myRangeY[2]*0.025, paste( "UP in Group '", left.label, "'", sep=""), cex=1, font=2)
	if ( !is.null(right.label)) text( myRangeX[2]*.75, myRangeY[2]*0.025, paste( "UP in Group '", right.label, "'", sep=""), cex=1, font=2)

	# and show the cell type color legend
	legend( 'topleft', names(cellColors), fill=cellColors, bg='white', cex=0.8)
	dev.flush()
	
	return( invisible( list( "x"=fold, "y"=y, "id"=genes )))
}


# modified version of a volcano plot, that makes one circle per cell type
`plotCellTypeClusters` <- function( file, geneColumn="GENE_ID", foldColumn="LOG2FOLD", pvalueColumn="AVG_PVALUE", 
					gene.pct=5.0, min.enrichment=1.2, max.Pvalue=0.01, label="", sep="\t", label.cex=1, pt.cex=0.65,
					left.label=NULL, right.label=NULL, forceXmax=NULL, forceYmax=NULL, 
					color.alpha=0.75, label.offset.cex=1, legend.cex=0.9, ...) {

	require( plotrix)

	if ( is.character(file)) {
		tmp <- read.delim( file, as.is=T, sep=sep)
		cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tmp))
	} else if (is.data.frame(file)) {
		tmp <- file
	} else {
		stop( "Argument 'file' must be a character string or a data frame.")
	}
	if ( !( all( c( geneColumn, foldColumn, pvalueColumn) %in% colnames(tmp)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, foldColumn, pvalueColumn, 
			"\n  \tFound: ", colnames(tmp))
		return()
	}

	# extract the parts we want
	genes <- tmp[[ geneColumn]]
	genes <- shortGeneName( genes, keep=1)
	fold <- as.numeric( tmp[[ foldColumn]])
	pval <- as.numeric( tmp[[ pvalueColumn]])
	celltype <- if ( "CellType" %in% colnames(tmp)) tmp$CellType else gene2CellType(genes, max.type=5)

	# Note:  With the new 2+ types per gene, with percentages, we have to do something different
	# For now, just keep/use the first (biggest %) one.  May be some way to prorate, but not now...
	topCellType <- sub( "\\:[0-9]+\\%.*", "", celltype)

	# always remove  non genes, etc.
	drops <- grep( "(ng)", genes, fixed=TRUE)
	if ( length( drops) > 0) {
		genes <- genes[ -drops]
		fold <- fold[ -drops]
		pval <- pval[ -drops]
		topCellType <- topCellType[ -drops]
	}
	NG <- length(genes)

	# make sure we are in DE order
	ord <- diffExpressRankOrder( fold, pval)
	genes <- genes[ord]
	fold <- fold[ord]
	pval <- pval[ord]
	topCellType <- topCellType[ord]
	# drawing order is different, to overlay the more DE on top
	ord2 <- rev( diffExpressRankOrder( abs(fold), pval))

	# get the cell type colors to use/show.  The color map is short names, but the data passed in
	# may have the longer human readable names
	cellColors <- getCellTypeColors()
	cellNames <- names(cellColors)
	# make some transparent colors too
	rgbCol <- col2rgb( cellColors)
	cellTransparentColors <- rgb( t(rgbCol)/256, alpha=color.alpha)
	geneCellColor <- cellColors[ debugPtrs <- match( topCellType, cellNames)]

	# we are plotting -log10(pval) on Y axis, Fold on X axis
	clip.fold <- 10
	clip.pvalue <- 1e-10
	# prevent absurd numbers from dominating
	fold[ fold > clip.fold] <- clip.fold
	fold[ fold < -clip.fold] <- -clip.fold
	pval[ pval < clip.pvalue] <- clip.pvalue
	y <- -( log10( pval))
	crop.x <- clip.fold
	crop.y <- -( log10( clip.pvalue))
	NG.use <- round( NG * (gene.pct/100) / 10) * 10

	# we wiil look at the top UP and DOWN genes
	# watch the fold change to make sure we don't cross zero
	UPgenes <- 1:NG.use
	DOWNgenes <- (NG-NG.use+1):NG
	UPgenes <- intersect( UPgenes, which(fold > 0))
	DOWNgenes <- intersect( DOWNgenes, which(fold < 0))

	# use cell type enrichment to decide who to highlight
	# as of new cell types allowing 2+ types, all cell type names are always 'short'
	# since we will decide later which balloons to draw, get all the enrichment details, not just the significant here...
	UPenrich <- cellTypeEnrichment( topCellType[UPgenes], mode="genes", minEnrich=1, maxPvalue=1, upOnly=F, verbose=F)
	DOWNenrich <- cellTypeEnrichment( topCellType[DOWNgenes], mode="genes", minEnrich=1, maxPvalue=1, upOnly=F, verbose=F)
	if ( is.null(UPenrich) || is.null(DOWNenrich)) return(NULL)

	# now let's calculate the clusters for each cell type, both UP and DOWN
	cellFac <- factor( topCellType)
	outCell <- outDir <- outCount <- outGenes <- outPct <- outFold <- outPval <- vector()
	outRadius <- outColor <- vector()
	nout <- 0

	# decide a scaling factor for turning percentages into a radius
	# lets say 100% should fill the area from Fold=0 to max X
	# will be applied to the sqrt of the Pct of Genes, so sqrt(100%) = 10%
	myBigX <- quantile( abs(fold), 0.999, na.rm=F)
	radius.scale.fac <- (myBigX * 0.5) * 0.10

	tapply( 1:NG, cellFac, function(k) {
		
		# given all the genes for one cell type, bail if not a real cell type
		ct <- topCellType[k[1]]
		if ( is.null(ct) || is.na(ct) || length(ct) < 1 || ct == "") return()
		ctColor <- cellTransparentColors[ match( ct, cellNames)]

		# see how many and where each group falls
		xUP <- yUP <- xDOWN <- yDOWN <- 0
		gUP <- gDOWN <- ""
		kUP <- intersect( k, UPgenes)
		nUP <- length(kUP)
		pctUP <- round( nUP * 100 / NG.use, digits=1)
		radUP <- sqrt(pctUP) * radius.scale.fac
		if (nUP) {
			xUP <- mean( fold[kUP], na.rm=T)
			yUP <- mean( y[kUP], na.rm=T)
			gUP <- paste( genes[kUP], collapse="; ")
		} 

		kDOWN <- intersect( k, DOWNgenes)
		nDOWN <- length(kDOWN)
		pctDOWN <- round( nDOWN * 100 / NG.use, digits=1)
		radDOWN <- sqrt(pctDOWN) * radius.scale.fac
		if (nDOWN) {
			xDOWN <- mean( fold[kDOWN], na.rm=T)
			yDOWN <- mean( y[kDOWN], na.rm=T)
			gDOWN <- paste( genes[kDOWN], collapse="; ")
		} 

		# save all the details
		now <- (nout+1) : (nout+2)
		outCell[now] <<- ct
		outColor[now] <<- ctColor
		outDir[now] <<- c( "UP", "DOWN")
		outCount[now] <<- c( nUP, nDOWN)
		outGenes[now] <<- c( gUP, gDOWN)
		outPct[now] <<- c( pctUP, pctDOWN)
		outRadius[now] <<- c( radUP, radDOWN)
		outFold[now] <<- c( xUP, xDOWN)
		outPval[now] <<- c( yUP, yDOWN)
		nout <<- nout + 2
	})

	out <- data.frame( "CellType"=outCell, "Direction"=outDir, "N_Genes"=outCount, "Pct_Genes"=outPct,
				"Color"=outColor, "Radius"=outRadius, "Log2Fold"=outFold, "Log10.Pvalue"=outPval,
				"GeneList"=outGenes, stringsAsFactors=F)

	# add extra room on X for the labels, and the cell types too
	# and retune the Y axis limits
	bigX <- max( 1, quantile( abs(fold), 0.999, na.rm=F), abs(out$Log2Fold)+out$Radius)
	myRangeX <- c( -bigX, bigX)
	if ( ! is.null( forceXmax)) {
		bigX <- as.numeric( forceXmax)
		myRangeX[1] <- -bigX
		myRangeX[2] <- bigX
	}
	# force all dots to be seen, whether we crop or not
	crop.x <- min( crop.x, quantile(fold, 0.999, na.rm=T))
	# don't let the crop be smaller than the balloons, or too small overall
	big.balloon <- max( abs(out$Log2Fold) + out$Radius)
	if (crop.x < big.balloon) crop.x <- big.balloon
	if (crop.x < 0.5) crop.x <- 0.5
	fold[ fold > crop.x] <- crop.x
	fold[ fold < -crop.x] <- -crop.x
	myRangeX[1] <- myRangeX[1] - diff(myRangeX)*0.175
	bigY <- max( 1, quantile( y, 0.999, na.rm=T), out$Log10.Pvalue+(out$Radius * 1.15))
	myRangeY <- c( 0, bigY)
	if ( ! is.null( forceYmax)) {
		bigY <- as.numeric(forceYmax)
		myRangeY[2] <- bigY
	}
	crop.y <- min( crop.y, max(y), bigY)
	y[ y > crop.y] <- crop.y

	# plot the volcano as a dust cloud, to downplay the genes
	mainText <- paste( "Volcano by Cell Type:   Enrichment of top ", gene.pct, 
			"% most DE Genes: (N=",NG.use," up, ", NG.use, " down)\n", label, sep="")
	plot( fold[ord2], y[ord2], type="p", main=mainText, xlab="Log2 Fold Change",
		ylab="-Log10 P", xlim=myRangeX, ylim=myRangeY, 
		pch=19, col=geneCellColor[ord2], cex=pt.cex, font.axis=2, font.lab=2, cex.main=0.8, ...)

	# show cropping lines?
	if ( pt.cex > 0.25) {
		lines( c(crop.x,crop.x), c(0,crop.y), col='grey40', lty=3, lwd=1)
		text( crop.x, min(1.6,bigY/2), "Crop Fold Change", col='grey40', srt=90, pos=4, cex=legend.cex)
		lines( c(-crop.x,-crop.x), c(0,crop.y), col='grey40', lty=3, lwd=1)
		text( -crop.x, min(2.6,bigY/2), "Crop Fold Change", col='grey40', srt=90, pos=2, cex=legend.cex)
		lines( c(-crop.x,crop.x), c(crop.y,crop.y), col='grey40', lty=3, lwd=1)
		text( -0, crop.y, "Crop -Log10 P", col='grey40', srt=0, pos=3, cex=legend.cex, offset=0.25)
		# then do the points again to emphasize
		redraw1 <- which( fold[ord2] >= bigX*0.9)
		redraw2 <- which( fold[ord2] <= -bigX*0.9)
		redraw3 <- which( y[ord2] >= bigY*0.9)
		redraw <- union( c( redraw1, redraw2), redraw3)
		if ( length(redraw)) {
			points( fold[ord2][redraw], y[ord2][redraw], pch=19, col=geneCellColor[ord2][redraw], cex=pt.cex)
		}
	}

	# now with all cell clusters known, we can draw their balloons
	ord <- order( out$Radius, decreasing=T)
	out <- out[ ord, ]
	out$Enrichment <- 1
	out$Enrichment.Pvalue <- 1
	toShow <- vector()
	
	# do in two passes, to see who should be drawn, and then repeat to draw them
	for ( i in 1:nrow(out)) {
		# use the enrichment and P-value to decide who to show
		myCellType <- out$CellType[i]
		myUpDown <- out$Direction[i]
		if ( myUpDown == "UP") {
			where <- match( myCellType, UPenrich$CellType)
			myEnrich <- UPenrich$Enrichment[where]
			myPval <- UPenrich$P_Value[where]
		} else {
			where <- match( myCellType, DOWNenrich$CellType)
			myEnrich <- DOWNenrich$Enrichment[where]
			myPval <- DOWNenrich$P_Value[where]
		}
		out$Enrichment[i] <- myEnrich
		out$Enrichment.Pvalue[i] <- myPval
		if ( ! is.na(where) && myEnrich >= min.enrichment && myPval <= max.Pvalue) {
			#draw.circle( out$Log2Fold[i], out$Log10.Pvalue[i], out$Radius[i], border=1, col=out$Color[i])
			toShow <- c( toShow, i)
		}
	}
	out$Drawn <- FALSE
	out$Drawn[toShow] <- TRUE
	
	# the to draw criteria was on enrichment alone.  But there may be cell types more UP and significant
	# that don't quite meet the enrichment criteria.  Find them and call them to be drawn to
	toShowUp <- intersect( toShow, which( out$Log2Fold > 0))
	toShowDown <- intersect( toShow, which( out$Log2Fold < 0))
	if ( length(toShowUp)) {
		minUpFold <- min( out$Log2Fold[toShowUp], na.rm=T)
		minUpPval <- min( out$Log10.Pvalue[toShowUp], na.rm=T)
		minUpRad <- min( out$Radius[toShowUp], na.rm=T)
		extraUp <- which( out$Log2Fold > minUpFold & out$Log10.Pvalue > minUpPval & out$Radius > minUpRad)
		if ( length(extraUp)) toShow <- sort( union( toShow, extraUp))
	}
	if ( length(toShowDown)) {
		minDownFold <- max( out$Log2Fold[toShowDown], na.rm=T)
		minDownPval <- min( out$Log10.Pvalue[toShowDown], na.rm=T)
		minDownRad <- min( out$Radius[toShowDown], na.rm=T)
		extraDown <- which( out$Log2Fold < minDownFold & out$Log10.Pvalue > minDownPval & out$Radius > minDownRad)
		if ( length(extraDown)) toShow <- sort( union( toShow, extraDown))
	}
	
	# draw every one from the best to the last good one worth showing by enrichment
	if ( length(toShow)) {
		for ( i in toShow) {
			draw.circle( out$Log2Fold[i], out$Log10.Pvalue[i], out$Radius[i], border=1, col=out$Color[i])
		}
		toShow <- sort( union( toShow, 1:max(toShow)))
	}
	out$Drawn[toShow] <- TRUE
	
	if ( length( toShow)) {
		P.text <- as.PvalueText( out$Enrichment.Pvalue[toShow], digits=3) 
		ctLabels <- paste( out$CellType[toShow], "\n(N=", out$N_Genes[toShow], ", P=", P.text, ")", sep="")
		if ( length(toShow) > 2) {
			thigmophobe.labels( out$Log2Fold[toShow], out$Log10.Pvalue[toShow], label=ctLabels, col=1, cex=label.cex,
					offset=(out$Radius[toShow]*label.offset.cex))
		} else {
			text( out$Log2Fold[toShow], out$Log10.Pvalue[toShow], label=ctLabels, col=1, cex=label.cex,
					offset=(out$Radius[toShow]*label.offset.cex), pos=c(2,4))
		}
	}

	# optional labels to remind which group is which
	if ( !is.null(left.label)) text( myRangeX[1]*.65, myRangeY[2]*0.025, paste( "UP in '", left.label, "'", sep=""), cex=1, font=2)
	if ( !is.null(right.label)) text( myRangeX[2]*.75, myRangeY[2]*0.025, paste( "UP in '", right.label, "'", sep=""), cex=1, font=2)

	# and show the cell type color legend
	legend( 'topleft', names(cellColors), fill=cellColors, bg='white', cex=legend.cex)
	dev.flush()
	return( invisible( out))
}


`writeCellTypeClusterExtras` <- function( tbl, resultsfile, resultsTbl=NULL, reference=getCellTypeReference()) {

	# given a data frame of details from the cell type volcano cluster plot tool (above), make some supporting files
	path <- dirname( resultsfile)
	file.root <- sub(  ".JOINED.txt$", "", basename( resultsfile))

	extras.path <- file.path( path, paste( reference, "SupplementalFiles", sep="."))
	if ( ! file.exists( extras.path)) dir.create( extras.path, recursive=T)

	# first write out the overall table of results
	outfile <- file.path( extras.path, paste( file.root, reference, "Enrichment.Overview.csv", sep="."))
	write.csv( tbl, outfile, row.names=F)

	# then visit each cell type cluster that did get drawn
	drawnRows <- which( tbl$Drawn == TRUE)
	if ( ! length( drawnRows)) return(0)

	# get the content, so we can extract subsets
	if ( is.null( resultsTbl)) resultsTbl <- read.delim( resultsfile, as.is=T)
	resultGenes <- shortGeneName( resultsTbl$GENE_ID, keep=1)

	for ( i in drawnRows) {
		celltype <- tbl$CellType[i]
		direct <- tbl$Direction[i]
		genes <- strsplit( tbl$GeneList[i], split="; ")[[1]]

		# part 1:  find those genes in the results and make a small table of just those genes in the cluster
		wh <- match( genes, resultGenes, nomatch=0)
		sml <- resultsTbl[ wh, ]
		sml <- data.frame( "CellType"=celltype, "Direction"=direct, sml, stringsAsFactors=F)
		outfile <- file.path( extras.path, paste( file.root, celltype, direct, "MemberGenes.txt", sep="."))
		write.table( sml, outfile, sep="\t", quote=F, row.names=F)

		# part 2: what pathways do these genes enrich for
		gsAns <- geneSetBestMatch( genes, nBest=20)
		outfile <- file.path( extras.path, paste( file.root, celltype, direct, "PathwayHits.csv", sep="."))
		write.csv( gsAns, outfile, row.names=F)

		# others as we add them...
	}
	return( length(drawnRows))
}


# show expression levels of cell sorting marker genes
`pipe.PlotCellSortGatingGenes` <- function( sampleID, annotationFile="Annotation.txt", optionsFile="Options.txt",
					results.path=NULL, label="", geneColumn="GENE_ID", intensityColumn=NULL, 
					sep="\t", col=NULL, sortingGenes=NULL, useLog=F, maxYshow=1000) {

	annT <- readAnnotationTable( annotationFile)
	optT <- readOptionsTable( optionsFile)
	if ( is.null(results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound="results", verbose=F)
	}
	speciesID <- getCurrentSpecies()
	prefix <- getCurrentSpeciesFilePrefix()

	# allow the ID to be a sample name, a expression data frame, a transcriptome file, or a cell type
	tmpDF <- NULL
	targetColumn <- NA
	if ( sampleID %in% annT$SampleID) {
		file <- file.path( results.path, "transcript", paste( sampleID, prefix, "Transcript.txt", sep="."))
		if ( ! file.exists(file)) {
			cat( "\nTranscriptome file not found: ", file)
			return(NULL)
		}
		tmpDF <- read.delim( file, as.is=T)
		if ( ! nrow(tmpDF)) tmpDF <- NULL
		if ( is.null(col) && "Color" %in% colnames(annT)) col <- annT$Color[ match( sampleID, annT$SampleID)]
		if ( label == "") label <- paste( "Sample=", sampleID)
	} else if (is.data.frame(sampleID)) {
		tmpDF <- sampleID
	} else if ( is.character(sampleID) && file.exists(sampleID)) {
		file <- sampleID
		tmpDF <- read.delim( file, as.is=T, sep=sep)
		cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tmpDF))
		if ( ! nrow(tmpDF)) tmpDF <- NULL
		if ( label == "") label <- paste( "File=", basename(sampleID))
	} else {
		targetM <- getCellTypeMatrix()
		targetColors <- getCellTypeColors()
		targetColumn <- match( sampleID, colnames(targetM), nomatch=NA) 
		if ( is.null(col)) col <- targetColors[ match( sampleID, names(targetColors))]
		if ( label == "") label <- paste( "Cell Type=", sampleID)
	}
	if ( is.null(col)) col="tan"

	# extract what we need, given what we were given
	genes <- inten <- NULL
	if ( is.null( intensityColumn)) intensityColumn <- getExpressionUnitsColumn( optionsFile, verbose=F)
	if ( ! is.null( tmpDF)) {
		if ( !( all( c( geneColumn, intensityColumn) %in% colnames(tmpDF)))) {
			cat( "\nMissing columns:  looked for: ", geneColumn, intensityColumn,
				"\n  \tFound: ", colnames(tmpDF))
			return(NULL)
		}
		genes <- shortGeneName( tmpDF[[ geneColumn]], keep=1)
		inten <- as.numeric( tmpDF[[ intensityColumn]])
	}
	if ( ! is.na( targetColumn)) {
		genes <- rownames( targetM)
		inten <- targetM[ , targetColumn]
	}
	if ( is.null(genes)) {
		cat( "\nError: unable to deduce 'sampleID'. Not a 'SampleID' or a filename or a Cell Type...")
		return(NULL)
	}

	# now deduce what genes to look at
	if ( is.null( sortingGenes)) {
		sortingGenes <- c( "CD3E", "CD4", "CD8A", "CD11c", "CD14", "CD16", "CD19", "CD25", "CD27", 
				"CD34", "CD38", "CD45", "CD56",  "CD63",   
				"CD74", "CD123", "CD127", "CD161", 
				"CCR4", "CCR6", "CCR7", "CXCR3", "CXCR5", 
				"HLA-DRA", "IGHD", "TRAV1-2", "TRD@", "TRG@", "TRDV2")
	}
	givenGenes <- sortingGenes
	sortingGenes <- alias2Gene( sortingGenes, speciesID="Hs_grc")
	if ( speciesID != "Hs_grc") {
		sortingGenes <- ortholog( sortingGenes, "Hs_grc", speciesID)
	}
	showNames <- paste( givenGenes, " (", sortingGenes, ")", sep="")
	noAlias <- which( givenGenes == sortingGenes)
	showNames[noAlias] <- sortingGenes[noAlias]
	drops <- which( sortingGenes == "")
	if ( length(drops)) {
		sortingGenes <- sortingGenes[ -drops]
		showNames <- showNames[ -drops]
	}
	NG <- length( sortingGenes)

	# OK, we are ready to gather and show those genes
	intenV <- rep.int( 0, NG)
	where <- match( sortingGenes, genes, nomatch=0)
	intenV[ where > 0] <- inten[where]
	names(intenV) <- showNames
	if ( ! is.null( maxYshow)) {
		intenV[ intenV > maxYshow] <- maxYshow
	}
	bigY <- max( intenV, 5)
	yRange <- c( 0, bigY)

	yLabel <- paste( "Gene Expression (", intensityColumn, ")", sep="")
	log <- ""
	if ( useLog) {
		intenV[ intenV < 0.02] <- 0.02
		yLabel <- paste( "Gene Expression (", intensityColumn, "  Log scale)", sep="")
		log <- "y"
		yRange <- c( 0.02, bigY)
	}
	barplot( intenV, col=col, main=paste( "Cell Sorting Genes:  ", label), xlab=NA, 
			ylab=yLabel, log=log, ylim=yRange, las=3)

	lines( c(-10,100), c( 0.5,0.5), col='grey30', lwd=1, lty=2)
	dev.flush()

	return( invisible( intenV))
}
