# cellTypeTools.R -- tools for working with our immune cell subset datasets, as per Gene Sets, etc.
#			we will/have combinded human and mouse data


# clean up the cell type names a bit
`cleanCellTypeName` <- function( cellTypes) {

	cells <- cellTypes
	cells <- sub( "^B1ab$", "B1a/B1b B cells", cells)
	cells <- sub( "^Bcell$", "B cells", cells)
	cells <- sub( "^B.CD5$", "CD5+ B cells", cells)
	cells <- sub( "^B.Memory$", "Memory B cells", cells)
	cells <- sub( "^B.Mem$", "Memory B cells", cells)
	cells <- sub( "^B.Naive$", "Naive B cells", cells)
	cells <- sub( "^CD4_Tcell$", "CD4+ T cells", cells)
	cells <- sub( "^CD4Tcell$", "CD4+ T cells", cells)
	cells <- sub( "^CD4.Cntrl.Mem$", "CD4+ Central Memory T cells", cells)
	cells <- sub( "^CD4.CntrlMem$", "CD4+ Central Memory T cells", cells)
	cells <- sub( "^CD4.CenMem$", "CD4+ Central Memory T cells", cells)
	cells <- sub( "^CD4.Tfh$", "CD4+ T Follicular Helper cells (Tfh)", cells)
	cells <- sub( "^CD4.Tfr$", "CD4+ T Follicular Regulatory cells (Tfr)", cells)
	cells <- sub( "^CD4.cTfh$", "CD4+ T Follicular Helper cells (Tfh)", cells)
	cells <- sub( "^CD4.cTfr$", "CD4+ T Follicular Regulatory cells (Tfr)", cells)
	cells <- sub( "^CD4.Eff.Mem$", "CD4+ Effector Memory T cells", cells)
	cells <- sub( "^CD4.EffMem$", "CD4+ Effector Memory T cells", cells)
	cells <- sub( "^CD4.Naive$", "CD4+ Naive T cells", cells)
	cells <- sub( "^CD4.Th1$", "CD4+ Type 1 Helper T cells", cells)
	cells <- sub( "^CD4.Th2$", "CD4+ Type 2 Helper T cells", cells)
	cells <- sub( "^CD4.Th17$", "CD4+ Type 17 Helper T cells", cells)
	cells <- sub( "^CD4.Treg$", "CD4+ Regulatory T cells (Treg)", cells)
	cells <- sub( "^CD8_Tcell$", "CD8+ T cells", cells)
	cells <- sub( "^CD8Tcell$", "CD8+ T cells", cells)
	cells <- sub( "^CD8.Cntrl.Mem$", "CD8+ Central Memory T cells", cells)
	cells <- sub( "^CD8.CntrlMem$", "CD8+ Central Memory T cells", cells)
	cells <- sub( "^CD8.CenMem$", "CD8+ Central Memory T cells", cells)
	cells <- sub( "^CD8.CytoxCD62Lh$", "CD8+ Cytotoxic CD62L+ T cells", cells)
	cells <- sub( "^CD8.CytoxCD62Ll$", "CD8+ Cytotoxic CD62L- T cells", cells)
	cells <- sub( "^CD8.Eff.Mem$", "CD8+ Effector Memory T cells", cells)
	cells <- sub( "^CD8.EffMem$", "CD8+ Effector Memory T cells", cells)
	cells <- sub( "^CD8.EffMem.preTRM$", "CD8+ Effector (pre Tissue Resident Memory) T cells", cells)
	cells <- sub( "^CD8.Naive$", "CD8+ Naive T cells", cells)
	cells <- sub( "^Dendritic$", "Dendritic cells", cells)
	cells <- sub( "^GD.nonVD2$", "Gamma Delta non-VD2 T cells", cells)
	cells <- sub( "^GD.VD2$", "Gamma Delta VD2+ T cells", cells)
	cells <- sub( "^Macrophage$", "Macrophages", cells)
	cells <- sub( "^Monocyte$", "Monocytes", cells)
	cells <- sub( "^MyeloidDC$", "Myeloid dendritic cells", cells)
	cells <- sub( "^mDendritic$", "Myeloid dendritic cells", cells)
	cells <- sub( "^Neutrophil$", "Neutrophils", cells)
	cells <- sub( "^NKcell$", "NK cells", cells)
	cells <- sub( "^NK$", "NK cells", cells)
	cells <- sub( "^NKTcell$", "NKT cells", cells)
	cells <- sub( "^RBC$", "Red Blood cells", cells)
	cells <- sub( "^TCRgd$", "Gamma delta T cells", cells)
	cells <- sub( "^Teff$", "Effector T cells (Teff)", cells)
	cells <- sub( "^Tfh$", "Follicular helper T cells (Tfh)", cells)
	cells <- sub( "^Tnaive$", "Naive T cells", cells)
	cells <- sub( "^Treg$", "Regulatory T cells (Treg)", cells)
	cells <- sub( "^WholeBlood$", "Whole Blood", cells)
	cells
}


`gene2CellType` <- function( genes, speciesID=getCurrentSpecies()) {

	out <- rep.int( "", N <- length(genes))

	genesIn <- shortGeneName( genes, keep=1)
	prefix <- getOtherSpeciesFilePrefix( speciesID)
	geneCellTypes <- NULL
	data( list=paste( prefix,"GeneCellTypes",sep="."), package="DuffyTools", envir=environment())
	if ( is.null(geneCellTypes)) {
		cat( "\nWarning:  No 'GeneCellTypes' object loaded for species: ", speciesID)
		return( out)
	}

	where <- match( genesIn, geneCellTypes$GENE_ID, nomatch=0)
	out[ where > 0] <- geneCellTypes$CellType[ where]
	if ( all( out == "")) warning( paste( "No genes mapped to cell types.  Verify current SpeciesID"))

	out
}


`cellTypeEnrichment` <- function( cellTypes, mode=c("geneSets", "genes"), 
				upOnly=TRUE, minEnrich=1.25, maxPvalue=0.01, 
				wt.enrich=1, wt.pvalue=2, speciesID=getCurrentSpecies(), 
				geneUniverse=NULL, correct=TRUE, verbose=T) {
	
	# grab the universe of cell type data
	mode <- match.arg( mode)
	if ( mode == "geneSets") {
		geneSetCellTypes <- NULL
		data( list="GeneSetCellTypes", package="DuffyTools", envir=environment())
		if ( is.null(geneSetCellTypes)) stop( "'GeneSetCellTypes' object not loaded.")
		cellTypeUniverse <- geneSetCellTypes$CellType
	}
	if ( mode == "genes") {
		# speciesID is needed, since orthologging may be required
		if (is.null(speciesID)) stop( "'speciesID' must not be NULL for gene cell types.")
		prefix <- getOtherSpeciesFilePrefix( speciesID)
		geneCellTypes <- NULL
		data( list=paste( prefix,"GeneCellTypes",sep="."), package="DuffyTools", envir=environment())
		if ( is.null(geneCellTypes)) stop( "'GeneCellTypes' object not loaded.")

		# this table often has more than one cell type per gene, to show the diversity
		# limit it to the first/best cell type choice for each gene
		dups <- which( duplicated( geneCellTypes$GENE_ID))
		if ( length(dups)) geneCellTypes <- geneCellTypes[ -dups, ]

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

	# ready to calculate enrichment
	givenTbl <- table( cellTypes)
	givenPcts <- givenTbl * 100 / sum(givenTbl)
	givenNames <- names( givenTbl)
	allTbl <- table( cellTypeUniverse)
	allPcts <- allTbl * 100 / sum(allTbl)
	allExpect <- allPcts * length(cellTypes) / 100
	allNames <- names( allTbl)
	
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
			pvals[i] <- sum( allPs[1:(x+1)])
		} else {
			# the probability of X up to K genes
			pvals[i] <- sum( allPs[(x+1):length(allPs)])
		}
		if (verbose) cat( "\n", i, both[i], enrich[i], pvals[i])
	}
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
	cellTypes <- sort( unique( cids))
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


`cellTypeGeneHeatmap` <- function( m, nColors=50, nGenesPerCellType=c("100","250","500","1000"), excludeCellTypePattern=NULL,
					includeCellTypePattern=NULL, speciesID=getCurrentSpecies(), verbose=T, ...) {

	# given a matrix of gene expression, keep only the top N genes from the selected cell types,
	# and pass that subset on to the heatmap tool.
	if ( is.null( rownames(m))) stop( "Matrix must have genes as rownames.")
	rownames(m) <- shortGeneName( rownames(m), keep=1)

	# step 1:  get the cell type data, and down-select the gene sets to use
	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	dataSetName <- paste( getCurrentSpeciesFilePrefix(), "HumanImmuneSubsets", sep=".")
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

