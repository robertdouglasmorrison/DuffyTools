# cellTypeTools.R -- tools for working with our immune cell subset datasets, as per Gene Sets, etc.
#			we will/have combinded human and mouse data

# circa 2020:   cloning all the LifeCycle data analysis & deconstruction tools to
#				apply them to mammalian immune cell type datasets.  See down below.


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


`cellTypeGeneHeatmap` <- function( m, nColors=50, nGenesPerCellType=c("25","50","100","250","500","1000"), 
					excludeCellTypePattern=NULL, includeCellTypePattern=NULL, 
					speciesID=getCurrentSpecies(), verbose=T, ...) {

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


# porting of Life Cycle tools below here...

`verifyCellTypeSetup` <- function() {

	isReady <- exists( "VectorSpace", envir=CellTypeEnv)
	isRightSpecies <- FALSE
	if ( isReady) {
		curSpecies <- get( "Species", envir=CellTypeEnv)
		if( ! is.null( curSpecies)) isRightSpecies <- ( curSpecies == getCurrentSpecies())
	}
	if ( !isReady || !isRightSpecies) CellTypeSetup( dataset="ImmuneCell.TargetMatrix", unitVectorMode="absolute")
	return()
}


`CellTypeSetup` <- function( dataset=c( "ImmuneCell.TargetMatrix", "Custom"), 
				unitVectorMode=c("absolute","relative","none"), 
				custom.file=NULL, custom.colors=NULL, preNormalize=TRUE, postNormalize=TRUE, 
				doRMA=FALSE, min.value=0.01, min.spread=2.0, verbose=FALSE) {

	dataset <- match.arg( dataset)
	unitVectorMode <- match.arg( unitVectorMode)
	ans <- NULL

	cat( "\nSetting up CellType dataset:  ", dataset)

	# get the data from the package
	if ( dataset == "ImmuneCell.TargetMatrix") {
	
		# make the dataset name from the current species
		prefix <- getCurrentSpeciesFilePrefix()
		dataFileName <- paste( prefix, dataset, sep=".")
		targetM <- targetColors <- NULL
		data( list=dataFileName, package="DuffyTools", envir=environment())
		if ( is.null( targetM)) {
			cat( "\nFailed to load Cell Type data object: ", dataFileName)
			return(NULL)
		}
		if ( exists( "immuneCellTargetColors")) targetColors <- immuneCellTargetColors

		ans <- buildCellTypeVectorSpace( file=NULL, tbl=targetM, 
			unitVectorMode=unitVectorMode, min.value=min.value, min.spread=min.spread, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=FALSE, verbose=verbose)
		ans2 <- buildCellTypeVectorSpace( file=NULL, tbl=targetM, 
			unitVectorMode="none", min.value=min.value, min.spread=NULL, 
			preNormalize=preNormalize, postNormalize=postNormalize,
			doRMA=FALSE, verbose=FALSE)
	}

	if ( dataset == "Custom") {

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

	CellTypeEnv[[ "VectorSpace" ]] <- ans
	CellTypeEnv[[ "IntensitySpace" ]] <- ans2
	CellTypeEnv[[ "Species" ]] <- getCurrentSpecies()
	if ( ! is.null( ans)) {
		CellTypeEnv[[ "N_STAGES" ]] <- ncol(ans) - 2   # geneID, Product come first
		CellTypeEnv[[ "STAGE_NAMES" ]] <- colnames(ans)[3:ncol(ans)]
		if ( ! is.null( targetColors)) CellTypeEnv[[ "STAGE_COLORS"]] <- targetColors
	}
	cat( "\nDone.\n")
}


`getCellTypeMatrix` <- function( mode=c("IntensitySpace", "VectorSpace")) {

	verifyCellTypeSetup()
	mode <- match.arg( mode)

	tbl <- CellTypeEnv[[ mode]]
	genes <- tbl$GENE_ID

	# get just the values, not the IDs and Products
	m <- as.matrix( tbl[ ,3:ncol(tbl)])
	rownames(m) <- genes
	return(m)
}


`buildCellTypeVectorSpace` <- function( file="", tbl=NULL,
					unitVectorMode=c( "absolute", "relative", "none"), 
					min.value=0.01, min.spread=2, preNormalize=FALSE, 
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
 

`calcCellTypeProfile` <- function( genes, inten, dropGenes=vector()) {

	verifyCellTypeSetup()

	# get the cell type data...  there is 'geneID, product', and then the intensities...
	vectorSpace <- CellTypeEnv[[ "VectorSpace"]]
	unitVectors <- vectorSpace[ , 3:ncol(vectorSpace)]
	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build up a tally of how much intensity goes to each cell type
	# do it for all genes, and the make stage histogram from the marker genes
	genes <- shortGeneName( genes, keep=1)
	# allow omitting some GeneIDs...
	if ( length( dropGenes)) {
		droppers <- shortGeneName( dropGenes, keep=1)
		toDrop <- which( genes %in% droppers)
		if ( length( toDrop)) {
			genes <- genes[ -toDrop]
			inten <- inten[ -toDrop]
		}
	}

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
	ans <- allSums * 100 / bigSum

	out <- list( "Profile"=ans, "IntensityVectors"=myIntens)
	return( out)
}


`calcCellTypeProfileFromFile` <- function( f, geneColumn="GENE_ID", intensityColumn="RPKM_M", 
					dropGenes=vector(), sep="\t") {

	verifyCellTypeSetup()

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

	return( calcCellTypeProfile( gset, inten, dropGenes=dropGenes))
}


`plotCellTypeProfileFromFileSet` <- function( fnames, fids, fcolors=NULL, geneColumn="GENE_ID", 
		intensityColumn="RPKM_M", yMax=NULL, legend.cex=0.8, max.labels=20, mask.low.pct=NULL,
		dropGenes=vector(), label="your label goes here...", sep="\t") {

	verifyCellTypeSetup()

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
				intensityColumn=intensityColumn, dropGenes=dropGenes, sep=sep)
		m[ i, ] <- ans$Profile
	}
	cat( "\n")

	# plot it
	plotCellTypeProfiles(m, col=fcolors, label=label, yMax=yMax, 
						legend.cex=legend.cex, max.labels=max.labels, mask.low.pct=mask.low.pct)

	return( m)
}


`plotCellTypeProfileFromMatrix` <- function( geneSet, intenMatrix, fids=colnames(intenMatrix), 
			fcolors=NULL, yMax=NULL, legend.cex=0.8, max.labels=20, mask.low.pct=NULL,
			dropGenes=vector(), label="your label goes here...") {

	verifyCellTypeSetup()

	N_STAGES <- CellTypeEnv[[ "N_STAGES"]]
	STAGE_NAMES <- CellTypeEnv[[ "STAGE_NAMES"]]

	# build the storage
	nColumns <- ncol( intenMatrix)
	m <- matrix( nrow=nColumns, ncol=N_STAGES)
	colnames(m) <- STAGE_NAMES
	rownames(m) <- fids

	# load each dataset in turn
	for( i in 1:nColumns) {
		ans <- calcCellTypeProfile( geneSet, intenMatrix[ ,i], dropGenes=dropGenes)
		m[ i, ] <- ans$Profile
	}

	# plot it
	plotCellTypeProfiles(m, col=fcolors, label=label, yMax=yMax, 
						legend.cex=legend.cex, max.labels=max.labels, mask.low.pct=mask.low.pct)

	return( m)
}


`plotCellTypeProfiles` <- function( m, col=NULL, yMax=NULL, label="", 
				legend.cex=0.8, max.labels=20, mask.low.pct=NULL) {

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
		 	m <- m[ , -drops]
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
		if ( N*NC >= 160) border <- NA
	}

	barSpace <- c( 0, N/4)

	if ( is.null( yMax)) yMax <- max(m) * 1.15
	mainText <- paste( "Cell Type Expression Profile Plot:\n", label)

	mp <- barplot(m, beside=T, col=col, border=border, main=mainText, 
		ylab="Percent of Total Gene Intensity", 
		space=barSpace, las=las, font.lab=2, font.axis=2, cex.lab=1, cex.axis=1, cex.names=0.8,
		xaxs="i", ylim=c(0,yMax), xlim=c( 1, (N*1.2)*(ncol(m)*1.2)))
	
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
}


`plotCellTypeProfileUnitVectors` <- function( gSet, col=1, lwd=1, legend=NA, plot=TRUE, yMax=1,
				legend.cex=1, label="", annotate=FALSE) {

	verifyCellTypeSetup()

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
	if ( ! is.na( legend)) {
		legend( x=legend, legend=gSet[who], col=colUse[who], lwd=3, cex=legend.cex)
	}
	
	# send back useful info
	units <- vectorSpace[ where, ]
	out <- list( "unitVectors"=units, "x"=xLocation)
	return( out)
}


`plotCellTypeProfileIntensity` <- function( gSet, col=1, lwd=1, pt.cex=1, legend=NA, plot=TRUE, 
				legend.cex=1, label="", minYmin=1, minYmax=NULL, threshold=NULL, annotate=FALSE) {

	verifyCellTypeSetup()

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

