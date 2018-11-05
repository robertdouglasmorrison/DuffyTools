# cellCodeTools.R -- functions to implement and use the CellCODE Surrogate Proportion Variables features


`getCellCODEtagData` <- function( speciesID=getCurrentSpecies(), targetName=NULL) {

	# load the 'standard' subsetting marker genes 'Tag Data' for a given species
	#  Cell Type for mammals;   Life Cycle for parasites;  etc.

	if ( speciesID != (prevID <- getCurrentSpecies())) {
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( prevID))
	}
	prefix <- getCurrentSpeciesFilePrefix()

	if ( is.null( targetName)) {
		if ( speciesID %in% MAMMAL_SPECIES) targetName <- "ImmuneCell"
		if ( speciesID %in% PARASITE_SPECIES) targetName <- "LifeCycle"
	}

	cellCodeFile <- paste( prefix, "CellCODE", targetName, "TagData", sep=".")
	tagM <- NULL
	data( list=cellCodeFile, envir=environment())
	if ( is.null( tagM)) {
		cat( "\nFailed to load CellCODE Tag Data object.  Tried: ", cellCodeFile, "\n")
		return(NULL)
	}

	return( tagM)
 }
 
 
 `buildCellCODEtagDataFromMatrix` <- function( targetM, speciesID=getCurrentSpecies(), ...) {

	# given a matrix of gene expression, that spans some set of features to be described,
	# as human cell types or parasite life cycle stage data, build a 'Tag Data' object 
	# using the CellCODE package default tool

	require( CellCODE)

	# we need the rownames to be the geneIDs
	if ( speciesID != (prevID <- getCurrentSpecies())) {
		setCurrentSpecies( speciesID)
		on.exit( setCurrentSpecies( prevID))
	}
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)

	if ( is.null( colnames( targetM))) {
		cat( "\nTarget Matrix must have column names that are CellTypes, LifeStages, etc.")
		return( NULL)
	}

	if ( is.null( rownames( targetM))) {
		cat( "\nTarget Matrix must have rownames that are genes.")
		return( NULL)
	}

	nFound <- length( which( rownames(targetM) %in% shortGeneName( gmap$GENE_ID, keep=1)))
	if ( nFound < (nrow(targetM)/2)) {
		cat( "\nMost Target Matrix gene names not found in Gene Map.")
		return( NULL)
	}

	out <- tagData( targetM, ...)
	return( out)
}


`calcCellCODE.SPVs` <- function( m, groups=colnames(m), tagM=getCellCODEtagData(), plot=FALSE, 
				method=c( "mixed", "raw", "residual", "SVA"), ...) {

	# submit a matrix of expression data to the Surrogate Proportion Varables tool

	if ( length(groups) != ncol(m)) {
		cat( "\nMust give a 'group' name for each matrix column")
		return( NULL)
	}

	# the SPV tool assumes replicates per group...  Catch if all are 1 each
	if ( length( unique( groups)) == ncol(m)) method <- "raw"

	tagNames <- rownames( tagM)
	gNames <- shortGeneName( rownames(m), keep=1)
	nBoth <- length( intersect( tagNames, gNames))
	if ( nBoth < length(tagNames)/4) {
		cat( "\nNot enough gene names in common.  Verify row names and current speciesID")
		return( NULL)
	}
		
	# things seem OK to proceed
	require( CellCODE)

	method <- match.arg( method)
	ans <- suppressMessages( getAllSPVs( m, grp=groups, dataTag=tagM, plot=plot, method=method, ...))
	rownames(ans) <- colnames(m)

	return( t( ans))
}


`SPVtoPctSPV` <- function( spv) {

	# given the SPV output that describes each dimension as a relative projection onto all samples,
	# with many negative values, transform the rows to have all positive entries that sum to unity
	spvNames <- rownames(spv)
	sampleIDs <- colnames(spv)
	NS <- ncol(spv)
	NV <- nrow(spv)

	spvRange <- t( apply( spv, 1, range, na.rm=T))
	spvDelta <- apply( spvRange, 1, diff)

	out <- spv

	# we don't want any sample to see 'none' of a dimension, so include some small floor percentage
	smallOffset <- 0.01 / NS

	# transform each row into positive values that sum to 1
	for ( i in 1:NV) {
		v <- spv[ i, ]
		v <- (v - spvRange[i,1]) + smallOffset
		out[ i, ] <- v / sum(v, na.rm=T)
	}
	out
}


`fileSet.CellCODE.Proportions` <- function( files, fids, groups=fids, tagM=getCellCODEtagData(), 
					geneColumn="GENE_ID", intensityColumn="RPKM_M", sep="\t", 
					method=c( "mixed", "raw", "residual", "SVA"), reduceByGroups=FALSE, 
					plot=FALSE, ...) {

	cat( "\nGathering transcriptomes..")
	m <- expressionFileSetToMatrix( files, fids, geneColumn=geneColumn, intensityColumn=intensityColumn,
					sep=sep)
	rownames(m) <- shortGeneName( rownames(m), keep=1)

	method <- match.arg( method)
	ans <- matrix.CellCODE.Proportions( m, groups=groups, tagM=tagM, method=method, plot=plot, 
			reduceByGroups=reduceByGroups, ...)

	return( ans)
}


`matrix.CellCODE.Proportions` <- function( m, groups=colnames(m), tagM=getCellCODEtagData(), plot=FALSE, 
				method=c( "mixed", "raw", "residual", "SVA"), reduceByGroups=FALSE, ...) {

	# submit a matrix of expression data to the Surrogate Proportion Varables tool
	# to estimate the relative proportions of each SPV dimension in each sample

	# 1)  Get the SPV's
	cat( "\nfinding CellCODE SPV's..  ")
	method <- match.arg( method)
	spv <- calcCellCODE.SPVs( m, groups=groups, tagM=tagM, plot=FALSE, method=method, ...)
	if ( is.null( spv)) stop( "Failed CellCODE SPV generation step..")

	# 2)  Transform those SPVs into relative percentages
	spv <- SPVtoPctSPV( spv)

	# 3)  Now partition the observed expression proportionally
	# result has same shape/size/names as the SPV data
	outPctM <- spv
	for ( i in 1:nrow(spv)) {
		thisSPV <- rownames(spv)[i]
		# grab the subset of sample expression that 'belongs' to this SPV dimension
		who <- which( tagM[ , i] > 0)
		myGenes <- rownames(tagM)[ who]
		where <- match( myGenes, rownames(m), nomatch=0)
		smlM <- m[ where, ]
		smlExpression <- sum( smlM, na.rm=T)
		outPctM[ i, ] <- spv[ i, ] * smlExpression
	}
	# now scale each column to sum to 100
	for (j in 1:ncol(outPctM)) {
		v <- outPctM[ , j]
		outPctM[ , j] <- v * 100 / sum(v)
	}

	if (reduceByGroups) {
		grpFac <- factor( groups)
		outPctM <- t( apply( outPctM, 1, function(x) tapply( x, grpFac, mean)))
	}


	if (plot) plotTranscriptProportions( outPctM)

	# Done
	outPctM
}
