# reducedMatrix.R -- condense down a matrix of "Gene X Sample" to "module X Trait"

reduceMatrixToModules <- function( m, geneModules, sampleTraits, gene.names=rownames(m),
				sample.names=colnames(m), average.FUN=mean, 
				baselineTrait=NULL, checkGenes=TRUE) {

	if ( is.null( names(geneModules))) stop( "'geneModules' must be a named list")

	# make sure all row group elements are in the matrix
	Rnames <- gene.names
	if (checkGenes) {
		cat( "\nChecking modules for named genes..")
		Rgroups <- lapply( geneModules, function(x) unique( x[ x %in% Rnames]))
		names(Rgroups) <- names(geneModules)
		RgroupLens <- sapply( Rgroups, length)
		Rgroups <- Rgroups[ RgroupLens > 1]
	} else {
		Rgroups <- geneModules
	}
	RgroupNames <- names( Rgroups)
	NR <- length( Rgroups)
	if ( NR < 2) {
		stop( "Not enough Modules left...  Check 'getCurrentSpecies()'")
	}

	cat( "\nReducing ", nrow(m), " gene rows down to ", NR, " modules..")
	Rptrs <- lapply( Rgroups, function(x) match( x, Rnames))
	if ( any( is.na( unlist( Rptrs)))) stop( "Module gene names to matrix gene names mapping error..")

	# make sure all column group elements are in the colnames
	if ( ! is.list(sampleTraits)) {
		mysamples <- sampleTraits
		mynames <- names(sampleTraits)
		sampleTraits <- as.list( mysamples)
		names(sampleTraits) <- mynames
	}
	if ( is.null( names( sampleTraits))) names(sampleTraits) <- sapply( sampleTraits, function(x) x[1])

	Cnames <- sample.names
	Cgroups <- lapply( sampleTraits, function(x) x[ x %in% Cnames])
	names(Cgroups) <- names(sampleTraits)
	CgroupLens <- sapply( Cgroups, length)
	Cgroups <- Cgroups[ CgroupLens > 0]
	CgroupNames <- names( Cgroups)
	NC <- length( Cgroups)

	if ( NC < ncol(m)) cat( "\nReducing ", ncol(m), " samples down to ", NC, " sample traits.")
	if ( NC < 2) {
		stop( "Not enough sample trait groups left...  Check 'SampleIDs'")
	}
	Cptrs <- lapply( Cgroups, function(x) match( x, Cnames))

	# as of 2020, adding in a PI Value result as well
	outV <- outP <- outN <- outPI <- matrix( NA, nrow=NR, ncol=NC)
	colnames(outN) <- colnames(outV) <- colnames(outP) <- colnames(outPI) <- CgroupNames
	#rownames(outN) <- rownames(outV) <- rownames(outP) <- RgroupNames
	outNames <- RgroupNames
	nGenes <- sapply( Rptrs, length)

	# we may have been asked to force a trait to be the baseline, i.e. hardwire to zero
	if ( ! is.null( baselineTrait)) {
		baselineColumn <- match( baselineTrait, colnames(outV), nomatch=0)
		if ( baselineColumn == 0) {
			cat( "\nError in 'ReduceMatrixToModules()':  baseline trait not a valid choice.")
			cat( "\nGiven: ", baselineTrait, "  \tChoices: ", colnames(outV))
			stop()
		}
		cat( "\nLinear shift to fix as baseline trait:  ", baselineTrait)

		# visit the data ahead of time to generate the median of the baseline trait for each module
		myJ <- Cptrs[[ baselineColumn]]
		shiftV <- vector( length=length( Rptrs))
		lapply( 1:NR, function(i) {
			myI <- Rptrs[[ i]]
			v <- as.vector( m[ myI, myJ])
			shiftV[ i] <<- average.FUN( v, na.rm=T)
		})
	}

	for ( j in 1:NC) {
		myJ <- Cptrs[[ j]]
		lapply( 1:NR, function(i) {
			myI <- Rptrs[[ i]]
			v <- as.vector( m[ myI, myJ])
			if ( ! is.null( baselineTrait)) v <- v - shiftV[i]
			outN[ i, j] <<- n <- sum( ! is.na(v))
			outV[ i, j] <<- average.FUN( v, na.rm=T)
			if( n > 1) {
				if ( diff( range( v, na.rm=T)) < 0.01) v <- jitter(v)
				pval <- t.test( v)$p.value
				if ( is.null(pval) || is.na( pval) || is.nan(pval)) pval <- 1
				# try to adjust for size of V
				pval <- pval * n
			} else {
				pval <- 1
			}
			outP[ i, j] <<- pval
			return()
		})
	}

	outP[ is.na(outP)] <- 1
	outP[ outP > 1] <- 1

	# there may be rows with no observed data
	minN <- apply( outN, 1, min)
	drops <- which( minN <= 1)
	if( length( drops) > 0) {
		outV <- outV[ -drops, ]
		outP <- outP[ -drops, ]
		outPI <- outPI[ -drops, ]
		outN <- outN[ -drops, ]
		outNames <- outNames[ -drops]
		nGenes <- nGenes[ -drops]
	}
	if ( nrow( outV) < 2) stop( "Not eough modules have data...")

	# we may have been asked to force a trait to be the baseline, i.e. hardwire to zero
	# now done ahead of time, to be in effect when we measure P values

	# given all the M values and P-values, make those PI values
	for ( j in 1:ncol(outV)) {
		outPI[ , j] <- piValue( outV[,j], outP[,j])
	}

	return( list( "moduleNames"=outNames, "matrix"=outV, "p.value"=outP, "pi.value"=outPI, "geneCounts"=nGenes))
}
