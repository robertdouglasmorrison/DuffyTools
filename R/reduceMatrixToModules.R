# reducedMatrix.R -- condense down a matrix of "Gene X Sample" to "module X Trait"
#			as of Dec 2023, removing all 'baseline adjustment' steps. Now handled
#			in previous steps when turning expression values to Mvalues

reduceMatrixToModules <- function( m, geneModules, sampleTraits, gene.names=rownames(m),
				sample.names=colnames(m), average.FUN=mean, checkGenes=TRUE) {

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

	for ( j in 1:NC) {
		myJ <- Cptrs[[ j]]
		lapply( 1:NR, function(i) {
			myI <- Rptrs[[ i]]
			v <- as.vector( m[ myI, myJ])
			outN[ i, j] <<- n <- sum( ! is.na(v))
			outV[ i, j] <<- average.FUN( v, na.rm=T)
			if( n > 1) {
				if ( diff( range( v, na.rm=T)) < 0.001) {
					pval <- 1
				} else {
					pval <- suppressWarnings( t.test( v)$p.value)
					if ( is.null(pval) || is.na( pval) || is.nan(pval)) pval <- 1
					# try to adjust for size of V
					pval <- pval * n
				}
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

	# given all the M values and P-values, make those PI values
	for ( j in 1:ncol(outV)) {
		outPI[ , j] <- piValue( outV[,j], outP[,j])
	}

	return( list( "moduleNames"=outNames, "matrix"=outV, "p.value"=outP, "pi.value"=outPI, "geneCounts"=nGenes))
}

