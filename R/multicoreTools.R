# multicoreTools.R - functions to help leverage the 'multicore' R package

# as of R 2.15, uses package 'parallel' instead...   almost all things stay the same


`multicore.currentCoreCount` <- function() {

	return( as.integer( getOption( "cores", default=1)))
}


`multicore.totalCoreCount` <- function() {
	if ( is.null(getOption( "cores"))) return( multicore.setup())
	return( detectCores())
}


`multicore.setup` <- function( max.cores=NULL) {

	if ( "parallel" %in% rownames( installed.packages())) {
		require( "parallel")
		options( "cores"=detectCores() )
		if ( ! is.null( max.cores)) options( "cores"=min( max.cores, getOption("cores")))
		cat( "\nMulticore support:   cores=", getOption("cores"),"\n")
	} else {
		options( "cores"=1)
	}
	return( getOption( "cores"))
}



`multicore.seqID.order` <- function( seqIDs) {

	# some organisms have the biggest chromosomes first, others last
	# re-order to put the biggest first, so long running chromosomes start sooner
	if ( getCurrentSpecies() %in% c( MAMMAL_SPECIES, INSECT_SPECIES)) {
		return( seqIDs)
	} else {
		return( rev( seqIDs))
	}
}


`multicore.samplePairs` <- function( sampleSet, groupSet=NULL, symmetric.pairs=FALSE) {

	# flatten a set (list?) of samples into a long list of 2-sample pairs 
	# that is directly usable by 'mclapply' to parallelize 'pipe.DiffExpression()'

	out <- vector( mode="list")
	nout <- 0

	if ( ! is.list( sampleSet)) {
		tmp <- vector( mode="list", length=1)
		tmp[[1]] <- as.vector(sampleSet)
		sampleSet <- tmp
	} else {
		if ( ! is.null(groupSet)) {
			cat( "\n'groupSet' only valid for a vector of samples, not a list")
			groupSet <- NULL
		}
	}

	for ( i in 1:length( sampleSet)) {
		thisSet <- sampleSet[[i]]
		N <- length(thisSet)
		if ( N < 2) next
		if ( symmetric.pairs) {
		    for (j in 1:N) {
		    	for (k in 1:N) {
			    if (j == k) next
			    if ( ! is.null(groupSet)) {
			    	if (groupSet[j] == groupSet[k]) next
			    }
			    nout <- nout + 1
			    out[[ nout]] <- c(thisSet[j], thisSet[k])
			}
		    }
		} else {
		    for (j in 1:(N-1)) for (k in (j+1):N) {
			    if ( ! is.null(groupSet)) {
			    	if (groupSet[j] == groupSet[k]) next
			    }
			    nout <- nout + 1
			    out[[ nout]] <- c(thisSet[j], thisSet[k])
		    }
		}
	}

	return( out)
}


`multicore.tapply` <- function( x, INDEX, FUN, ...) {

	# do the subsetting of X into a list of vectors, to pass to 'pvec'

	nCores <- getOption("cores")
	if ( !is.null(nCores) && nCores > 1) {

		# get all the subsets
		setOfSets <- base::tapply( x, INDEX, FUN=function(x) {x}, simplify=FALSE)

		# transform to a list of non-empty chunks
		listOfSets <- vector( mode="list")
		Nin <- length(setOfSets)
		Nout <- 0
		ptr <- rep( 0, Nin)

		for ( i in 1:Nin) {
			thisvec <- unlist( setOfSets[ i])
			if ( !is.null(thisvec) && length(thisvec) > 0) {
				Nout <- Nout + 1
				listOfSets[[ Nout]] <- unlist(thisvec)
				ptr[i] <- Nout
			}
		}
		
		# call the parallel op
		ans <- mclapply( listOfSets, FUN=FUN, ..., 
				mc.cores=min( length(listOfSets), nCores),
				mc.preschedule=FALSE, mc.allow.recursive=FALSE)
		MCLAPPLY_DEBUG <<- ans

		out <- vector( mode="list", length=Nin)
		out[ ptr > 0] <- ans[ ptr]

		return( out)
	} else {
		return( base::tapply( x, INDEX, FUN, ...))
	}
}


`multicore.lapply` <- function( x, FUN, ..., preschedule=FALSE) {

	#if ( length(x) == 1) return( lapply( FUN(x, ...)))
	if ( length(x) == 1) return( FUN(x, ...))

	nCores <- getOption("cores")
	if ( !is.null(nCores) && nCores > 1) {

		# call the parallel op
		ans <- mclapply( x, FUN=FUN, ..., 
				mc.cores=min( length(x), nCores), mc.preschedule=preschedule,
				mc.allow.recursive=FALSE)
		MCLAPPLY_DEBUG <<- ans
		return( ans)
	} else {
		return( base::lapply( x, FUN, ...))
	}
}


`multicore.by` <- function( data, INDICES, FUN, ...) {

	# do the subsetting of X into a list of data frames

	nCores <- getOption("cores")
	if ( !is.null(nCores) && nCores > 1) {

		# get all the subsets
		rowPtrs <- base::tapply( 1:nrow(data), INDICES, FUN=NULL)

		# transform to a list of non-empty chunks
		listOfDFs <- vector( mode="list")
		Nin <- max( rowPtrs)
		Nout <- 0
		ptr <- rep( 0, Nin)
		for ( i in 1:Nin) {
			who <- which( rowPtrs == i)
			if (length(who) > 0) {
				Nout <- Nout + 1
				sml <- data[ who, ]
				listOfDFs[[ Nout]] <- sml
				ptr[ i ] <- Nout
			}
		}

		# call the parallel op
		ans <- mclapply( listOfDFs, FUN=FUN, ..., 
				mc.cores=min( length(listOfDFs), nCores),
				mc.preschedule=FALSE, mc.allow.recursive=FALSE)
		MCLAPPLY_DEBUG <<- ans

		out <- vector( mode="list", length=Nin)
		out[ ptr > 0] <- ans[ ptr]
		names(out) <- levels(INDICES)

		return( out)
	} else {
		return( by( data, INDICES, FUN, ...))
	}
}

