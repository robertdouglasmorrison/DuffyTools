# transcriptBlend.R -- code to implement transcript deconvolution, by modeling sum of component 'dimensions'


# create one transcript (named gene intensities) from a blend of column weights of a matrix of reference transcriptomes
`transcriptBlend` <- function( m, wts) {

	NT <- ncol(m)
	NWT <- length( wts)
	if (NT != NWT) stop( "transcriptBlend:  size of matrix & weights disagree")

	wts[ is.na(wts)] <- 0
	wts[ wts < 0] <- 0

	mTmp <- m
	for ( i in 1:NT) mTmp[ , i] <- m[ , i] * wts[i]

	out <- apply( mTmp, 1, sum, na.rm=T)
	names(out) <- rownames(m)
	out
}


# small function to generate synthetic transcriptome from blend of files
`fileSet.syntheticTranscriptome` <- function( fset, sids, wts=rep.int( 1, length(fset)), noise=0.1, 
				geneColumn="GENE_ID", intensityColumn="RPKM_M", sep="\t", rescale=TRUE) {
	
	if ( ! all( found <- file.exists( fset))) {
		cat( "\nSome named transcript files not found: \n", fset[ !found])
		return( NULL)
	}

	m <- expressionFileSetToMatrix( fset, sids, geneColumn=geneColumn, intensityColumn=intensityColumn,
					sep=sep, verbose=F)
	NG <- nrow(m)
	NC <- ncol(m)

	# add random noise to the data
	if ( noise > 0) {
		for ( i in 1:NC) {
			v <- m[ ,i]
			myNoise <- runif( NG, (1-noise), (1+noise))
			m[ ,i] <- v * myNoise
		}
	}

	wts <- as.numeric( wts)
	if ( any( is.na(wts))) stop( "Blending weights must be numeric and not NA")
				

	# do that blend
	ans <- transcriptBlend( m, wts)
	ans[ is.na(ans)] <- 0
	ans[ ans < 0] <- 0

	if ( sum( wts) != 1 && rescale) {
		scaleFac <- 1 / sum(wts)
		ans <- ans * scaleFac
	}

	# OK, make our result file, using an input file as a guide
	oldDF <- read.delim( fset[1], as.is=T, sep=sep)

	gids <- names(ans)
	wh <- match( gids, oldDF[[geneColumn]], nomatch=0)
	prods <- rep.int( "", NG)
	prods[wh > 0] <- oldDF$PRODUCT[ wh]

	out <- data.frame( "GENE_ID"=gids, "PRODUCT"=prods, "RPKM_M"=as.numeric(ans), stringsAsFactors=F)
	ord <- order( out$RPKM_M, decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	if ( geneColumn != "GENE_ID") names(out)[1] <- geneColumn
	if ( intensityColumn != "RPKM_M") names(out)[3] <- intensityColumn

	out
}


# the NLS fit function:  find the optimal blend
`fit.transcriptBlend` <- function( x, m, geneColumn="GENE_ID", intensityColumn="RPKM_M", useLog=FALSE, 
				normalize=TRUE, minIntensity=0, arrayFloorIntensity=NULL, 
				dropLowVarianceGenes=NULL, 
				algorithm=c("port", "default", "plinear", "LM", "GenSA"), 
				startFractions=NULL, verbose=TRUE) {

	# argument can be a data frame or a file name
	if ( is.data.frame(x)) {
		tbl <- x
	} else {
		f <- x[1]
		if ( ! file.exists(f)) {
			cat( "\nTranscriptome file not found:  ", f)
			return(NULL)
		}
		tbl <- read.delim( f, as.is=T)
	}

	if ( ! all( c( geneColumn, intensityColumn) %in% colnames(tbl))) {
		cat( "\nSome needed column names not in transcriptome.   Needed: ", geneColumn, intensityColumn)
		cat( "\nFound:  ", colnames(tbl))
		stop()
	}

	# get gene IDs and make sure we see most in both
	genes <- tbl[[ geneColumn]]
	where <- match( genes, rownames(m), nomatch=0)

	# if not enough found, try reducing the target matrix gene names to the short form...
	if ( sum( where > 0) < length(genes)/2) {
		rownames2 <- shortGeneName( rownames(m), keep=1)
		genes2 <- shortGeneName( genes, keep=1)
		where2 <- match( genes2, rownames2, nomatch=0)
		if ( sum( where2 > 0) > sum(where > 0)) where <- where2
	}

	tblUse <- tbl[ where > 0, ]
	NG <- nrow(tblUse)
	if (NG < 2*ncol(m)) {
		cat( "\n\nError:  Not enough gene IDs in common.  Check speciesID, transcriptome, & target matrix.\n")
		stop( "Unable to continue..")
	}

	genesUse <- tblUse[[ geneColumn]]
	intenUse <- tblUse[[ intensityColumn]]
	mUse <- m[ where, ]
	NC <- ncol( mUse)

	# the given transcriptome may be a microarray with non-zero baseline...  If so, clip and shift it toward zero
	if ( ! is.null( arrayFloorIntensity)) {
		floor <- as.numeric( arrayFloorIntensity)
		intenUse <- pmax( rep.int(0,NG), intenUse - floor)
	}

	# we also need to not use genes that can't help with the fit, e.g. those that are zero
	drops <- union( which( intenUse < minIntensity), which( apply( mUse, MARGIN=1, max, na.rm=T) < minIntensity))
	if ( length( drops)) {
		intenUse <- intenUse[ -drops]
		mUse <- mUse[ -drops, ]
		genesUse <- genesUse[ -drops]
	}

	# we also need to drop any NA's, as they will be dropped from NLS and break our results creation...
	drops <- which( is.na( intenUse))
	if ( length( drops)) {
		intenUse <- intenUse[ -drops]
		mUse <- mUse[ -drops, ]
		genesUse <- genesUse[ -drops]
	}

	# we may be asked to drop low variance genes from the target, i.e. those that don't differ by dimension
	if ( ! is.null( dropLowVarianceGenes)) {
		pctDrop <- 0
		if ( is.logical(dropLowVarianceGenes)) pctDrop <- 0.25
		if ( is.numeric( dropLowVarianceGenes)) pctDrop <- dropLowVarianceGenes
		pctDrop <- max( 0, pctDrop)
		pctDrop <- min( 0.9, pctDrop)
		if ( pctDrop > 0) {
			# generate the CVs for each gene
			myM <- apply( mUse, 1, mean, na.rm=T)
			myM <- ifelse( myM < 1, 1, myM)
			mySD <- apply( mUse, 1, sd, na.rm=T)
			myCV <- mySD / myM
			ord <- order( myCV, decreasing=T)
			nKeep <- round( nrow(mUse) * (1.0 -pctDrop))
			whoToKeep <- sort( ord[1:nKeep])
			if (verbose) cat( "\nDropping some low variance target genes..  Keeping ", nKeep, "of", nrow(mUse))
			intenUse <- intenUse[ whoToKeep]
			mUse <- mUse[ whoToKeep, ]
			genesUse <- genesUse[ whoToKeep]
		}
	}

	# ready to do the fit
	if (verbose) {
		cat( "\nFitting", length(genesUse), "genes as a blend of", NC, "transcriptomes:\n")
		cat( colnames(mUse), "\n")
	}

	# log scale and/or normalize ??
	if ( useLog) {
		intenUse <- log2( intenUse + 1)
		mUse <- log2( mUse + 1)
	}
	if ( normalize) {
		# we may use Quankenbush (total intensity normalize)
		mSums <- apply( mUse, MARGIN=2, sum, na.rm=T)
		tSum <- sum( intenUse)
		medSum <- median( c( mSums, tSum))
		mScale <- medSum / mSums
		tScale <- medSum / tSum
		for (j in 1:ncol(mUse)) mUse[,j] <- mUse[ ,j] * mScale[j]
		intenUse <- intenUse * tScale
		# we may force use of quantile normalize
		#mPreNorm <- cbind( intenUse, mUse)
		#mPostNorm <- rankEquivIntensity( mPreNorm, verbose=FALSE)
		#intenUse <- mPostNorm[ ,1]
		#mUse <- mPostNorm[ ,2:ncol(mPostNorm)]
	}

	# set the controls used by NLS
	controlList <- nls.control( tol=1e-05, maxiter=300, minFactor=1/1024, warnOnly=T)

	# one call for any number of dimensions
	# starting guess is 1/K if we were not given anything explicit
	if ( is.null( startFractions)) {
		startFractions <- rep.int( 1/NC, NC)
	} else {
		if ( length(startFractions) != NC) stop( "'startFractions' length incompatible with target matrix")
	}
	starts <- list( "WTi"=startFractions)

	algorithm <- match.arg( algorithm)
	if (algorithm == "GenSA") require( GenSA)
	if (algorithm == "LM") require( minpack.lm)
	if (verbose) cat( "\nFitting Algorithm:  ", algorithm)

	if ( algorithm == "port") {
		lowerBound <- rep.int( 0, NC)
		upperBound <- rep.int( 100, NC)
		fitAns <- try( nls( intenUse ~ transcriptBlend( mUse, WTi), start=starts,
				control=controlList, algorithm="port", lower=lowerBound, 
				upper=upperBound))
	} else if (algorithm == "GenSA") {
		WTi <- startFractions
		lowerBound <- rep.int( 0, NC)
		upperBound <- rep.int( 100, NC)
		fitAns <- try( do.TranscriptBlend.GenSA( intenUse, mUse, WTi, start=starts,
				lower=lowerBound, upper=upperBound))
	} else if (algorithm == "LM") {
		lowerBound <- rep.int( 0, NC)
		upperBound <- rep.int( 100, NC)
		fitAns <- try( nlsLM( intenUse ~ transcriptBlend( mUse, WTi), start=starts,
				control=controlList, algorithm="LM"))  #, lower=lowerBound, 
				# upper=upperBound))
	} else {
		fitAns <- try( nls( intenUse ~ transcriptBlend( mUse, WTi), start=starts,
				control=controlList, algorithm=algorithm))
	}

	if ( class( fitAns) == "try-error") {
		cat( "\nFitting of transcriptome failed...")
		return(NULL)
	} 

	# the simulated annealling answer is all set, but the others need the summary done explicitly
	if ( algorithm == "GenSA") { 
		out <- fitAns
	} else {
		fractions <- coef( fitAns)
		names(fractions) <- colnames(mUse)
		resids <- residuals(fitAns)
		fits <- fitted(fitAns)
		obs <- intenUse
		names(obs) <- names(resids) <- names(fits) <- shortGeneName( genesUse, keep=1)
		rms <- sqrt( mean( resids^2))
		# let's do a few other stats...
		cor.ans <- cor.test( obs, fits)
		r2.pearson <- cor.ans$estimate ^ 2
		pv <- cor.ans$p.value
		# also do the more general coefficient of determination
		meanI <- mean( obs, na.rm=T)
		SStotal <- sum( (obs - meanI) ^ 2)
		SSresid <- sum( resids ^ 2)
		r2.cod <- 1.0 - (SSresid / SStotal)

		out <- list( "BestFit"=fractions, "Observed"=obs, "Residuals"=resids, 
				"RMS.Deviation"=rms, "R2.CoD"=r2.cod, "R2.Pearson"=r2.pearson, "Pvalue"=pv)
	}

	return( invisible(out))
}


`do.TranscriptBlend.GenSA` <- function( inten, m, wts, start, lower, upper, seed=NULL) {

	# wrapper to implement fitting transcriptome by Generalize Simulated Annealing

	# GenSA as implemented has a hard coded seed for RNG.  We want random behavior each time thru
	# GEnSA docs suggest using a negative value
	my.seed <- -(as.integer( Sys.time()))


	# the penalty function that GenSA will minimize
	genSA.intensity.residual <- function( wts, obs, m) {

		model <- transcriptBlend( m, wts)

		# use mean absolute difference
		#resid <- abs( obs - model)
		#return( mean( resid, na.rm=T))

		# switching to RSS
		resid2 <- (obs - model) ^ 2
		return( sum( resid2, na.rm=T))
	}


	# set up to call GenSA

	# say 'good enough' if we explain 95% of the intensity
	stopValue <- mean( inten, na.rm=T) * 0.01
	control.list <- list( "maxit"=5000, "threshold.stop"=stopValue, 
				"temperature"=6000, "smooth"=FALSE, "max.call"=10000000,
				"max.time"=100, "trace.mat"=TRUE, "seed"=my.seed)

	#make sure the starts are above zero
	cat( "\nDebug GenSA: starts: ", wts)
	wts[ wts < 0.001] <- 0.001

	ans <- GenSA( par=wts, lower=lower, upper=upper, fn=genSA.intensity.residual, 
			control=control.list, obs=inten, m=m)

	# extract the answers
	fractions <- ans$par
	GenSA.Trace <<- ans$trace.mat
	names( fractions) <- colnames(m)
	model <- transcriptBlend( m, fractions)
	resids <- inten - model
	rms <- sqrt( mean( resids^2))
	# let's do a few other stats...
	cor.ans <- cor.test( inten, model)
	r2.pearson <- cor.ans$estimate ^ 2
	pv <- cor.ans$p.value
	# also do the more general coefficient of determination
	meanI <- mean( inten, na.rm=T)
	SStotal <- sum( (inten - meanI) ^ 2)
	SSresid <- sum( resids ^ 2)
	r2.cod <- 1.0 - (SSresid / SStotal)

	out <- list( "BestFit"=fractions, "Observed"=inten, "Residuals"=resids, 
			"RMS.Deviation"=rms, "R2.CoD"=r2.cod, "R2.Pearson"=r2.pearson, "Pvalue"=pv)
	return( out)
}

