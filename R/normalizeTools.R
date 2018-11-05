# normalizeTools.R    - do RMA (Robust Multi-chip Average) and other normalization techniques



`duffyRMA` <- function( m, targetBGlevel=NULL, magnitudeScale=NULL, columnSets=list( "all"=1:ncol(m)),
			verbose=FALSE) {

	m_tmp <- m

	if ( ! is.null( magnitudeScale)) {
		globalSum <- nrow(m) * as.numeric( magnitudeScale)
		if (verbose) cat( "\nMagnitude Scaling to global intensity: ", globalSum)
		m_tmp <- duffyMagnitudeNormalize( m_tmp, globalSum, verbose=verbose)
	}

	# step 1: background subtraction
	m_tmp <- duffyRMA.bgSubtract( m_tmp, targetBGlevel=targetBGlevel, verbose=verbose)

	# step 2:  QN
	m_out <- duffyRMA.qn( m_tmp, columnSets=columnSets, verbose=verbose)

	return( m_out)
}


# simple global magnitude normalize to have all slides have a roughly similar scaling
`duffyMagnitudeNormalize` <- function( m, globalSum=sum(m,na.rm=T)/ncol(m), verbose=FALSE ) {

	if (verbose) cat( "\nGlobal Magnitude scaling...")
	mOut <- m
	avgValue <- globalSum
	for( i in 1:ncol(m)) {
		v <- m[ ,i]
		tmpV <- sum( v, na.rm=T)
		scaleFac <- avgValue / tmpV
		newV <- v * scaleFac
		mOut[ ,i] <- newV
	}
	return( mOut)
}



`duffyRMA.bgSubtract` <- function( m, targetBGlevel=NULL, verbose=FALSE) {

	nc <- ncol( m)
	nr <- nrow( m)

	if (verbose) cat( "\nRMA: background shift...   DKM")
	if ( ! is.null( targetBGlevel))  {
		if (verbose) cat( "   Preset =", targetBGlevel)
		target <- targetBGlevel
	} else {
		# get the current background level for each column
		# using "mode" of the density kernel
		curBG <- vector()
		for ( i in 1:nc) {
			curBG[ i] <- getDensityKernelMode( m[ , i])
		}
		medbg <- median( curBG, na.rm=TRUE)
		medbg <- 10 * round( medbg / 10)
		if (verbose) cat( "   Median =", medbg)
		target <- medbg
	}

	return( DKMshift( m, target=target))
}


`duffyRMA.qn` <- function( m, columnSets=list( "all"=1:ncol(m)), verbose=FALSE) {

	if (verbose) cat( "\nDoing RMA quantile normalize...")
	# do the quantile normalize...
	tmp <- rankEquivIntensity( m, columnSets=columnSets, verbose=verbose)
	
	return( tmp)
}


`getDensityKernelMode` <- function( v) {

	# mode of density kernal
	tmp <- density( v, kernel="epanechnikov", n=round( length(v)*0.05))
	big <- which.max( tmp$y)
	vMode <- tmp$x[ big]
	return( vMode)
}


`getGlobalMetrics` <- function(m, rows=1:nrow(m)) {

	# every column is a valid data set, so average them
	meanSet <- apply(m[rows,],2,mean.default, na.rm=TRUE)
	medSet <- apply(m[rows,],2,median, na.rm=TRUE)
	minSet <- apply(m[rows,],2,min, na.rm=TRUE)
	maxSet <- apply(m[rows,],2,max, na.rm=TRUE)
	sumSet <- apply(m[rows,],2,sum, na.rm=TRUE)
	myMean<- mean.default(meanSet)
	myMed <- mean.default(medSet)
	myMin <- mean.default(minSet)
	myMax <- mean.default(maxSet)
	mySum <- mean.default(sumSet)
	tmp <- log2(m[rows,])
	logmeanSet <- apply( tmp, 2, mean.default, na.rm=TRUE)
	myLogmean <- median( 2 ^ logmeanSet)

	dkmSet <- apply( m[rows,], 2, getDensityKernelMode)
	myDKM <- median( dkmSet)

	# baseline may need to be estimated...
	if ( exists( "spotMap") && nrow( m) == nrow( spotMap)) {
		myBaseline <- median( as.vector( m[ whoBaselineSpots(), ]))
	} else {
		# fall back is the mode of density.
		myBaseline <- myDKM
	}

	out <- list("mean"=myMean, "median"=myMed, "floor"=myMin, "max"=myMax, "sum"=mySum, "logmean"=myLogmean, 
			"baseline"=myBaseline, "dkm"=myDKM)

	return(out)
}


`globalNormalize` <- function( m, intNormMethod=c("sum", "median", "none", "logmean", 
				"baseline", "dkm", "rank"), rows=1:nrow(m), metricSet=NULL, target=NULL) {

	intNormMethod <- tolower( intNormMethod)
	intNormMethod <- match.arg( intNormMethod)

	# the target could be almost anything...
	if ( ! is.null( target)) {
		if ( typeof(target) == "character") {
			evalTarget <- eval( parse( text=base::paste( "as.integer( ", target, ")")))
			target2 <- range( evalTarget)[2]
			target <- range( evalTarget)[1]
		}
	}
	
	# normalize all columns of intensity to the same Total
	method <- match.arg(intNormMethod)
	if (method == "none") return(m)
	
	cat("\n  ScalingMethod:", intNormMethod)

	# make metrics from this data, if none given
	if (is.null(metricSet)) metricSet <- getGlobalMetrics(m, rows=rows)
		
	nspots <- nrow(m)

	# sum, each column gets scaled to same total summation
	if (method == "sum") {
		if ( is.null( target)) target <- metricSet$sum
		sums <- apply(m[rows,],2,sum,na.rm=TRUE)
		n <- length(sums)
		for (i in 1:n) {
			fac <- rep(target/sums[i],times=nspots)
			m[,i] <- m[,i] * fac
		}
		return(m)
	} else if (method == "median") {
		# median, each column gets scaled to have same median value
		if ( is.null( target)) target <- metricSet$median
		mids <- apply(m[rows,],2,median,na.rm=TRUE)
		n <- length(mids)
		for (i in 1:n) {
			fac <- rep(target/mids[i],times=nspots)
			m[,i] <- m[,i] * fac
		}
		return(m)
	} else if (method == "logmean") {
		# each column gets scaled to have same mean of its log2 values
		if ( is.null( target)) target <- metricSet$logmean
		log2Targ <- log2( target)
		tmp <- log2( m)
		tmpScaling <- tmp[rows,]
		l2means <- apply( tmpScaling, 2, mean.default, na.rm=TRUE)
		n <- length(l2means)
		for (i in 1:n) {
			fac <- rep(log2Targ/l2means[i],times=nspots)
			tmp[,i] <- tmp[,i] * fac
		}
		# undo the log step...
		m <- rep(2,times=(ncol(m)*nrow(m))) ^ tmp
		return(m)
	} else if (method == "baseline") {
		# baseline, each column gets a linear shift to have same median baseline value
		if ( is.null( target)) target <- metricSet$baseline
		if ( exists( "spotMap")) {
			who <- whoBaselineSpots()
		} else {
			who <- 1:nrow(m)
		}
		curBaseMedian <- apply( m[ who, ], 2, median, na.rm=TRUE)
		n <- length(curBaseMedian)
		globalFloor <- metricSet$floor
		for (i in 1:n) {
			fac <- rep( (target - curBaseMedian[i]), times=nspots)
			m[,i] <- m[,i] + fac
		}
		m <- clipMinIntensity( m, globalFloor)
		return(m)
	} else if (method == "dkm") {
		if ( is.null( target)) target <- metricSet$dkm
		m <- DKMshift( m, target)
		return(m)
	} else if (method == "rank") {
		# rank, each column gets scaled to a uniform non-parametric scale
		if ( is.null( target)) target <- metricSet$floor
		if ( target2 == target) target2 <- metricSet$maximum
		rnkValues <- seq( log(target,2), log(target2,2), length.out=nspots)
		rnkValues <- 2 ^ rnkValues
		n <- ncol(m)
		for (i in 1:n) {
			ord <- base::order( m[ ,i], na.last=FALSE)
			m[ ord,i] <- rnkValues
		}
		return(m)
	} else {
		# unknown...
		stop( paste( "unknown global normalize mode:  ", method, "\nKnown choices: ", intNormMethod))
	}
	return(NULL)
}


`clipMinIntensity` <- function( m, target=getGlobalMetrics(m)$floor) {

	# get the current mini9imum fro the metric set
	out <- ifelse(  m > target, m, target)
	return( out)
}


`DKMshift` <- function( m, target=getGlobalMetrics(m)$dkm) {

	NC <- ncol(m)
	mOut <- m

	targetDKM <- target

	# apply a 'smart' linear shift to each slide
	for ( i in 1:NC) {
		thisV <- m[ , i]
		thisDKM <- getDensityKernelMode( thisV)
		thisShift <- targetDKM - thisDKM

		if ( thisShift > 0) {
			# the slide needs to move right to higher magnitudes, so
			# use a linear shift (or multiplicitive scale)
			newV <- thisV + thisShift
			#scalefac <- targetDKM / thisDKM
			#newV <- thisV * scalefac

		} else {
			# negative shifs are much tougher...
			# we can linearly shift to the left, but the low values can go negative

			# step 1: linear shift
			tmpV <- thisV + thisShift

			# step 2: compress the lower tail to maintain its relative shape
			oldMin <- min( thisV)
			tmpMin <- min( tmpV)
			targetMin <- min( oldMin, targetDKM * 0.25)
			whoToFix <- which( tmpV < targetDKM)

			# make a scaled transformation to those pts to put them 
			# back into the range from "targetMin  to  targetDKM"
			maxCorrection <- targetMin - tmpMin
			maxDistance <- targetDKM - tmpMin

			myDistances <- targetDKM - tmpV[ whoToFix]

			# instead of linear, lets weight by distance away
			myFactors <- ( myDistances / maxDistance)
			myCorrections <- maxCorrection * myFactors
			tmpV[ whoToFix] <- tmpV[ whoToFix] + myCorrections
			newV <- tmpV
		}
		mOut[ , i] <- newV
	}
	return( mOut)
}



`rankEquivIntensity` <- function(m, columnSets=list( "all"=1:ncol(m)), verbose=TRUE) {

	# turn an object of multiple slide intensities
	# into a Rank Equivalent Intensity (REI) object
	nr <- nrow(m)
	nc <- ncol(m)

	# start with a couple copies
	newInt <- rankInt <- ranks <- m

	# first rank each spot in its slide, to build an ordered set of slides
	for (i in 1:nc) {
		ord <- base::order( m[ ,i], decreasing=TRUE, na.last=TRUE)
		ranks[ord ,i] <- 1:nr
		rankInt[ ,i] <- m[ord ,i]
	}

	# now do each group of columns
	for ( l in 1:length( columnSets)) {
	  
		useCols <- applyCols <- columnSets[[ l]]
		setName <- names( columnSets)[l]
		if (verbose) cat( "   ", setName, "(N=", length( useCols), ")...", sep="")
		if ( ! all( useCols %in% 1:nc)) {
			cat( "RMA error:  invalid slide column(s):   given: ", useCols, "   possible: ", 1:nc)
			next
		}
		if ( length( useCols) < 2) next

		# now average each row to get a REI for each rank
		rei <- apply( rankInt[ , useCols], 1, sqrtmean, na.rm=TRUE)
	
		# make sure they're still ordered
		rei <- base::sort( rei, decreasing=TRUE)

		# now build the new set of Rank Equiv Intensities
		for (i in applyCols) {
			newInt[ ,i] <- rei[ ranks[ ,i] ]
		}

	}

	# pass this back
	return( newInt)
}

