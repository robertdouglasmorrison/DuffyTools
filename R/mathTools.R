# mathTools.R	-- assorted math functions


# simple moving average
# get an 'average' value for each overlapping subset of 'window' values
movingAverage <- function( x, window=9, at=names(x), by=1, FUN=mean.default, do.ends=TRUE) {

	# apply the smoothing function to a window of K adjacent points
	if ( window < 3 || window > length(x)) stop( "Window size must be in 3..length(x)")
	K <- floor((window-1)/2) * 2 + 1
	N <- length(x)

	froms <- 1:(N-K+1)
	tos <- froms + K - 1

	if ( do.ends) {
		tosE1 <- ceiling(K/2) : (tos[1]-1)
		fromsE1 <- rep( 1, times=length(tosE1))
		fromsE2 <- ( 1:floor(K/2)) + froms[length(froms)] 
		tosE2 <- rep( N, times=length( fromsE2))
		froms <- c( fromsE1, froms, fromsE2)
		tos <- c( tosE1, tos, tosE2)
	}

	if ( by > 0) {
		use <- seq( 1, length(froms), by=round(by))
		froms <- froms[ use]
		tos <- tos[use]
	}

	values <- base::mapply( froms, tos, MoreArgs=list( "x"=x, "FUN"=FUN), FUN=function( i,j,x,FUN) {
			return( FUN( x[i:j]))
		})

	if ( is.numeric( at)) {
		at <- base::mapply( froms, tos, MoreArgs=list( "x"=as.numeric(at), "FUN"=mean.default), FUN=function( i,j,x,FUN) {
			return( FUN( x[i:j]))
		})
	}

	names( values) <- at
	return( values)
}


# takes the result of 'movingAverage' and forces the elements to have uniform distance between points
unitIntervalMA <- function( x, step=1) {

	# get the locations from the names
	atIn <- as.numeric( names(x))
	yIn <- x

	# get the integer range of X
	atRange <- range( round( atIn))
	atOut <- seq.int( atRange[1], atRange[2], by=step)

	ans <-spline( x=atIn, y=yIn, xout=atOut)

	xOut <- ans$x
	yOut <- ans$y
	# note that the spline at the very end points is a bit suspect
	yOut[1] <- yIn[1]
	names( yOut) <- xOut

	if ( xOut[ length(xOut)] < atRange[2]) {
		N <- length( yOut) + 1
		yOut[N] <- yIn[ length(yIn)]
		names(yOut)[N] <- atRange[2]
	}

	return( yOut)
}
	


# mean of a set of 'not normally distributed' values, like microarray probe intensities,
# by using a log transform, then mean, then unlog the result...
logmean <- function( x, na.rm=FALSE) {

	# zeros and negative break log
	# trap at smallest positive seen
	useX <- x
	useX[ x < 0] <- NA
	if ( sum( usable <- (! is.na(useX))) < 1) return( NA)
	if ( all( useX[ usable] == 0)) return( 0)
	isZero <- which( useX == 0)
	if ( length( isZero)) useX[ isZero] <- min( c( useX[ -isZero], 1e-300), na.rm=T)
	logx <- log2( useX)
	m <- mean.default( logx, na.rm=na.rm)
	return( 2 ^ m)
}


logMeanPlusSD <- function( x, nSD=1, na.rm=FALSE) {

	x[ x == 0] <- 1e-300
	x <- x[ x > 0]
	if ( length(x) < 1) return(0)
	logx <- log2( x)
	m <- mean.default( logx, na.rm=na.rm)
	mysd <- sd( logx, na.rm=na.rm)
	return( 2 ^ ( m + ( nSD * mysd)))
}


logMedianPlusSD <- function( x, nSD=1, na.rm=FALSE) {

	x <- x[ x > 0]
	if ( length(x) < 1) return(0)
	logx <- log2( x)
	m <- median( logx, na.rm=na.rm)
	mysd <- sd( logx, na.rm=na.rm)
	return( 2 ^ ( m + ( nSD * mysd)))
}


# square of the mean of the sqrt of X,  (i.e. generalized mean with p = 0.5)
sqrtmean <- function( x, na.rm=FALSE) {

	m <- mean.default( sqrt(x), na.rm=na.rm)
	return(m * m)
}


## Tukey biweight described here: http://www.affymetrix.com/support/technical/whitepapers/sadd_whitepaper.pdf
## "Statistical Algorithms Description Document"

# Input: a vector of values
# Output: a vector of weights based on Tukey biweight algorithm in the same order as the input

tukey.biweight = function(x, c=5, e=0.001, na.rm=FALSE){

	# c -  Tuning constant (fudge factor?)
	# e -  Small factor to prevent division by zero

	n = length( x)
	M = median( x, na.rm=na.rm)
	dM = x - M
	S = median( abs( dM), na.rm=na.rm)
	tuneConstant = rep( (c*S+e), times=n)
	u = dM / tuneConstant
	w = rep( 0, times=n)
	wTerm = (1.0 - (u * u))
	wTerm = (wTerm * wTerm)
	idx = which( abs(u)<=1);
	w[ idx] = wTerm[ idx]
	return(w);
}


tukey.mean <- function( x,  c=5,  e=0.001, na.rm=FALSE, plot=FALSE) {


	# allow a bit wider window when the number of observations is
	# too small for median to be robust...
	if (length( x) < 5) c <- c * 2

	wts <- tukey.biweight( x, c, e, na.rm=na.rm)
	m <- weighted.mean( x, wts, na.rm=na.rm)

	if ( ! plot) return( m)

	ord <- order(x)
	a <- hist( x, main="Tukey's Biweight Mean", ylab="Histogram Frequency",
			xlab="Observed Values to be Mean Averaged", col='gray90')
	xShow <- x[ord]
	yShow <- wts[ord]
	yMax <- max( a$counts)
	yTick <- yMax * 0.15
	yAt <- seq( 0, yMax, length.out=5)
	axis( side=4, at=yAt, label=round( max(yShow)*yAt/yMax, digits=2))
	mtext( "Weights calculated by 'Tukey.biweight()'", side=4, line=2)
	yShow <- yShow * yMax
	lines( xShow, yShow, col=2, lwd=2)
	points( jitter(xShow), jitter( yShow), pch=19, col=1, cex=1.2)

	lines( c(m,m), c(0,yTick), col=4, lwd=3, lty=3)
	mreg <- mean(x, na.rm=T)
	lines( c(mreg,mreg), c(0,yTick), col=3, lwd=3, lty=3)
	pos <- c(2,4)
	if ( mreg < m) pos <- c(4,2)
	text( c( m, mreg), c(yTick,yTick), c( "Tukey Mean", "Regular Mean"), col=c(4,3), pos=pos)

	return( m)
}


tukey.logmean <- function( x,  c=5,  e=0.0001, na.rm=FALSE) {

	logx <- log( x, 2)
	logm <- tukey.mean( logx, c, e, na.rm=na.rm)
	m <- 2^logm
	return( m)
}


errorBar <- function( x, mode=c("se", "sd", "mad", "ci", "gm.ci"), average.FUN=mean, plot=TRUE, at=1, whisker=0.2, 
			horiz=FALSE, error.col=1, error.lty=1, error.lwd=1) {

	mn <- average.FUN( x, na.rm=T)
	mode <- match.arg(mode)
	se1 <- se2 <- s <- sd( x, na.rm=T)
	if ( mode == "se") {
		se1 <- se2 <- s / sqrt( length(x))
	}
	if ( mode == "mad") {
		mn <- median( x, na.rm=T)
		dx <- abs( x - mn)
		se1 <- se2 <- mad <- median( dx)
	}
	if ( mode == "ci") {
		tt <- t.test( x, na.rm=T)
		ci <- tt$conf.int
		# the confidence interval is absolute, so turn these back to 'relative to mean' values
		se1 <- mn - ci[1]
		se2 <- ci[2] - mn
	}
	if ( mode == "gm.ci") {
		tt <- t.test( log10(x), na.rm=T)
		ci <- 10 ^ tt$conf.int
		# the confidence interval is absolute, so turn these back to 'relative to mean' values
		se1 <- mn - ci[1]
		se2 <- ci[2] - mn
	}

	if (plot && horiz) {
		lo <- xlo <- mn - se1
		hi <- xhi <- mn + se2
		lines( y=c( at,at), x=c( xlo, xhi), col=error.col, lty=error.lty, lwd=error.lwd)
		if ( whisker > 0) {
			lines( y=c( at-whisker,at+whisker), x=c( xlo, xlo), col=error.col, lty=error.lty, lwd=error.lwd)
			lines( y=c( at-whisker,at+whisker), x=c( xhi, xhi), col=error.col, lty=error.lty, lwd=error.lwd)
		}
	} else {
		lo <- ylo <- mn - se1
		hi <- yhi <- mn + se2
		lines( c( at,at), c( ylo, yhi), col=error.col, lty=error.lty, lwd=error.lwd)
		if ( whisker > 0) {
			lines( c( at-whisker,at+whisker), c( ylo, ylo), col=error.col, lty=error.lty, lwd=error.lwd)
			lines( c( at-whisker,at+whisker), c( yhi, yhi), col=error.col, lty=error.lty, lwd=error.lwd)
		}
	}

	return( invisible( c("min"=lo, "mid"=mn, "max"=hi)))
}
		

averageLineWithErrorBars <- function( x, y, average.FUN=mean, col=1, error.col=1, error.lty=1, whisker=0.2, ...) {

	xfac <- factor( x)
	ptrs <- tapply( x, xfac, FUN=NULL)
	avgY <- tapply( y, xfac, FUN=average.FUN, na.rm=T)

	lines( levels(xfac), avgY, col=col, ...)

	out <- matrix( NA, 3, nlevels(xfac))
	colnames(out) <- levels(xfac)
	rownames(out) <- c( "min", "mid", "max")
	nOut <- 0

	tapply( ptrs, xfac, function(x) {
			at <- as.numeric( levels(xfac)[ x[1]])
			myYs <- y[ ptrs == x[1]]
			ans <- errorBar( myYs, "se", average.FUN=average.FUN, at=at, whisker=whisker, 
					error.col=error.col, error.lty=error.lty)
			nOut <<- nOut + 1
			out[ , nOut] <<- ans
			return(NULL)
		})

	return( invisible( out))
}

