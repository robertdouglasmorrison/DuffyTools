# statsTools.R -- various statistics functions


`cv` <- function( x, y=NULL, na.rm=FALSE) {


	# do either the 1-sample or 2-sample Coefficient fo Variance
	if ( is.null( y)) {

		# simple case
		return( sd( x, na.rm=na.rm) / mean( x, na.rm=na.rm))
	}

	# the 2-sample CV:  using the 'within subject standard deviation method (Jones & Payne, 1997)
	if ( length(y) != length(x)) stop( "X and Y must be of equal length")

	twoN <- ( 2 * length(x))
	myMean <- sum( x + y) / twoN
	mySD <- sqrt( sum( (x-y)^2) / twoN)

	return( mySD / myMean)
}


`z.score` <- function(x, na.rm=FALSE) {
	return( (x - mean(x,na.rm=na.rm)) / sd(x,na.rm=na.rm))
}


`sparse.t.test` <- function( x, y=NULL, min.obs=3, ..., min.random.value=NULL) {

	# make extra observed values using sensible constraints
	myX <- x[ ! is.na(x)]
	meanX <- mean( myX)
	if ( length(myX) < min.obs) {
		sdX <- min( sqrt( abs( meanX/2)), abs( meanX/2))
		extraX <- rnorm( min.obs, mean=meanX, sd=sdX)
		if ( ! is.null( min.random.value)) extraX[ extraX < min.random.value] <- min.random.value
		x <- c( x, extraX)
	}
	if ( is.null(y)) {
		ans <- t.test( x, ...)
		# put the true estimate back in
		names(meanX) <- names(ans$estimate)
		ans$estimate <- meanX
		return( ans)
	}

	myY <- y[ ! is.na(y)]
	meanY <- mean( myY)
	if ( length(myY) < min.obs) {
		sdY <- min( sqrt( abs( meanY/2)), abs( meanY/2))
		extraY <- rnorm( min.obs, mean=meanY, sd=sdY)
		if ( ! is.null( min.random.value)) extraY[ extraY < min.random.value] <- min.random.value
		y <- c( y, extraY)
	}
	# put the true estimate back in
	ans <- t.test( x, y, ...)
	myMeans <- c( meanX, meanY)
	names(myMeans) <- names(ans$estimate)
	ans$estimate <- myMeans
	return( ans)
}


`two.means.difference` <- function( m1=0, m2=0, sd1=1, sd2=1, n1=5, n2=5,
					alternative=c("two.sided","less","greater")) {

	alternative <- match.arg( alternative)
	pvalue.factor <- if (alternative == "two.sided") 2 else 1
	LOWER.TAIL <- (alternative == "less")

	# calculate the difference in means when we don't have actual values as inputs to t.test
	mDiff <- m1 - m2
	sdTerm1 <- sd1^2 / n1
	sdTerm2 <- sd2^2 / n2
	sdJoint <- sqrt( sdTerm1 + sdTerm2)
	z <- mDiff / sdJoint
	pval <- min( pnorm( z, mean=0, sd=1, lower.tail=LOWER.TAIL) * pvalue.factor, 1.0)

	return( list( "mean.diff"=mDiff, "z"=z, "p.value"=pval))
}


`MAD.outliers` <- function( x, alpha=3, beta=1.4826) {

	med <- median( x, na.rm=T)
	abdev <- abs( x - med)
	mad <- beta * median( abdev)
	cutmin <- med - (mad * alpha)
	cutmax <- med + (mad * alpha)
	outliers <- sort( which( x < cutmin | x > cutmax))
	if ( length(outliers) && ! is.null(names(x))) names(outliers) <- names(x)[outliers]
	return( outliers)
}


`sd.pop` <- function(x) { sd(x)*sqrt((length(x)-1)/length(x))}


`biserial.cor.test` <- function( x, y, use = c("all.obs", "complete.obs"), y.level = 2) {

	# version of point bi-serial correlation
	if (!is.numeric(x)) stop("'x' must be a numeric variable.\n")
	y <- as.factor(y)
	if (length(levs <- levels(y)) != 2) stop("'y' must be a dichotomous variable.\n")
	if (length(x) != length(y)) stop("'x' and 'y' do not have the same length")
	use <- match.arg(use)
	if (use == "complete.obs") {
		cc.ind <- complete.cases(x, y)
		x <- x[cc.ind]
		y <- y[cc.ind]
	}
	ind <- (y == levs[y.level])
	N <- length(y)
	diff.mu <- mean(x[ind]) - mean(x[!ind])
	prob <- mean(ind)

	# biserial from 'ltm' package uses the sample SD instead of the Population SD
	#ans <- diff.mu * sqrt(prob * (1 - prob))/sd.pop(x)
	ans <- diff.mu * sqrt(prob * (1 - prob))/sd(x)

	t <- (ans*sqrt(N-2)) / sqrt( 1 - ans^2) 
	pv <- min( pt( q=t, df=N-2, lower.tail=F) * 2, 1.0)
	return( list( "estimate"=ans, "p.value"=pv))
}

