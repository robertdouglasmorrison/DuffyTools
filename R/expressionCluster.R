# expressionCluster.R

# run a matrix of expression intensities  through a clustering tool

`expressionCluster` <- function( x, useLog=FALSE, normalize=FALSE, FUN=diana) {


	require( cluster)

	if ( useLog) {
		x <- log2( x + 1)
	}

	if (normalize) {
		sums <- apply( x, MARGIN=2, sum, na.rm=T)
		scaleFac <- median(sums) / sums
		for (j in 1:ncol(x)) {
			x[,j] <- x[,j] * scaleFac[j]
		}
	}

	clusterAns <- FUN( t(x), diss=F, metric="manhattan", stand=F, keep.diss=F, keep.data=F)

	return( clusterAns)
}

