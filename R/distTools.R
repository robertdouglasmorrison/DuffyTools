# distTools.R -- manipulate distance objects, as from hclust, etc.


`distToMatrix` <- function( dist, labels=NULL) {

	attr <- attributes(dist)
	N <- attr$Size
	hasDiag <- attr$Diag
	isUpper <- attr$Upper

	# build the storage
	m <- matrix( 0, nrow=N, ncol=N)
	if ( is.null( labels)) {
		colnames(m) <- rownames(m) <- 1:N
	} else {
		lbls <- rep( labels, length.out=N)
		colnames(m) <- rownames(m) <- lbls
	}

	# fill the storage
	k <- 0
	for ( i in 1:(N-1)) 
	for ( j in (i+1):N) {
		k <- k + 1
		m[ j, i] <- m[i,j] <- dist[k]
	}

	m
}


`closestDistMatch` <- function( m, n=3, decreasing=FALSE) {

	# given a distance matrix, return a data frame of the closest N entries
	nams <- rownames(m)
	N <- nrow(m)
	if ( N != ncol(m)) stop( "Expected a square symmetric distance matrix")
	if ( m[1,N] != m[N,1]) stop( "Expected a square symmetric distance matrix")

	# the diagonal is not a selectable choice
	for ( i in 1:N) m[i,i] <- NA
	
	distOut <- matrix( NA, nrow=n, ncol=N)
	nameOut <- matrix( "", nrow=n, ncol=N)
	for ( i in 1:N) {
		closest <- order( m[,i], decreasing=decreasing)[1:n]
		distOut[ ,i] <- m[ closest, i]
		nameOut[ ,i] <- nams[ closest]
	}
	
	out <- data.frame( "Name"=nams, stringsAsFactors=FALSE)
	for ( j in 1:n) {
		smallDF <- data.frame( "Name"=nameOut[j,], "Dist"=distOut[j, ], stringsAsFactors=FALSE)
		colnames(smallDF) <- paste( c("Match","Dist"), j, sep="_")
		out <- cbind( out, smallDF, stringsAsFactors=FALSE)
	}
	out
}

