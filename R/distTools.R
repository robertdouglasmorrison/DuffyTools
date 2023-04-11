# distTools.R -- manipulate distance objects, as from hclust, etc.


`distToMatrix` <- function( dist, labels=NULL) {

	# turn a R 'dist' object (that is just one halve of the full N x N square)
	# into a full matrix

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


`fastADist` <- function( seqs) {

	# try to do a faster 'adist' call by only using one copy of each
	NS <- length( seqs)
	if ( !NS) return( matrix( 0, 0, 0))

	uSeqs <- unique.default( seqs)
	NU <- length( uSeqs)
	if ( NU == NS) return( adist(seqs))

	# get the distances of just the unique sequences
	udm <- adist( uSeqs)

	# now figure out how to reconstruct the full matix
	whU <- match( seqs, uSeqs)
	mOut <- matrix( 0, nrow=NS, ncol=NS)
	colnames(mOut) <- rownames(mOut) <- names(seqs)
	for ( i in 1:NS) {
		myU <- whU[i]
		myV <- udm[ myU, whU]
		mOut[ , i] <- myV
	}
	return( mOut)
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


`alignedSeqDist` <- function( x, delete.temp.files=TRUE) {

	# use MAFFT MSA alignment to call edit distance, not just 'adist' or 'stringDist'

	# could be given several types of argument
	fastaFile <- alnFile <- NULL
	makeFasta <- TRUE
	if ( is.character(x) && length(x) == 1 && grepl( "\\.fa(sta)?$", tolower(x))) {
		fastaFile <- x
		alnFile <- sub( "\\.fa(sta)?$", ".aln", x)
		makeFasta <- deleteFasta <- deleteAln <- FALSE
	} else if ( is.Fasta(x)) {
		seqs <- x$seq
		desc <- x$desc
	} else {
		seqs <- x
		desc <- names(x)
		if ( is.null(desc)) {
			desc <- 1:length(seqs)
		}
	} 

	# make a FASTA for MAFFT?
	if ( makeFasta) {
		fastaFile <- "Temp.AlignedSeqDist.fasta"
		alnFile <- "Temp.AlignedSeqDist.aln"
		file.delete( alnFile)
		writeFasta( as.Fasta( desc, seqs), fastaFile, line=100)
		deleteFasta <- deleteAln <- delete.temp.files
	} else {
		tmpFA <- loadFasta( fastaFile, verbose=F)
		desc <- tmpFA$desc
	}

	# do we do the MAFFT call
	if ( file.exists( alnFile)) {
		aln <- readALN( alnFile)
	} else {
		aln <- mafft( fastaFile, alnFile)
	}

	# extract the aligned matrxi, and calculate the distances
	aaM <- aln$alignment
	N <- nrow(aaM)
	NC <- ncol(aaM)
	dm <- matrix( 0, nrow=N, ncol=N)
	rownames(dm) <- colnames(dm) <- desc
	for ( i in 1:(N-1)) {
		ch1 <- aaM[ i, ]
		for (j in (i+1):N) {
			ch2 <- aaM[ j, ]
			dm[i,j] <- dm[j,i] <- sum( ch1 != ch2)
		}
	}

	# clean up
	if (deleteAln) file.delete( alnFile)
	if (deleteFasta) file.delete( fastaFile)

	return( dm)
}


`expressionDist` <- function(x) {

	# given a matrix of gene expression data, return a distance matrix of transcriptome similarity
	if ( ! is.matrix(x)) stop( "'expressionDist()' expects a matrix of gene expression values")
	N <- ncol(x)
	out <- matrix( 0, nrow=N, ncol=N)
	rownames(out) <- colnames(out) <- colnames(x)

	# use Euclidian Distance in N dimensions
	for ( i in 1:(N-1)) {
		x1 <- as.numeric( x[ ,i])
		for (j in (i+1):N) {
			x2 <- as.numeric( x[ ,j])
			dx <- x1 - x2
			dd <- sum( dx * dx, na.rm=T)
			out[ i, j] <- out[ j, i] <- round( sqrt(dd), digits=3)
		}
	}
	return(out)
}

