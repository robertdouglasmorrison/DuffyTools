# phyloTreeTools.R -- make phylo tree plots


`plotPhyloTree` <- function( seqs=NULL, seqNames=names(seqs), dm=NULL,
			label.offset=1, tree.font=1, col=NULL, rotate.tree=0, tree.type="u",
			main=NULL, ...) {

	require( Biostrings)
	require( ape)
	checkX11()
	saveMAI <- par( "mai")
	on.exit( par( "mai"=saveMAI))
	par( mai=c(0.5,0.5,0.5,0.5))

	# one of 'sequences' or 'distance matrix' must be given
	if ( ! is.null(seqs)) {
		ids <- seqNames
		N <- length(seqs)
	} else if ( ! is.null( dm) && is.matrix(dm)) {
		ids <- colnames(dm)
		N <- ncol(dm)
	} else {
		cat( "\nError:  one of 'seqs' or 'dm' must be non-NULL")
		return(NULL)
	}

	if ( is.null( dm)) {
		cat( "\nMeasuring inter-sequence distances (N=", N, ")", sep="")
		names(seqs) <- ids
		dm <- stringDist( seqs)
	} else {
		# catch if we were given a full symmetric matrix
		if ( is.matrix(dm) && (nrow(dm) == ncol(dm))) {
			dm <- as.dist( dm)
		}
	}
	
	# catch/cleanup the label colors
	if ( is.null( col)) {
		tip.color <- rainbow( N, end=0.76)
	} else {
		tip.color <- col
	}
	if ( length(tip.color) != N) tip.color <- rep( tip.color, length.out=N)

	# draw it
	plot.phylo( as.phylo( hclust( dm)), type=tree.type, lab4ut="a", label.offset=label.offset, 
			font=tree.font, tip.color=tip.color, rotate.tree=rotate.tree, ...)
	if ( ! is.null( main)) mtext( main, side=3, font=2, cex=1.1)

	return( invisible( dm))
}

