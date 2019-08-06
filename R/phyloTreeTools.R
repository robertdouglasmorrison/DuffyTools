# phyloTreeTools.R -- make phylo tree plots


`plotPhyloTree` <- function( seqs, seqNames=names(seqs), dm=NULL,
			label.offset=1, tree.font=1, col=NULL, rotate.tree=0, tree.type="u",
			main=NULL, ...) {


	ids <- seqNames

	require( Biostrings)
	require( ape)
	checkX11()
	par( mai=c(0.5,0.5,0.5,0.5))

	if ( is.null( dm)) {
		cat( "\nMeasuring inter-sequence distances (N=", length(seqs),")", sep="")
		names(seqs) <- ids
		dm <- stringDist( seqs)
		#dm <- as.matrix( dm)
	} else {
		# catch if we were given a full symmetric matrix
		if ( is.matrix(dm) && (nrow(dm) == ncol(dm))) {
			dm <- as.dist( dm)
		}
	}
	
	# catch/cleanup the label colors
	if ( is.null( col)) {
		tip.color <- rainbow(length(seqs),end=0.76)
	} else {
		tip.color <- col
	}
	if ( length(seqs) != length(tip.color)) tip.color <- rep( tip.color, length.out=length(seqs))

	# draw it
	plot.phylo( as.phylo( hclust( dm)), type=tree.type, lab4ut="a", label.offset=label.offset, 
			font=tree.font, tip.color=tip.color, rotate.tree=rotate.tree, ...)
	if ( ! is.null( main)) mtext( main, side=3, font=2, cex=1.1)

	return( invisible( dm))
}

