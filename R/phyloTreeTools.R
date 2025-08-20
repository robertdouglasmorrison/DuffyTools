# phyloTreeTools.R -- make phylo tree plots


`plotPhyloTree` <- function( seqs=NULL, seqNames=names(seqs), dm=NULL, newick=NULL, 
			tree.font=1, col=NULL, rotate.tree=0, tree.type="u",
			main=NULL, label.offset=if (is.null(newick)) 1 else 0.01, ...) {

	require( Biostrings)
	if (version$major == "4" && as.numeric( version$minor) >= 4) require( pwalign)
	require( ape)
	checkX11()
	saveMAI <- par( "mai")
	on.exit( par( "mai"=saveMAI))
	par( mai=c(0.5,0.5,0.5,0.5))

	# one of 'sequences' or 'distance matrix' or a Newick format file/string must be given
	isNEWICK <- FALSE
	if ( ! is.null(seqs)) {
		ids <- seqNames
		N <- length(seqs)
	} else if ( ! is.null( dm) && is.matrix(dm)) {
		ids <- colnames(dm)
		N <- ncol(dm)
	} else if ( ! is.null( newick)) {
		isNEWICK <- TRUE
	} else {
		cat( "\nError:  one of 'seqs' or 'dm' or 'newick' must be non-NULL")
		return(NULL)
	}

	if ( is.null( dm) && ! is.null(seqs)) {
		cat( "\nMeasuring inter-sequence distances (N=", N, ")", sep="")
		names(seqs) <- ids
		dm <- stringDist( seqs)
	} else if ( ! is.null(dm)) {
		# catch if we were given a full symmetric matrix
		if ( is.matrix(dm) && (nrow(dm) == ncol(dm))) {
			dm <- as.dist( dm)
		}
	} else {
		# read in a Newick format tree
		newickAns <- NULL
		if ( is.character(newick) && file.exists(newick)) newickAns <- read.tree( file=newick)
		if ( is.character(newick) && grepl( "\\(.+\\)", newick)) newickAns <- read.tree( text=newick)
		if ( is.list(newick) && "tip.label" %in% names(newick)) newickAns <- newick
		if ( is.null(newickAns)) {
			cat( "\nError:  failed to load a Newick tree from arguement: ", newick)
			return(NULL)
		}
		# clean the numeric prefixes from the labels
		ids <- newickAns$tip.label <- sub( "^[0-9]+[_ ]", "", newickAns$tip.label) 
		N <- length(ids)
	}
	# catch/cleanup the label colors
	if ( is.null( col)) {
		tip.color <- rainbow( N, end=0.76)
	} else {
		tip.color <- col
	}
	if ( length(tip.color) != N) tip.color <- rep( tip.color, length.out=N)

	# draw it
	if ( isNEWICK) {
		plot.phylo( newickAns, type=tree.type, lab4ut="a", label.offset=label.offset, 
				font=tree.font, tip.color=tip.color, rotate.tree=rotate.tree, ...)
	} else {
		plot.phylo( as.phylo( hclust( dm)), type=tree.type, lab4ut="a", label.offset=label.offset, 
				font=tree.font, tip.color=tip.color, rotate.tree=rotate.tree, ...)
	}
	if ( ! is.null( main)) mtext( main, side=3, font=2, cex=1.1)

	return( invisible( dm))
}

