# geneUniverse.R -- specify a subset of genes for comparison/analysis tools

`as.GeneUniverse` <- function(genes) {

	# let the subset of genes be given by a variety of ways
	out <- NULL
	if ( is.null(genes)) return(out)

	# a data.frame may have a GeneID column
	if ( is.data.frame(genes) && "GENE_ID" %in% colnames(genes)) {
		out <- genes$GENE_ID
	} else if ( is.character( genes)) {
		out <- genes
	} else {
		stop( paste( "Unknown/unsupported 'genes' argument to as.GeneUniverse() function: ", typeof(genes)))
	}

	return( shortGeneName( genes, keep=1))
}

		
