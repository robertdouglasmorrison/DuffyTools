# significantGenes.R  -- get the subset of genes that are UP and Significant


`significantGenes` <- function( file, geneColumn="GENE_ID", foldColumn="LOG2FOLD", pvalueColumn="PVALUE", 
			keepIntergenics=FALSE, min.fold=1, max.pvalue=0.05, min.genes=NULL, 
			shortNames=TRUE, sep="\t", verbose=TRUE) {

	tbl <- read.delim( file, as.is=T, sep=sep)
	if (verbose) cat( "\nRead file: ", file, "\nN_Genes: ", nrow(tbl))

	if ( !( all( c( geneColumn, foldColumn, pvalueColumn) %in% colnames(tbl)))) {
		cat( "\nMissing columns:  looked for: ", geneColumn, foldColumn, pvalueColumn, 
				"\n  \tFound: ", colnames(tbl))
		return()
	}

	# extract the parts we want
	genes <- tbl[[ geneColumn]]
	if (shortNames) genes <- shortGeneName( genes, keep=1)
	fold <- as.numeric( tbl[[ foldColumn]])
	pval <- as.numeric( tbl[[ pvalueColumn]])
	if ( "PRODUCT" %in% colnames(tbl)) {
		prods <- as.character( tbl[[ grep( "PRODUCT", colnames(tbl))[1] ]])
	} else {
		prods <- NULL
	}

	# allow the removal of non genes, etc.
	drops <- vector()
	if ( ! keepIntergenics) {
		drops <- grep( "(ng)", genes, fixed=TRUE)
		if (verbose) cat( "\nDrop as non-genes: ", length(drops))
	}
	drops2 <- which( fold < min.fold)
	if (verbose) cat( "\nDrop due to fold-change: ", length(drops2))
	drops3 <- which( pval > max.pvalue)
	if (verbose) cat( "\nDrop due to P-value: ", length(drops3))
	drops <- sort( union( drops, union( drops2, drops3)))
	keep <- setdiff( 1:nrow(tbl), drops)
	if ( ! is.null( min.genes)) {
		keep <- union( 1:min.genes, keep)
		if (verbose) cat( "\nKeeping at least 'min.genes'=", min.genes)
	}
	genes <- genes[ keep]
	fold <- fold[ keep]
	pval <- pval[ keep]
	if ( is.null( prods)) {
		prods <- gene2Product( genes)
	} else {
		prods <- prods[ keep]
	}

	out <- data.frame( "GENE_ID"=genes, "PRODUCT"=prods, "LOG2FOLD"=fold, 
				"PVALUE"=pval, stringsAsFactors=F)

	if ( nrow(out)) rownames(out) <- 1:nrow(out)
	if (verbose) cat( "\nN_Significant_Genes: ", nrow(out))
	return( out)
}

