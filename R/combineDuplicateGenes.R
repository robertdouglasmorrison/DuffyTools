# conbineDuplicateGenes.R -- merge (average) rows with the same geneID


`combineDuplicateGenes` <- function( x, geneColumn="GENE_ID", intensityColumn=NULL, 
				notesColumn=NULL) {

	# get the factoring of all the genes
	genePtr <- match( geneColumn, colnames(x), nomatch=0)
	if (genePtr < 1) {
		cat( "\n'geneColumn' not found.  Tried: ", geneColumn, "  Found: ", colnames(x))
		return( x)
	}
	
	genesIn <- x[[ genePtr]]
	nGenesIn <- length( genesIn)
	geneFac <- factor( genesIn)
	nGenesOut <- nlevels(geneFac)
	if ( nGenesIn == nGenesOut) return( x)

	cat( "  Combining Duplicate Genes:  ", nGenesIn, " -> ", nGenesOut)

	keeperRows <- tapply( 1:nGenesIn, geneFac, function(y) y[1])

	# make the result from the first copy of each gene
	out <- x[ keeperRows, ]

	# visit each column and perhaps do an averaging...
	cat( "   Averaging columns:  ")
	for ( j in 1:ncol(x)) {
		v <- x[[j]]
		if ( is.character( v)) next
		if ( is.logical( v)) next
		if ( is.numeric(v)) {
			avg.FUN <- mean
			col.text.in <- colnames(x)[j]
			col.text <- toupper( gsub( "_", "", col.text.in))
			col.text <- gsub( ".", "", col.text, fixed=T)
			cat( "  ", col.text.in)
			if ( regexpr( "PVAL",col.text) > 0) avg.FUN <- p.combine   #logmean
			vv <- tapply( v, geneFac, function(y) if ( length(y) < 2) y[1] else avg.FUN(y))
			if ( is.integer(v)) vv <- round(vv)
			out[[j]] <- vv
			# note the TPM units must sum to 1 million.  Re-assert that explicitly
			if ( grepl( "TPM", colnames(x)[j])) {
				vv <- round( vv * 1000000 / sum(vv,na.rm=T), digits=4)
				out[[j]] <- vv
			}
		}
	}

	# show what we did if given a 'notes' column
	if ( ! is.null( notesColumn)) {
		notesPtr <- match( notesColumn, colnames(x), nomatch=0)
		if ( notesPtr > 0) {
			notesText <- tapply( x[[notesColumn]], geneFac, function(y) { 
						if (length(y) < 2) y[1] else 
							paste( "AVG(", paste(y,collapse=","),")",sep="")
					})
			out[[notesPtr]] <- notesText
		}
	}

	# perhaps resort into the given ordering
	if ( ! is.null( intensityColumn)) {
		intensityPtr <- match( intensityColumn, colnames(x), nomatch=0)
		if ( intensityPtr > 0) {
			cat( "\nReordering genes by: ", intensityColumn)
			ord <- order( out[[ intensityColumn]], decreasing=T)
			out <- out[ ord, ]
			rownames(out) <- 1:nrow(out)
		} else {
			cat( "\nIntensity column for reordering genes not found: ", intensityColumn)
		}
	}
	cat( "  Done.\n")

	return( out)
}

