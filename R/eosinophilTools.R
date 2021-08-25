# eosinophilTools.R -- try to capture info from eosinophil marker genes



# extract expression values from a variety of object types
`eosinophilExpression` <- function( x, value.mode=c("absolute","relative"), sep="\t") {

	eosinMarkerGenes <- c(  "FRRS1", "ADORA3", "ALOX15", "CCL23", "AFF2", 
				"CEBPE", "AOC1", "TFF3", "PRSS33", "CCL15", 
				"KHDRBS3", "VSTM1", "IDO1", "CD24", "OLIG2")
	eosinRPKMvalues <- c(    61.4777, 100.8599, 284.4524, 94.6696, 20.9992, 
				128.5605, 24.3029, 32.6940, 86.2324, 8.9697, 
				6.6057, 128.8078, 206.893, 88.6310, 21.5762)

	isMAT <- isDF <- isVEC <- FALSE
	out <- NULL

	if ( is.character(x)) {
		x  <- read.delim( x[1], as.is=T, sep=sep)
	}

	if ( is.matrix(x)) {
		geneNames <- shortGeneName( rownames(x), keep=1)
		values <- x
		isMAT <- TRUE
	} else if ( is.data.frame( x)) {
		gidColumn <- which( toupper(colnames(x)) %in% c( "GENE_ID", "GENEID", "GENE"))[1]
		if ( is.na( gidColumn)) {
			cat( "\nUnable to find 'GENE_ID' column in data frame.  \nEncountered: ", colnames(x))
			return( NULL)
		}
		intenColumn <- which( toupper(colnames(x)) %in% c( "RPKM_M", "INTENSITY", "TPM_M"))[1]
		if ( is.na( intenColumn)) {
			cat( "\nUnable to find 'RPKM_M' column in data frame.  \nEncountered: ", colnames(x))
			return( NULL)
		}
		geneNames <- shortGeneName( x[[gidColumn]], keep=1) 
		values <- as.numeric( x[[intenColumn]])
		isDF <- TRUE
	} else {
		geneNames <- shortGeneName( names(x), keep=1)
		values <- as.numeric(x)
		isVEC <- TRUE
	}

	where <- match( eosinMarkerGenes, geneNames, nomatch=NA)
	if ( all( is.na( where))) {
		cat( "\nWarning:  failed to find any named eosinophil genes..")
		return( NULL)
	}

	value.mode <- match.arg( value.mode)
	if ( isVEC || isDF) {
		out <- values[ where]
		names(out) <- eosinMarkerGenes
		if ( value.mode == "relative") out <- out / eosinRPKMvalues
	} else if ( isMAT) {
		out <- values[ where, , drop=F]
		rownames(out) <- eosinMarkerGenes
		if ( value.mode == "relative") for ( j in 1:ncol(out)) out[,j] <- out[,j] / eosinRPKMvalues
	}
	return( out)
}

