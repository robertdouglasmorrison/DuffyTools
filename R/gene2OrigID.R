# gene2OrigID.R

`gene2OrigID` <- function( gNames) {

	gNames <- as.character( gNames)

	# default behavior is to return the genes as is
	out <- gNames

	# make sure set up
	geneMap <- getCurrentGeneMap()
	if ( nrow( geneMap) < 1) {
		cat( "\nWarning:  no gene map is current...")
		return( out)
	}
	if ( ! ("ORIG_ID" %in% colnames(geneMap))) return( out)

	# find those genes, and extract out the OrigID
	where <- base::match( gNames, geneMap$GENE_ID, nomatch=0)
	found <- where > 0
	out[ found] <-  geneMap$ORIG_ID[ where]

	return(out)
}


`orig2GeneID` <- function( gNames) {

	gNames <- as.character( gNames)

	# default behavior is to return the genes as is
	out <- gNames

	# make sure set up
	geneMap <- getCurrentGeneMap()
	if ( nrow( geneMap) < 1) {
		cat( "\nWarning:  no gene map is current...")
		return( out)
	}
	if ( ! ("ORIG_ID" %in% colnames(geneMap))) return( out)

	# find those genes, and extract out the OrigID
	where <- base::match( gNames, geneMap$ORIG_ID, nomatch=0)
	found <- where > 0
	out[ found] <-  geneMap$GENE_ID[ where]

	return(out)
}


addOrigIDterms <- function( mydf, idColumn="GENE_ID") {

	if ( ! idColumn %in% colnames(mydf)) {
		cat( "\nParasite GeneID column not found: ", idColumn, "\nFound: ", colnames(mydf))
		return( mydf)
	}

	origID <- gene2OrigID( mydf[[ idColumn]])

	out <- cbind( "ORIG_ID"=origID, mydf, stringsAsFactors=FALSE)
	return( out)
}
