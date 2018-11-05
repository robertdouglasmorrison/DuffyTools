# humanIDterms.R - turn human GENE_IDs into two extra columns ENTREZ_ID and common NAME


`getHumanIDterms` <- function( geneIDs) {

	geneIDs <- as.character( geneIDs)

	# format is {commonname:GInumber:chromosome:location}
	ginum <- sub( "(^.+:GI)([0-9]+)(:?.*$)", "\\2", geneIDs)
	nam <- sub( "(^.+?)(:.*$)", "\\1", geneIDs)

	# verify the Entrez ID is valid...
	suppressWarnings( ginum[ is.na( as.integer( ginum))] <- "" )

	out <- data.frame( "GENE_NAME"=nam, "ENTREZ_ID"=ginum)
	return( out)
}


`addHumanIDterms` <- function( mydf, idColumn="GENE_ID", aliasesToo=FALSE, speciesID=getCurrentSpecies()) {

	if ( ! idColumn %in% colnames(mydf)) {
		cat( "\nGeneID column not found: ", idColumn, "\nFound: ", colnames(mydf))
		return( mydf)
	}

	humanTerms <- getHumanIDterms( mydf[[ idColumn]])

	if ( aliasesToo) { 
		on.exit( setCurrentSpecies( getCurrentSpecies()))
		setCurrentSpecies( speciesID)
		humanTerms$ALIASES <- gene2Alias( humanTerms$GENE_NAME)

		# special case:  if given only the gene names, then the names are the IDs, and GI terms are all empty
		if ( all( humanTerms$ENTREZ_ID == "")) humanTerms <- humanTerms[ , -grep("ENTREZ_ID",colnames(humanTerms)), drop=F]
		if ( all( humanTerms$GENE_NAME == mydf[[idColumn]])) humanTerms <- humanTerms[ , -grep("GENE_NAME",colnames(humanTerms)), drop=F]
	}

	out <- cbind( humanTerms, mydf, stringsAsFactors=FALSE)
	return( out)
}


`addNameIDterms` <- function( mydf, idColumn="GENE_ID", nameColumn="GENE_NAME") {

	if ( ! idColumn %in% colnames(mydf)) {
		cat( "\nGeneID column not found: ", idColumn, "\nFound: ", colnames(mydf))
		return( mydf)
	}

	gmap <- getCurrentGeneMap()
	if ( ! "NAME" %in% colnames(gmap)) {
		cat( "\n'NAME' column not found in current species map..")
		return( mydf)
	}

	nameTerms <- rep.int( "", nrow(mydf))
	where <- match( mydf[[idColumn]], gmap$GENE_ID, nomatch=0)
	nameTerms[ where > 0] <- gmap$NAME[ where]

	# allow catching of Alternate Gene Map naming...
	stillBlank <- which( nameTerms == "")
	where2 <- match( sub( "::.+", "", mydf[[idColumn]]), gmap$GENE_ID, nomatch=0)
	use <- which( where == 0 & where2 > 0)
	nameTerms[ use] <- gmap$NAME[ where2[use]]

	if ( nameColumn %in% colnames(mydf)) {
		out <- mydf
		out[[ nameColumn]] <- nameTerms
	} else {
		out <- cbind( "GENE_NAME"=nameTerms, mydf, stringsAsFactors=FALSE)
		colnames(out)[1] <- nameColumn
	}
	return( out)
}


`addAltGeneIdentifierColumn` <- function( mydf, idColumn="GENE_ID", altColumn="ALT_GENE_ID") {

	# Alt Gene Maps use a compound name like GENE::extra
	# turn this into 2 columns for things that need the GENE_ID to be a real gene ID.
	if ( ! idColumn %in% colnames(mydf)) {
		cat( "\nGeneID column not found.  Expected: ", idColumn, "\nFound: ", colnames(mydf))
		return( mydf)
	}

	nameIn <- mydf[[idColumn]]
	nameOut <- sub( "::.+", "", nameIn)

	# let's put the new field right after the gene field
	whereGID <- match( idColumn, colnames(mydf))
	leftSide <- mydf[ , 1:whereGID]
	rightSide <- mydf[ , (whereGID+1):ncol(mydf)]

	out <- cbind( leftSide, "ALT_GENE_ID"=nameIn, rightSide, stringsAsFactors=FALSE)

	# put the shorter read GeneID back in where it goes
	out[[idColumn]] <- nameOut
	colnames(out)[ whereGID+1] <- altColumn
					
	if ( "ORIG_ID" %in% colnames(out)) out$ORIG_ID <- gene2OrigID( nameOut)

	# we are done
	out
}
