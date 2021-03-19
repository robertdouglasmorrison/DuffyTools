# gene2Product.R

`gene2Product` <- function( gNames, speciesID=getCurrentSpecies()) {

	# default behavior is to return empty character strings
	out <- rep( "", times=length( gNames))

	if ( speciesID != (curSpeciesID <- getCurrentSpecies())) {
		setCurrentSpecies(speciesID)
		on.exit( setCurrentSpecies(curSpeciesID))
	}

	# make sure set up
	geneMap <- getCurrentGeneMap()
	if ( nrow( geneMap) < 1) {
		cat( "\nWarning:  no gene map is current...")
		return( out)
	}

	# find those genes, and extract out the Product
	where <- base::match( gNames, geneMap$GENE_ID, nomatch=0)
	found <- where > 0
	out[ found] <-  geneMap$PRODUCT[ where]

	# if not yet found, try by common name too
	try2 <- which( !found)
	where <- base::match( gNames[try2], geneMap$NAME, nomatch=0)
	out[ try2[ where > 0]] <- geneMap$PRODUCT[ where]

	if ( "ORIG_ID" %in% colnames(geneMap)) {
		found <- (out != "")
		try3 <- which( !found)
		where <- base::match( gNames[try3], geneMap$ORIG_ID, nomatch=0)
		out[ try3[ where > 0]] <- geneMap$PRODUCT[ where]
	}

	return(out)
}


`gene2ProductAllSpecies` <- function( gNames, hints=NULL) {

	# default behavior is to return empty character strings
	NG <- length(gNames)
	out <- rep.int( "", NG)

	# do for all species
	saveSpecies <- getCurrentSpecies()
	allSpecies <- getCurrentTargetSpecies()

	# check against the VSA first, as they give 'fuller' descriptions for the VAR genes
	outvsa <- vsaGeneProduct( gNames)
	goodVSA <- which( outvsa != "")
	out[ goodVSA] <- outvsa[ goodVSA]
	found <- ( out != "")

	for ( spec in allSpecies) {
		setCurrentSpecies( spec)
		geneMap <- getCurrentGeneMap()

		# find those genes, and extract out the Product

		# pass 1: perfect hits...
		try1 <- which( !found)
		where <- base::match( gNames[try1], geneMap$GENE_ID, nomatch=0)
		out[ try1[ where > 0]] <-  geneMap$PRODUCT[ where]
		if ( all( found <- (out != ""))) break

		if ( "ORIG_ID" %in% colnames(geneMap)) {
			try2 <- which( !found)
			where <- base::match( gNames[try2], geneMap$ORIG_ID, nomatch=0)
			out[ try2[ where > 0]] <-  geneMap$PRODUCT[ where]
			if ( all( found <- (out != ""))) break
		}

		# pass 2:  human only short names
		if ( spec %in% MAMMAL_SPECIES) {
			try2b <- which( !found)
			gNamesTry <- shortGeneName( gNames[ try2b], keep=1)
			gNamesMap <- shortGeneName( geneMap$GENE_ID, keep=1)
			where <- base::match( gNamesTry, gNamesMap, nomatch=0)
			out[ try2b[ where > 0]] <-  geneMap$PRODUCT[ where]
			if ( all( found <- (out != ""))) break
		}

		# pass 3: partial matching
		try3 <- which( !found)
		where <- pmatch( gNames[try3], geneMap$GENE_ID, nomatch=0, duplicates.ok=TRUE)
		out[ try3[ where > 0]] <-  geneMap$PRODUCT[ where]
		if ( all( found <- (out != ""))) break

		# if not yet found, try by common name too
		try4 <- which( !found)
		where <- base::match( gNames[try4], geneMap$NAME, nomatch=0)
		out[ try4[ where > 0]] <- geneMap$PRODUCT[ where]
		if ( all( found <- (out != ""))) break
	}

	# those still blank could be JOSE's?
	stillBlank <- which( ! found)
	if ( length( stillBlank)) {
		outjos <- josGeneProduct( gNames[ stillBlank])
		out[ stillBlank] <- outjos
		stillBlank <- which( out == "")
	}

	# if given a set of hints, use that to fill in 'still blank' ones
	if ( length( stillBlank) > 0) {
	    if ( !is.null(hints) && (length(hints) == length( gNames))) {
		out[ stillBlank] <- hints[ stillBlank]
	    }
	}

	setCurrentSpecies( saveSpecies)

	return(out)
}


`extractGeneSymbolFromProduct` <- function( products) {

	out <- rep.int( "", N <- length(products))

	lapply( 1:N, function(x) {

			# the gene symbol, if present, is in parentheses at the end
			txt <- sub( " +$", "", products[x])
			nc <- nchar(txt)
			if ( substr(txt, nc,nc) != ")") return()
			# find the very last '('
			opens <- gregexpr( "(", txt, fixed=T)[[1]]
			openp <- opens[ length(opens)]
			if ( openp < 1) return()

			symbol <- substr( txt, openp+1, nc-1)
			out[x] <<- symbol
		})
	
	return(out)
}
