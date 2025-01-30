# aliasTools.R


`loadAliasTable` <- function() {

	curSpecies <- getCurrentSpecies()
	if ( (!exists( "AliasTable", envir=AliasEnv)) || (AliasEnv[[ "AliasSpecies"]] != curSpecies)) {
		AliasTable <- NULL
		toLoad <- paste( getCurrentSpeciesFilePrefix(), "AliasTable", sep=".")
		data( list=toLoad, envir=environment())
		if ( is.null( AliasTable)) {
			# before we complain, allow looking in the current path
			localFile <- paste( toLoad, "rda", sep=".")
			if ( file.exists(localFile)) {
				cat( "\nInfo: no Alias Table in package. Using local file:  ", localFile, "\n")
				who <- load( localFile)
			} else {
				cat( "\nWarn: Failed to find/load Alias Table:  ", toLoad, "\n")
				AliasEnv[[ "AliasTable"]] <- NULL
				AliasEnv[[ "AliasSpecies"]] <- ""
				return( data.frame())
			}
		}
		
		# prep it a bit...
		#AliasTable$Alias <- toupper( AliasTable$Alias)
		AliasTable$Version <- as.numeric( AliasTable$Version)
		ord <- order( AliasTable$Alias, -(AliasTable$Version))
		AliasTable <- AliasTable[ ord, ]
		rownames(AliasTable) <- 1:nrow(AliasTable)
		AliasTable$UpperAlias <- toupper( AliasTable$Alias)
		AliasEnv[[ "AliasTable"]] <- AliasTable
		AliasEnv[[ "AliasSpecies"]] <- curSpecies
	}

	aliasTable <- AliasEnv[[ "AliasTable"]]
	return( aliasTable)
}


`alias2Gene` <- function( genes, speciesID=getCurrentSpecies()) {

	if ( speciesID != (curSpeciesID <- getCurrentSpecies())) {
		setCurrentSpecies(speciesID)
		on.exit( setCurrentSpecies(curSpeciesID))
	}

	aliasTable <- loadAliasTable()

	N <- length( genes)
	genesOut <- genes
	if ( ! N) return( genesOut)

	genesOutUP <- toupper(genes)
	versionOut <- rep( 0, times=N)

	# don't allow changes if the gene is already in the current gene map
	geneMap <- getCurrentGeneMap()
	goodGenes <- shortGeneName( geneMap$GENE_ID, keep=1)
	noChange <- which( genesOut %in% goodGenes)

	repeat {
		where <- base::match( genesOut, aliasTable$Alias, nomatch=0)
		whereUp <- base::match( genesOutUP, aliasTable$UpperAlias, nomatch=0)
		hitsUpOnly <- which( where == 0 & whereUp > 0)
		if ( length( hitsUpOnly)) where[ hitsUpOnly] <- whereUp[ hitsUpOnly]
		# don't allow genes already good to get changed by an alias entry
		where[ noChange] <- 0
		hasEntry <- which( where > 0)
		if ( length( hasEntry) < 1) break
		newGene <- aliasTable$GeneID[ where]
		newVersion <- aliasTable$Version[ where]
		okToChange <- which( versionOut[hasEntry] < newVersion)
		if ( length( okToChange) < 1) break

		genesOut[ hasEntry[ okToChange]] <- newGene[ okToChange]
		genesOutUP[ hasEntry[ okToChange]] <- toupper( newGene[ okToChange])
		versionOut[ hasEntry[ okToChange]] <- newVersion[ okToChange]
	}

	# some aliases may have revision suffixes.  Strip and try one more time
	notYet <- setdiff( which( genesOut == genes), noChange)
	if ( length( notYet)) {
		hasSuffix <- grep( "\\.[1-9]$", genes[notYet])
		genes2 <- sub( "\\.[1-9]$", "", genes[notYet])
		where <- base::match( genes2, aliasTable$Alias, nomatch=0)
		if (any( where > 0)) genesOut[ notYet[ where > 0]] <- aliasTable$GeneID[ where]
	}

	# one more thing to try:  many GeneIDs have various special characters, that the given gene symobls may not have
	# try to test the given names against genes (not aliases)
	notYet <- setdiff( which( genesOut == genes), noChange)
	if ( length( notYet)) {
		uniqGenes <- unique.default( aliasTable$GeneID)
		testAliasGenes <- gsub( "[^A-Z0-9]", "", toupper(uniqGenes))
		testGenes <- gsub( "[^A-Z0-9]", "", genesOutUP[notYet])
		where <- base::match( testGenes, testAliasGenes, nomatch=0)
		if (any( where > 0)) genesOut[ notYet[ where > 0]] <- uniqGenes[ where]
	}

	# one more thing to try:  many aliases have various special characters, that can break the full equality test
	# try to turn the aliases in the table into just the alphanumeric core.
	notYet <- setdiff( which( genesOut == genes), noChange)
	if ( length( notYet)) {
		# the grep substitution on the full alias table may be super slow, can we only do those with a chance...
		firstLetters <- unique( substr( genes[notYet], 1, 1))
		myPattern <- paste( "^[", paste( firstLetters, collapse=""), "]", sep="")
		toTest <- grep( myPattern, aliasTable$Alias)
		if ( length( toTest)) {
			testAlias <- gsub( "[^A-z0-9]", "", aliasTable$Alias[toTest])
			testGenes <- gsub( "[^A-z0-9]", "", genesOutUP[notYet])
			where <- base::match( testGenes, toupper(testAlias), nomatch=0)
			if (any( where > 0)) genesOut[ notYet[ where > 0]] <- aliasTable$GeneID[ toTest[ where]]
		}
	}

	return( genesOut)
}


`gene2Alias` <- function( genes, coreSubsetOnly=TRUE) {

	aliasTable <- subset.data.frame( loadAliasTable(), GeneID %in% genes)
	aFac <- factor( aliasTable$GeneID)
	aLevels <- levels(aFac)

	N <- length( genes)
	genesOut <- rep.int( "", N)

	aliaiOut <- tapply( aliasTable$Alias, aFac, function(x) {
			aliai <- x
			aliai <- aliai[ ! is.na( aliai)]
			if ( ! length(aliai)) return( "")

			# we may want to prune out some overly generic aliases...  This is provisional
			if (coreSubsetOnly) {
				drops <- vector()
				drops <- c( drops, grep( "^CCDS", aliai))
				drops <- c( drops, grep( "^ENSG", aliai))
				drops <- c( drops, grep( "^OTTHUMG", aliai))
				drops <- c( drops, grep( "^GI[0-9]+$", aliai))
				drops <- c( drops, grep( "^GI:", aliai))
				drops <- c( drops, grep( "^HGNC:", aliai))
				drops <- c( drops, grep( "^HPRD:", aliai))
				drops <- c( drops, grep( "^Hs\\.", aliai))
				drops <- c( drops, grep( "^MGI:", aliai))
				drops <- c( drops, grep( "^RGD:", aliai))
				drops <- c( drops, grep( "^UniSTS:", aliai))
				drops <- c( drops, grep( "[0-9]{5}$", aliai))
				drops <- c( drops, grep( "^MGC[0-9]{4,6}$", aliai))
				drops <- c( drops, grep( "^KIAA[0-9]{4,6}$", aliai))
				drops <- c( drops, grep( "^Q[0-9]{1,3}[A-Z]{1,3}[0-9]{1,3}[A-Z]?[0-9]?$", aliai))

				drops <- sort( unique( drops))
				if ( length(drops)) aliai <- aliai[ -drops]
				if ( ! length(aliai)) return("")
			}
			return( paste( aliai, collapse=", "))
		})
	
	where <- match( genes, aLevels, nomatch=0)
	genesOut[ where > 0] <- aliaiOut[where]

	return( genesOut)
}


`aliasGREP` <- function( genes) {

	aliasTable <- loadAliasTable()
	aliasTable$UpperAlias <- toupper( aliasTable$Alias)

	N <- length( genes)
	genesOut <- rep.int( "", N)

	lapply( 1:N, function(i) {
			hitsA <- grep( genes[i], aliasTable$Alias)
			hitsAup <- grep( toupper(genes[i]), aliasTable$UpperAlias)
			hitsG <- grep( genes[i], aliasTable$GeneID)
			if ( length( c( hitsA, hitsAup, hitsG)) < 1) return()
			ansA <- aliasTable$Alias[ union(hitsA,hitsAup)]
			ansG <- aliasTable$GeneID[hitsG]
			genesOut[i] <<- paste( sort( unique( c( ansA, ansG))), collapse=", ")
		})

	return( genesOut)
}


`addGeneSymbolsToAliasTable` <- function( speciesID, version) {

	setCurrentSpecies( speciesID)
	gmap <- subset( getCurrentGeneMap(), REAL_G == TRUE)

	aliasTable <- loadAliasTable()
	if ( nrow(aliasTable) < 1) {
		cat( "\nFound no defined alias data for species:  ", speciesID)
		cat( "\nNo Gene Symbols added...")
		return()
	}

	# we will make sure new symbols don't overlap with any already in the table
	aliasTable$UpperAlias <- toupper( aliasTable$Alias)

	# extract the symbols from the gene product terms, ignoring those with no defined symbol
	gmap$SYMBOL <- extractGeneSymbolFromProduct( gmap$PRODUCT)
	drops <- which( gmap$SYMBOL == "")
	if ( length( drops)) gmap <- gmap[ -drops, ]
	cat( "\nFound ", nrow(gmap), " gene symbols")

	# we insist on unique symbols, so drop any that are too generic
	redundantSymbols <- unique( gmap$SYMBOL[ duplicated( gmap$SYMBOL)])
	drops <- which( gmap$SYMBOL %in% redundantSymbols)
	if ( length(drops)) {
		cat( "\nDropping for being non-unique: ", length(drops))
		gmap <- gmap[ -drops, ]
	}

	# also drop any that are already in the table
	drops <- which( gmap$SYMBOL %in% aliasTable$Alias)
	dropsUP <- which( toupper(gmap$SYMBOL) %in% aliasTable$UpperAlias)
	drops <- sort( union( drops, dropsUP))
	# lets be more precise...  drop if it already maps to the current geneID, but keep if it maps to older IDs
	newKey <- paste( gmap$GENE_ID[drops], toupper(gmap$SYMBOL[drops]), sep="::")
	aliasKey <- paste( aliasTable$GeneID, toupper(aliasTable$Alias), sep="::")
	reallyDrop <- which( newKey %in% aliasKey)
	drops <- drops[ reallyDrop]
	if ( length(drops)) {
		cat( "\nAlready in Alias Table: ", length(drops))
		gmap <- gmap[ -drops, ]
	}

	cat( "\nNew Symbol aliases to add:  ", nrow(gmap))
	if ( nrow(gmap) < 1) return()
	cat( "\nAssigning new aliases version number: ", version)
	oldVersion <- max( aliasTable$Version, na.rm=T)
	if ( version < oldVersion) {
		cat( "\nGive a new version number no less than highest previous version: ", oldVersion)
		return()
	}

	smlDF <- data.frame( "GeneID"=gmap$GENE_ID, "Alias"=gmap$SYMBOL, "Version"=version, stringsAsFactors=F)

	# strip our Upper case column off, and make sure we have just the columns/order we want
	aliasTable <- loadAliasTable()
	wantCols <- c( "GeneID", "Alias", "Version")
	aliasTable <- aliasTable[ , wantCols]

	# make the new bigger alias table
	AliasTable <- rbind( aliasTable, smlDF)
	outfile <- paste( getCurrentSpeciesFilePrefix(), "AliasTable.rda", sep=".")
	save( AliasTable, file=outfile)
	cat( "\nWrote new alias table:  :", outfile, "\n")

}
