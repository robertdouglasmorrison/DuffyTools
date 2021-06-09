# orthologTools.R -- tools for handling orthologging between species


`ortholog` <- function( genes, from="Pf3D7", to="Py17X") {

	noOrtho <- ""
	gOut <- rep.int( noOrtho, length(genes))

	# do we need to reload ?
	if ( ! exists( "OrthoTable", envir=OrthoEnv)) {
		orthoTable <- NULL
		toLoad <- "OrthologTable"
		data( list=list( toLoad), envir=environment())
		if ( is.null( orthoTable)) {
			cat( "\nFailed to load Ortholog Table:  ", toLoad)
			OrthoEnv[[ "OrthoTable"]] <- NULL
			return( gOut)
		}
		
		# prep it a bit...
		OrthoEnv[[ "OrthoTable"]] <- orthoTable
		OrthoEnv[[ "OrthoStrains"]] <- colnames(orthoTable)
	}

	orthoTable <- OrthoEnv[[ "OrthoTable"]]
	orthoStrains <- OrthoEnv[[ "OrthoStrains"]]
	ifrom <- match( from, orthoStrains, nomatch=0)
	if ( ifrom == 0) {
		cat( "\n'from' strain not in Ortholog table: ", from)
		cat( "\nKnown strains: ", orthoStrains)
		return( gOut)
	}
	ito <- match( to, orthoStrains, nomatch=0)
	if ( ito == 0) {
		cat( "\n'to' strain not in Ortholog table: ", to)
		cat( "\nKnown strains: ", orthoStrains)
		return( gOut)
	}
	
	genes <- gsub( " ", "", genes)
	where <- match( genes, orthoTable[ , ifrom], nomatch=0)
	gOut[ where > 0] <- orthoTable[ where, ito]

	# the ortholog table can see the same GeneID in multiple rows, and there is no guarantee that
	# the first occurance is the 'best' one...  try to resolve...
	isOrthGenes <- which( orthoTable[ , ifrom] %in% genes)
	orthGenes <- orthoTable[ isOrthGenes, ifrom]
	isDupliGene <- which( duplicated( orthGenes))
	dupGenes <- orthGenes[ isDupliGene]
	for ( dg in dupGenes) {
		myRows <- which( orthoTable[ , ifrom] == dg)
		altGene <- as.GeneID( dg, speciesID=to)
		# if any of the choices happen to be the exact match, use that instead
		myChoices <- orthoTable[ myRows, ito]
		if ( any( myChoices == altGene)) {
			where <- match( dg, genes)
			gOut[ where] <- altGene
		}
	}

	# try to extend our set for the mammalian systems that use conserved naming
	if ( all( c( from, to) %in% MAMMAL_SPECIES)) {
		isBlank <- which( gOut == noOrtho)
		if ( length( isBlank)) {

			# we have some that still need an ortholog
			# apply the most typical name notation differences
			alternateOrthoGenes <- genes[ isBlank]
			alternateOrthoGenes <- as.GeneID( alternateOrthoGenes, speciesID=to)

			# grab the set of known genes in the 'to' species
			curSpecies <- getCurrentSpecies()
			setCurrentSpecies( to)
			gmap <- getCurrentGeneMap()
			knownGenes <- gmap$NAME
			setCurrentSpecies( curSpecies)

			# see which of those we have
			validGenes <- which( alternateOrthoGenes %in% knownGenes)
			if ( length( validGenes)) {
				gOut[ isBlank[ validGenes]] <- alternateOrthoGenes[ validGenes]
			}
		}
	}

	return( gOut)
}


`orthologSet` <- function( genes, from="Pf3D7", to="Py17X") {

	noOrtho <- ""
	gOut <- noOrtho

	# rather than at most one ortholog per gene, find all that match any in the given set

	# do we need to reload ?
	if ( ! exists( "OrthoTable", envir=OrthoEnv)) {
		orthoTable <- NULL
		toLoad <- "OrthologTable"
		data( list=list( toLoad), envir=environment())
		if ( is.null( orthoTable)) {
			cat( "\nFailed to load Ortholog Table:  ", toLoad)
			OrthoEnv[[ "OrthoTable"]] <- NULL
			return( gOut)
		}
		
		# prep it a bit...
		OrthoEnv[[ "OrthoTable"]] <- orthoTable
		OrthoEnv[[ "OrthoStrains"]] <- colnames(orthoTable)
	}

	orthoTable <- OrthoEnv[[ "OrthoTable"]]
	orthoStrains <- OrthoEnv[[ "OrthoStrains"]]
	ifrom <- match( from, orthoStrains, nomatch=0)
	if ( ifrom == 0) {
		cat( "\n'from' strain not in Ortholog table: ", from)
		cat( "\nKnown strains: ", orthoStrains)
		return( gOut)
	}
	ito <- match( to, orthoStrains, nomatch=0)
	if ( ito == 0) {
		cat( "\n'to' strain not in Ortholog table: ", to)
		cat( "\nKnown strains: ", orthoStrains)
		return( gOut)
	}

	genes <- gsub( " ", "", genes)

	# find all rows with a hit
	hits <- which( orthoTable[ , ifrom] %in% genes)
	if( ! length( hits)) return( gOut)

	gOut <- orthoTable[ hits, ito]

	# try to extend our set for the mammalian systems that use conserved naming
	if ( all( c( from, to) %in% MAMMAL_SPECIES)) {
		genesFound <- orthoTable[ hits, ifrom]
		genesNotFound <- setdiff( genes, genesFound)
		if ( length( genesNotFound)) {

			# we have some that still need an ortholog
			# apply the most typical name notation differences
			alternateOrthoGenes <- genesNotFound
			if ( from %in% c( "Mmu_grc", "Rnor_grc")) {
				alternateOrthoGenes <- toupper( alternateOrthoGenes)
			}
			if ( to %in% c( "Mmu_grc", "Rnor_grc")) {
				alternateOrthoGenes <- tolower( alternateOrthoGenes)
				substr( alternateOrthoGenes, 1,1) <- toupper( substr( alternateOrthoGenes,1,1))
			}

			# grab the set of known genes in the 'to' species
			curSpecies <- getCurrentSpecies()
			setCurrentSpecies( to)
			gmap <- getCurrentGeneMap()
			knownGenes <- gmap$NAME
			setCurrentSpecies( curSpecies)

			# see which of those we have
			validGenes <- which( alternateOrthoGenes %in% knownGenes)
			if ( length( validGenes)) {
				gOut <- c( gOut, alternateOrthoGenes[ validGenes])
			}
		}
	}

	# clean redundants and garbage
	gOut <- sort( unique( gOut))
	gOut <- setdiff( gOut, c( "", " "))
	return( gOut)
}


`allOrthologs` <- function( genes) {

	noOrtho <- ""
	N <- length(genes)
	gOut <- vector( mode="list", length=N)
	for (i in 1:N) gOut[[i]] <- noOrtho
	names(gOut) <- genes

	# do we need to reload ?
	if ( ! exists( "OrthoTable", envir=OrthoEnv)) {
		orthoTable <- NULL
		toLoad <- "OrthologTable"
		data( list=list( toLoad), envir=environment())
		if ( is.null( orthoTable)) {
			cat( "\nFailed to load Ortholog Table:  ", toLoad)
			OrthoEnv[[ "OrthoTable"]] <- NULL
			return( gOut)
		}
		
		# prep it a bit...  the table may have a 'Product' column we will ignore
		OrthoEnv[[ "OrthoTable"]] <- orthoTable
		OrthoEnv[[ "OrthoStrains"]] <- setdiff( colnames(orthoTable), "Product")
	}

	orthoTable <- OrthoEnv[[ "OrthoTable"]]
	orthoStrains <- OrthoEnv[[ "OrthoStrains"]]

	# to find all orthologs, we must check all species
	genes <- gsub( " ", "", genes)
	columnsToVisit <- setdiff( 1:ncol(orthoTable), which( colnames(orthoTable) %in% c( "Product")))
	for ( i in 1:N) {
		gset <- vector()
		for ( j in columnsToVisit) {
			wh <- which( orthoTable[ ,j] == genes[i])
			if ( length(wh)) gset <- c( gset, unlist(orthoTable[ wh, columnsToVisit]))
		}
		# strip out the gene itself and any blank entries
		if ( length(gset)) gOut[[i]] <- setdiff( sort( unique( gset)), c( genes[i], noOrtho))
	}
	return( gOut)
}


`symmetricOrthologPairs` <- function( species1, species2, verbose=T) {

	allSpecies <- getAllSpecies()
	if ( ! all( c( species1, species2) %in% allSpecies)) {
		cat( "\nGiven at least one unknown SpeciesIDs:  ", species1, species2)
		cat( "\nKnown choices:  ", allSpecies)
		return( NULL)
	}

	curSpecies <- getCurrentSpecies()

	setCurrentSpecies( species1)
	genes1 <- subset( getCurrentGeneMap(), REAL_G == TRUE)$GENE_ID
	genes1 <- shortGeneName( genes1, keep=1)
	if (verbose) cat( "\nN_Genes in ", species1, " = ", length(genes1))

	setCurrentSpecies( species2)
	genes2 <- subset( getCurrentGeneMap(), REAL_G == TRUE)$GENE_ID
	genes2 <- shortGeneName( genes2, keep=1)
	if (verbose) cat( "\nN_Genes in ", species2, " = ", length(genes2))

	setCurrentSpecies( curSpecies)

	ortho1 <- ortholog( genes1, from=species1, to=species2)
	where1in2 <- match( ortho1, genes2, nomatch=0)
	ortho2 <- ortholog( genes2, from=species2, to=species1)
	where2in1 <- match( ortho2, genes1, nomatch=0)

	# any with no match in the other species is out
	drop1 <- which( ortho1 == "" | where1in2 == 0)
	drop2 <- which( ortho2 == "" | where2in1 == 0)
	if ( length(drop1)) {
		genes1 <- genes1[ -drop1]
		ortho1 <- ortho1[ -drop1]
		if (verbose) cat( "\nN_Drop from ", species1, "for no ortholog in ", species2, " = ", length(drop1))
	}
	if ( length(drop2)) {
		genes2 <- genes2[ -drop2]
		ortho2 <- ortho2[ -drop2]
		if (verbose) cat( "\nN_Drop from ", species2, "for no ortholog in ", species1, " = ", length(drop2))
	}

	# to be a valid pair, the orthologs must be symmetric
	repeat {
		where1in2 <- match( ortho1, genes2, nomatch=0)
		where2in1 <- match( ortho2, genes1, nomatch=0)
		if ( all( c( where1in2, where2in1) > 0)) break
		drop1 <- which( where1in2 == 0)
		drop2 <- which( where2in1 == 0)
		if ( length(drop1)) {
			genes1 <- genes1[ -drop1]
			ortho1 <- ortho1[ -drop1]
			if (verbose) cat( "\nN_Drop from ", species1, "for no ortholog in ", 
					species2, " = ", length(drop1))
		}
		if ( length(drop2)) {
			genes2 <- genes2[ -drop2]
			ortho2 <- ortho2[ -drop2]
			if (verbose) cat( "\nN_Drop from ", species2, "for no ortholog in ", 
					species1, " = ", length(drop2))
		}
	}
	keeper1 <- which( where2in1[ where1in2] == 1:length( where1in2))
	keeper2 <- which( where1in2[ where2in1] == 1:length( where2in1))
	genes1 <- genes1[ keeper1]
	ortho1 <- ortho1[ keeper1]
	genes2 <- genes2[ keeper2]
	ortho2 <- ortho2[ keeper2]

	# the two answers should be identical, so merge just to be sure...
	out1 <- data.frame( genes1, ortho1, stringsAsFactors=F)
	colnames(out1) <- c( species1, species2)
	out2 <- data.frame( ortho2, genes2, stringsAsFactors=F)
	colnames(out2) <- c( species1, species2)
	out <- rbind( out1, out2)
	keys <- paste( out[[1]], out[[2]], sep=":")
	dups <- which( duplicated( keys))
	if ( length( dups)) out <- out[ -dups, ]
	rownames(out) <- 1:nrow(out)
	if (verbose) cat( "\nN_Symmetric Ortholog Genes: ", nrow(out),"\n")
	return(out)
}


`orthologTransform` <- function( x, from, to="Pf3D7", geneColumn="GENE_ID", 
				intensityColumn="RPKM_M", productColumn="PRODUCT", 
				symmetricOnly=FALSE, sep="\t", keepGenes=NULL) {

	# get what we need to know about both species
	curSpecies <- getCurrentSpecies()
	on.exit( setCurrentSpecies(curSpecies))
	fromSpecies <- from
	setCurrentSpecies(fromSpecies)
	fromPrefix <- getCurrentSpeciesFilePrefix()
	toSpecies <- to
	setCurrentSpecies(toSpecies)
	toPrefix <- getCurrentSpeciesFilePrefix()

	# we may have been given a data frame or a file name
	if ( is.data.frame(x)) {
		tbl <- x
		isFILE <- FALSE
	} else { 
		if ( length(x) > 1) cat( "\nWarning:  only using first given filename.")
		infile <- x[1]
		tbl <- read.delim( infile, as.is=T, sep=sep)
		isFILE <- TRUE
		outpath <- dirname( infile)
		fromPattern <- paste( ".", fromPrefix, ".", sep="")
		toPattern <- paste( ".", toPrefix, ".", sep="")
		outfile <- gsub( fromPattern, toPattern, basename( infile), fixed=T)
		outfile <- file.path( outpath, outfile)
	}

	gidPtr <- match( geneColumn, colnames(tbl), nomatch=0)
	if ( gidPtr == 0) {
		cat( "\nRequired GeneID column not found:  ", geneColumn)
		cat( "\nFound: ", colnames(tbl))
		return(NULL)
	}
	prodPtr <- match( productColumn, colnames(tbl), nomatch=0)

	# get the orthologs for these genes
	fromGenes <- tbl[[gidPtr]]
	shortGenes <- shortGeneName( fromGenes, keep=1)

        # allow keeping some genes by name no matter if there is no known ortholog
        if ( ! is.null( keepGenes)) {
                keepIDs <- shortGeneName( keepGenes, keep=1)
                whoKeep <- which( shortGenes %in% keepIDs)
        } else {
                whoKeep <- vector()
        }

	if (symmetricOnly) {
		orthoAns <- symmetricOrthologPairs( species1=fromSpecies, species2=toSpecies)
		if ( is.null(orthoAns)) return(NULL)
		where <- match( shortGenes, orthoAns[ , fromSpecies], nomatch=0)
		drops <- which( where == 0)
		toGenes <- orthoAns[ where, toSpecies]
	} else {
		toGenes <- ortholog( shortGenes, from=fromSpecies, to=toSpecies)
		drops <- which( toGenes == "")
                if ( length( whoKeep)) {
			drops <- setdiff( drops, whoKeep)
                	toGenes[whoKeep] <- shortGenes[whoKeep]
		}
		toGenes <- toGenes[ -drops]
	}
	# only keep those with orthologs
	if ( length( drops)) {
		out <- tbl[ -drops, ]
		cat( "\nDrop some genes due to ortholog transform:  ",nrow(tbl), " down to ", nrow(out))
	}

	# if we needed to shorten the gene names to find orthologs, expand them back to full names
	if (any( shortGenes != fromGenes)) {
		gmap <- getCurrentGeneMap()
		whereShort <- match( toGenes, shortGeneName( gmap$GENE_ID, keep=1), nomatch=0)
		toGenes[ whereShort > 0] <- gmap$GENE_ID[ whereShort]
	}

	# put these orthologs in as the true gene IDs
	out$ORIG_ID <- out[[gidPtr]]
	out[[gidPtr]] <- toGenes
	rownames(out) <- 1:nrow(out)

	# grab the species products if we can
	if (prodPtr > 0) {
		toProd <- gene2Product( toGenes)
		hasProd <- which( toProd != "")
		out[[prodPtr]][ hasProd] <- toProd[ hasProd]
	}

	# the table may have duplicate gene IDs after the ortholog conversion
	# catch and combine
	out <- combineDuplicateGenes( out, geneColumn=geneColumn, intensityColumn=intensityColumn, 
					notesColumn="ORIG_ID")

	if ( isFILE) {
		write.table( out, file=outfile, sep=sep, quote=(sep == ","), row.names=F)
		cat( "\nWrote ortholog transformed file:  ", outfile)
		return( invisible(out))
	} else {
		return( out)
	}
}

