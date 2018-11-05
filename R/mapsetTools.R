# mapsetTools.R

# utilities to manage the map sets inside the DuffyRNAseq package


# quick wrappers to get one map
`getCurrentGeneMap` <- function() { return( MapSetEnv[[ "geneMap" ]] ) }
`getCurrentExonMap` <- function() { return( MapSetEnv[[ "exonMap" ]] ) }
`getCurrentCdsMap` <- function() { 
				if ( exists( "cdsMap", envir=MapSetEnv)) {
					return( MapSetEnv[[ "cdsMap" ]] ) 
				} else {
					return( MapSetEnv[[ "exonMap" ]] )
				}
			}
`getCurrentRrnaMap` <- function() { return( MapSetEnv[[ "rrnaMap" ]] ) }
`getCurrentSeqMap` <- function() { return( MapSetEnv[[ "seqMap" ]] ) }


# get the entire mapset at once
`getCurrentMapSet` <- function() { 
	prefix <- getCurrentSpeciesFilePrefix()
	storageName <- base::paste( prefix, "MapSet", sep="")
	obj <- MapSetEnv[[ storageName]]
	return( obj)
}


# write the current mapset to flat text files...
`exportCurrentMapSet` <- function( path=".") {

	species <- getCurrentSpecies()
	if ( is.null( species)) {
		warning( "No Current Species is set...")
		return()
	}
	prefix <- getCurrentSpeciesFilePrefix()
	mapset <- getCurrentMapSet()

	if ( ! file.exists( path)) dir.create( path, recursive=TRUE)

	cat( "\nExporting MapSet for:  ", species, "...")
	outfile <- file.path( path, paste( prefix, "SeqMap", "txt", sep="."))
	write.table( getCurrentSeqMap(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	cat( outfile, "...")
	outfile <- file.path( path, paste( prefix, "GeneMap", "txt", sep="."))
	write.table( getCurrentGeneMap(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	cat( outfile, "...")
	outfile <- file.path( path, paste( prefix, "ExonMap", "txt", sep="."))
	write.table( getCurrentExonMap(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	cat( outfile, "...")
	outfile <- file.path( path, paste( prefix, "CdsMap", "txt", sep="."))
	write.table( getCurrentCdsMap(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	cat( outfile, "...")
	outfile <- file.path( path, paste( prefix, "RrnaMap", "txt", sep="."))
	write.table( getCurrentRrnaMap(), file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	cat( outfile, "...")
	outfile <- file.path( path, paste( prefix, "MapSetNames", "txt", sep="."))
	tmpdf <- data.frame( "speciesID"=mapset$speciesID, "speciesFilePrefix"=mapset$speciesFilePrefix, 
			"speciesText"=mapset$speciesText, stringsAsFactors=FALSE)
	write.table( tmpdf, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	cat( outfile, "\n")
}


# import the current mapset from a set of flat text files...
`importMapSet` <- function( speciesFilePrefix=NULL, path=".", checkOverlaps=FALSE) {

	prefix <- speciesFilePrefix
	if ( is.null( prefix)) {
		warning( "Explicit Species file prefix is required...")
		return()
	}

	cat( "\nImporting MapSet for Prefix:  ", prefix, "...")
	infile <- file.path( path, paste( prefix, "MapSetNames", "txt", sep="."))
	tmpdf <- read.delim( file=infile, as.is=TRUE)
	cat( "\n")
	print( tmpdf)
	infile <- file.path( path, paste( prefix, "SeqMap", "txt", sep="."))
	seqmap <- read.delim( file=infile, as.is=TRUE)
	cat( infile, "...")
	infile <- file.path( path, paste( prefix, "GeneMap", "txt", sep="."))
	genemap <- read.delim( file=infile, as.is=TRUE)
	cat( infile, "...")
	infile <- file.path( path, paste( prefix, "ExonMap", "txt", sep="."))
	exonmap <- read.delim( file=infile, as.is=TRUE)
	cat( infile, "...")
	infile <- file.path( path, paste( prefix, "CdsMap", "txt", sep="."))
	if ( ! file.exists( infile)) {
		cat( "  using Exons as CDS..")
		infile <- file.path( path, paste( prefix, "ExonMap", "txt", sep="."))
	}
	cdsmap <- read.delim( file=infile, as.is=TRUE)
	cat( infile, "...")
	infile <- file.path( path, paste( prefix, "RrnaMap", "txt", sep="."))
	rrnamap <- read.delim( file=infile, as.is=TRUE)
	cat( infile, "...")
	
	mapset <- list( "speciesID"=tmpdf$speciesID, "speciesFilePrefix"=tmpdf$speciesFilePrefix,
			"speciesText"=tmpdf$speciesText, "seqMap"=seqmap, "geneMap"=genemap, 
			"exonMap"=exonmap, "cdsMap"=cdsmap, "rrnaMap"=rrnamap)

	cat( "\nValidating...")
	validateMapSet( mapset, checkOverlaps=checkOverlaps)

	cat( "\nSaving new 'MapSet.rda' file...")
	mapSet <- mapset
	save( mapSet, file=file.path( path, paste( prefix, "MapSet.rda", sep=".")))

	cat( "\nAdding MapSet...")
	addMapSet( mapset)

	cat( "\nDone.")
	return( invisible( mapset))
}



# add a complete mapset to the MapSet environment:  user level
`addMapSet` <- function( mapset) {

	if ( ! all( MAPSET_NAMES %in% names( mapset))) stop( paste( "addMapset:  invalid mapSet names.  Expected: ", 
		paste( MAPSET_NAMES, collapse="  ")))

	# add/overwrite the given mapSet to the map environment
	species <- mapset$speciesID
	prefix <- mapset$speciesFilePrefix

	# find where this species goes
	speciesSet <- prefixSet <- vector()

	# first see if we have any species already
	if ( "SpeciesSet" %in% ls( envir=MapSetEnv)) {
		speciesSet <- MapSetEnv[[ "SpeciesSet"]]
		prefixSet <- MapSetEnv[[ "PrefixSet"]]
	}
	
	# we either overwrite or append
	if ( species %in% speciesSet) {
		loc <- base::match( species, speciesSet)
		prefixSet[loc] <- prefix
	} else {
		loc <- length( speciesSet) + 1
		speciesSet[loc] <- species
		prefixSet[loc] <- prefix
	}

	# verify that this speices and prefix are truly unique
	allHitsPrefix <- which( prefixSet == prefix)
	allHitsSpecies <- which( base::substr( speciesSet, 1, base::nchar(species)) == species)
	if ( any( allHitsPrefix != loc) || any( allHitsSpecies != loc)) {
		stop( paste( "Non-unique SpeciesID and/or FilePrefix:  ", species,  prefix, "\nCurrent:", 
			paste( speciesSet, prefixSet, sep=" - ", collapse="\t")))

	}

	# good, so store the names and that mapset
	MapSetEnv[[ "SpeciesSet" ]] <- speciesSet
	MapSetEnv[[ "PrefixSet" ]] <- prefixSet

	obj <- mapset
	storageName <- base::paste( prefix, "MapSet", sep="")
	MapSetEnv[[ storageName ]] <- obj

	return( species)
}


# load a resident mapset as the current one, by species name or by direct pointer:  low level
`loadMapSet` <- function( speciesID=NULL, pos=NULL) {

	# give a position pointer
	if ( ! is.null( pos)) {
		prefix <- MapSetEnv[[ "PrefixSet"]][pos]
		speciesID <- MapSetEnv[[ "SpeciesSet"]][pos]
		storageName <- base::paste( prefix, "MapSet", sep="")
		storedMapset <- MapSetEnv[[ storageName]]
		for ( imap in 1:length( MAPSET_NAMES)) {
			nam <- MAPSET_NAMES[ imap]
			storageFrom <- nam
			storageTo <- nam
			MapSetEnv[[ storageTo ]] <- storedMapset[[ storageFrom ]]
		}
		# currently, a 'cdsMap' is an optional member of the MapSet
		if ( "cdsMap" %in% names( storedMapset)) {
			MapSetEnv[[ "cdsMap" ]] <- storedMapset[[ "cdsMap" ]]
		} else {
			# if no CDS explicitly given, make sure to remove any from the 'current' storage!!
			if ( "cdsMap" %in% ls( MapSetEnv)) {
				rm( "cdsMap", envir=MapSetEnv)
			}
		}
		MapSetEnv[[ "currentSpecies"]] <- speciesID
		return( speciesID)
	}

	# given a speciesID
	if ( ! is.null( speciesID)) {
		pos <- base::match( speciesID, MapSetEnv[[ "SpeciesSet"]], nomatch=0)
		if ( pos == 0) stop( paste( "loadMapSet:  invalid speciesID: ", speciesID, "  Not a known speciesID"))
		prefix <- MapSetEnv[[ "PrefixSet"]][pos]
		storageName <- base::paste( prefix, "MapSet", sep="")
		storedMapset <- MapSetEnv[[ storageName]]
		for ( imap in 1:length( MAPSET_NAMES)) {
			nam <- MAPSET_NAMES[ imap]
			storageFrom <- nam
			storageTo <- nam
			MapSetEnv[[ storageTo ]] <- storedMapset[[ storageFrom ]]
		}
		# currently, a 'cdsMap' is an optional member of the MapSet
		if ( "cdsMap" %in% names( storedMapset)) {
			MapSetEnv[[ "cdsMap" ]] <- storedMapset[[ "cdsMap" ]]
		} else {
			# if no CDS explicitly given, make sure to remove any from the 'current' storage!!
			if ( "cdsMap" %in% ls( MapSetEnv)) {
				rm( "cdsMap", envir=MapSetEnv)
			}
		}
		MapSetEnv[[ "currentSpecies"]] <- speciesID
		return( speciesID)
	}
	stop( "loadMapset:  required argument 'pos' or 'speciesID' not given")
}
	

# reset to factory defaults
`mapset.defaults` <- function() {

	# load the package mapSets, and stuff them into the mapSet environment
	cat( "\nLoading MapSets:    \tP.falciparum..")
	findAndLoadMapSet( "Pf.MapSet")

	cat( " P.knowlesi..")
	findAndLoadMapSet( "Pk.MapSet")

	cat( " P.berghei..")
	findAndLoadMapSet( "Pb.MapSet")

	cat( " P.yoelii..")
	findAndLoadMapSet( "Py.MapSet")

	cat( " H.sapians..")
	findAndLoadMapSet( "Hs.MapSet")

	cat( " M.musculus..")
	findAndLoadMapSet( "Mmus.MapSet")

	cat( " M.mulatta..")
	findAndLoadMapSet( "MacMu.MapSet")

	cat( " T.brucei..")
	findAndLoadMapSet( "Tb.MapSet")

	cat( " L.major..")
	findAndLoadMapSet( "Lmj.MapSet")

	cat( " M.tuberculosis..")
	findAndLoadMapSet( "MTb.MapSet")

	cat( " M.smegmatis..")
	findAndLoadMapSet( "Msmeg.MapSet")

	cat( " E.coli..")
	findAndLoadMapSet( "Eco.MapSet")

	cat( " S.typhimurium..")
	findAndLoadMapSet( "Sty.MapSet")

	cat( " A.gambiae..")
	findAndLoadMapSet( "Agam.MapSet")

	cat( "Done.\n")

	# start in Pf
	setCurrentSpecies( "Pf3D7")
	return()
}


`findAndLoadMapSet` <- function( mapsetName, mapset.path=NULL) {

	mapSet <- NULL
	mapsetFile <- paste( mapsetName, "rda", sep=".")

	# first look in the explicit path, CWD, and then HOME
	if ( ! is.null( mapset.path)) {
		tryFile <- file.path( mapset.path[1], mapsetFile)
		if ( file.exists( tryFile)) {
			who <- load( tryFile)
			if ( who == "mapSet") {
				cat( " ")
				addMapSet( mapset=mapSet)
				return( mapSet)
			}
		}
	}

	tryFile <- file.path( ".", mapsetFile)
	if ( file.exists( tryFile)) {
		who <- load( tryFile)
		if ( who == "mapSet") {
			cat( "(CWD)")
			addMapSet( mapset=mapSet)
			return( mapSet)
		}
	}
	tryFile <- file.path( Sys.getenv( "HOME"), mapsetFile)
	if ( file.exists( tryFile)) {
		who <- load( tryFile)
		if ( who == "mapSet") {
			cat( "(HOME)")
			addMapSet( mapset=mapSet)
			return( mapSet)
		}
	}

	# lastly, try the installed DuffyTools data area
	# this dataset is called 'mapSet' and has all the annotation maps
	data( list=mapsetName,  envir=environment())
	if ( is.null( mapSet)) {
		cat( "\nFailed to load: ", mapsetFile, "  Not found in CWD, HOME, or DuffyTools library")
	} else {
		addMapSet( mapset=mapSet)
	}
	return( mapSet)
}


`dropAntiSenseGenes` <- function( geneMap) {

	dropNR <- which( geneMap$REAL_G == FALSE)
	dropAS <- grep( "_AS$", geneMap$GENE_ID)
	drops <- intersect( dropNR, dropAS)
	if ( length( drops) > 0) {
		geneMap <- geneMap[ -drops, ]
		rownames(geneMap) <- 1:nrow(geneMap)
	}
	return( geneMap)
}


`updateMapSetGeneProductTerms` <- function( geneProductFile, speciesID="Pf3D7", 
			geneColumn="GENE_ID", productColumn="PRODUCT", sep="\t") {

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	geneMap <- getCurrentGeneMap()
	exportCurrentMapSet()

	tbl <- read.delim( geneProductFile, sep=sep, as.is=T)
	tblIDs <- tbl[[geneColumn]]
	if ( is.null( tblIDs)) stop( paste( "'GeneID' column not found: ", geneColumn))
	tblProds <- tbl[[productColumn]]
	if ( is.null( tblProds)) stop( paste( "'Product' column not found: ", productColumn))

	# only change the product term that we do find...
	where <- match( geneMap$GENE_ID, tblIDs, nomatch=0)
	newProds <- oldProds <- geneMap$PRODUCT
	newProds[ where > 0] <- tblProds[ where]
	changed <- which( newProds != oldProds)
	if ( length(changed)) {
		cat( "\nUpdated Product terms: ", length( changed), "\n")
		print( paste( geneMap$GENE_ID[changed], ": ", oldProds[changed], " -> ", newProds[changed], sep=""))
	}

	geneMap$PRODUCT <- newProds
	geneFile <- paste( prefix, "GeneMap.txt", sep=".")
	write.table( geneMap, geneFile, sep="\t", quote=F, row.names=F)

	importMapSet( prefix)
}


`test.mapsetTools` <- function() {

	# load a mapSet into the local enviroment, and see if matches via MapSet environment lookup
	data( Pf.MapSet, envir=environment())
	setCurrentSpecies( "Pf3D7")
	checkEquals( mapSet$seqMap, getCurrentSeqMap())
	loadMapSet( pos=2)
	checkEquals( MapSetEnv[[ "currentSpecies" ]], "Hs_grc" )
	loadMapSet( speciesID="Pf3D7")
	checkEquals( MapSetEnv[[ "currentSpecies" ]], "Pf3D7" )
	checkEquals( mapSet$seqMap, getCurrentSeqMap())
}
