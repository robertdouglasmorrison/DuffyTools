# speciesTools.R


`setCurrentSpecies` <-
function( speciesID=NULL, seqID=NULL) {

	if ( all( is.null( c( speciesID, seqID)))) stop( "setCurrentSpecies:  one of 'speciesID', 'seqID' required")

	# get what maps we have stored
	speciesSet <- MapSetEnv[[ "SpeciesSet"]]
	prefixSet <- MapSetEnv[[ "PrefixSet"]]

	# they gave us a seqID
	if ( ! is.null( seqID)) {
		# look thru all mapSets till we find this SeqID
		imap <- 0
		for( i in 1:length( prefixSet)) {
			nam <- base::paste( prefixSet[i], "MapSet", sep="")
			mapset <- MapSetEnv[[ nam ]]
			smap <- mapset$seqMap
			if ( seqID %in% smap$SEQ_ID) {
				imap <- i
				break
			}
		}
		if ( imap == 0) stop( paste( "setCurrentSpecies:  no species contains SEQ_ID: ", seqID))
		loadMapSet( pos=imap)
	}

	# they gave us a speciesID
	if ( ! is.null( speciesID)) {
		# look thru all mapSets till we find this SeqID
		imap <- base::match( speciesID, speciesSet, nomatch=0)
		if ( imap == 0) stop( paste( "setCurrentSpecies:  unknown speciesID: ", speciesID))
		loadMapSet( pos=imap)
	}
	return( getCurrentSpecies() )
}


`getCurrentSpecies` <-
function() { return( MapSetEnv[[ "currentSpecies"]] ) }


`getAllSpecies` <-
function() { return( MapSetEnv[[ "SpeciesSet"]] ) }


`getAllSpeciesFilePrefixes` <-
function() { out <- MapSetEnv[[ "PrefixSet"]];  names( out) <- MapSetEnv[[ "SpeciesSet"]];  return( out) }


`getCurrentSpeciesFilePrefix` <-
function() { 
	curspecies <- MapSetEnv[[ "currentSpecies"]]
	who <- base::match( curspecies, MapSetEnv[[ "SpeciesSet"]], nomatch=0)
	return( if( who == 0) NA else MapSetEnv[[ "PrefixSet"]][ who] )
}

`getOtherSpeciesFilePrefix` <- function( speciesID) { 

	allPrefixes <- getAllSpeciesFilePrefixes()
	who <- base::match( speciesID, names( allPrefixes), nomatch=0)
	if ( any( who == 0)) stop( paste( "SpeciesID not from a known species: ", 
					base::sort( unique.default(speciesID))))
	return( allPrefixes[ who])
}


# find the species that owns this seqID
`getSpeciesFromSeqID` <- function( seqID=NULL, verbose=TRUE) {

	if ( is.null( seqID)) stop( "getSpeciesFromSeqID:  required 'seqID' argument is missing")

	# get what maps we have stored
	speciesSet <- MapSetEnv[[ "SpeciesSet"]]
	prefixSet <- MapSetEnv[[ "PrefixSet"]]

	# they gave us a seqID
	# look thru all mapSets till we find this SeqID
	out <- rep( NA, times=length( seqID))

	for( i in 1:length( prefixSet)) {
		nam <- base::paste( prefixSet[i], "MapSet", sep="")
		mapset <- MapSetEnv[[ nam ]]
		smap <- mapset$seqMap
		hits <- base::match( seqID, smap$SEQ_ID, nomatch=0)
		out[ hits > 0] <- speciesSet[ i]
	}
	if ( any( is.na( out)) && verbose) {
		cat( paste( "getSpeciesFromSeqID:  no species contains SEQ_ID: ", 
			paste( unique.default( seqID[ is.na(out)]), collapse=" | ")))
	}
	return( out)
}


`getSpeciesText` <- function( speciesID=NULL) {

	if ( is.null( speciesID)) stop( "getSpeciesText:  required 'speciesID' argument is missing")

	# get what maps we have stored
	speciesSet <- MapSetEnv[[ "SpeciesSet"]]
	prefixSet <- MapSetEnv[[ "PrefixSet"]]

	who <- match( speciesID, speciesSet, nomatch=0)
	if ( who == 0) {
		cat( "SpeciesID not found: ", speciesID, "\nKnown: ", speciesSet)
		stop()
	}

	nam <- base::paste( prefixSet[who], "MapSet", sep="")
	mapset <- MapSetEnv[[ nam ]]
	txt <- mapset$speciesText
	return( txt)
}


`test.speciesTools` <- function() {

	setCurrentSpecies( speciesID="Hs_grc")
	checkEquals( getCurrentSpecies(), "Hs_grc")
	setCurrentSpecies( seqID="Pf3D7_14_v3")
	checkEquals( getCurrentSpecies(), "Pf3D7")
	checkEquals( getCurrentSpeciesFilePrefix(), "Pf")
}

