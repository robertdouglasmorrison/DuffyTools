# as.GeneID.R -- standardize the upper/lower case and capitalization of GeneIDs

`as.GeneID` <- function( genes, speciesID=getCurrentSpecies()) {

	# default is to return things unchanged
	out <- genes
	NG <- length(genes)

	#some expect all upper case
	allCapitalSpecies <- c( PARASITE_SPECIES, "Hs_grc", "MacMu", "MacFas", "Anan", "Msmeg_mc2_155",
				"Mchel", "Mabsc")
	firstLetterCapitalSpecies <- c( "Mmu_grc", "Rnor_grc", "MT_H37")
	allLowerSpecies <- c( "Drerio")

	# first do the general case
	if ( speciesID %in% allCapitalSpecies) {
		out <- toupper( genes)
	}
	if ( speciesID %in% allLowerSpecies) {
		out <- tolower( genes)
	}
	if ( speciesID %in% firstLetterCapitalSpecies) {
		nc <- nchar( genes)
		first <- toupper( substr( genes, 1, 1))
		rest <- tolower( substr( genes, rep.int(2,NG), nc))
		out <- paste( first, rest, sep="")
	}

	# then some may have special needs
	lowerCase.C.suffix <- c( "Mabsc")
	if ( speciesID %in% lowerCase.C.suffix) {
		out <- sub( "C$", "c", out)
	}

	if ( speciesID == "MT_H37") {
		out <- sub( "a$", "A", out)
		out <- sub( "b$", "B", out)
		out <- sub( "d$", "D", out)
	}

	return( out)
}
