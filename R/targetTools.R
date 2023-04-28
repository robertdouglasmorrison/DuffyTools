# targetTools.R

# utilities to manage the Target Organisms and Indexes inside the DuffyRNAseq package


# quick wrappers to get one from the targets environment
`getCurrentTarget` <- function() { 
	
	ans <- TargetEnv[[ "CurrentTarget" ]]
	speciesSet <- getCurrentTargetSpecies()
	txt <- sapply( speciesSet, FUN=getSpeciesText)
	txt <- paste( txt, collapse="; ")
	ans$SpeciesText <- txt
	return( ans)
}


`setCurrentTarget` <- function( targetID=NULL, optionsFile=NULL) {

	allTargets <- TargetEnv[[ "AllTargets"]]

	if ( is.null( targetID)) {
		targetID <- getOptionValue( optionsFile, "targetID", notfound="HsPf")
	}
	where <- base::match( targetID, allTargets$TargetID, nomatch=0)
	if ( where == 0) {
		stop( paste( "setCurrentTarget:  Error:  targetID not found.  Given: ", 
				targetID, "Current choices:  ", 
				paste( allTargets$TargetID, collapse=" | ")))
	}
	
	assign( "CurrentTarget", value=allTargets[ where, ],  envir=TargetEnv)

	# force the current species to the first in this target
	speciesSet <- getCurrentTargetSpecies()
	setCurrentSpecies( speciesID=speciesSet[1])

	return( allTargets$TargetID[ where])
}


`getAndSetTarget` <- function( optionsFile="Options.txt", sampleID=NULL, annotationFile="Annotation.txt", verbose=F) {

	target <- getOptionValue( optionsFile, "targetID", notfound="HsPf", verbose=verbose)

	# let the annotation file override...
	if ( ! is.null( sampleID)) {
		sampleID <- originalSamplePairID( sampleID[1], annotationFile)
		target <- getAnnotationValue( annotationFile, key=sampleID, columnArg="TargetID", 
						notfound=target, verbose=verbose)
	}

	setCurrentTarget( target)
	if (verbose) cat( "\nTarget Species Set:    ", target, "\t", getCurrentTargetSpecies(), "\n")
	return( target)
}


# get the entire Target Set at once
`getAllTargets` <- function() { 

	ans <- TargetEnv[[ "AllTargets"]]
	ans$SpeciesText <- ""
	for ( i in 1:nrow(ans)) {
		speciesSet <- strsplit( ans$SpeciesSet[i], split=",")[[1]]
		txt <- sapply( speciesSet, FUN=getSpeciesText)
		txt <- paste( txt, collapse="; ")
		ans$SpeciesText[i] <- txt
	}
	return(ans)
}


# get the set of speciesIDs for the current target
`getCurrentTargetSpecies` <- function() { 

	thisTarget <- TargetEnv[[ "CurrentTarget" ]]
	if (nrow( thisTarget) < 1) return( vector())
	speciesSet <- strsplit( gsub( " ", "", thisTarget$SpeciesSet), split=",", fixed=T)[[1]]
	return( speciesSet)
}


# get the set of speciesIDs for the current target
`getCurrentTargetFilePrefix` <- function() { 

	thisTarget <- TargetEnv[[ "CurrentTarget" ]]
	if (nrow( thisTarget) < 1) return( vector())
	prefixSet <- strsplit( gsub( " ", "", thisTarget$PrefixSet), split=",", fixed=T)[[1]]
	return( prefixSet)
}


# add an explicit target to the targets envirn
`addTarget` <- function( targetID, speciesSet=targetID, prefixSet=speciesSet, mapset.path=NULL) {

	speciesVec <- strsplit( gsub( " ", "", speciesSet), split=",", fixed=T)[[1]]
	prefixVec <- strsplit( gsub( " ", "", prefixSet), split=",", fixed=T)[[1]]
	if ( length(speciesVec) != length( prefixVec)) {
		cat( "\nNeed the same number of species IDs and file prefixes")
		cat( "\nSpecies:  ", speciesVec)
		cat( "\nPrefixes: ", prefixVec)
		stop()
	}

	# make a one line data.frame for this given target
	thisTarget <- data.frame( targetID, speciesSet, prefixSet, stringsAsFactors=FALSE)
	colnames( thisTarget) <- TARGET_NAMES

	# get the currently known targets and see if its already in
	allTargets <- TargetEnv[[ "AllTargets"]]
	who <- (-1)
	if (nrow( allTargets) > 0) {
		who <- base::match( targetID, allTargets$TargetID, nomatch=0)
	}

	# either overwrite, append, or start fresh
	if ( who > 0) {
		allTargets[ who, ] <- thisTarget[ 1, ]
		myRow <- who
	} else if ( who == 0) {
		allTargets <- rbind( allTargets, thisTarget)
		myRow <- nrow( allTargets)
	} else {
		allTargets <- thisTarget
		myRow <- 1
	}
	rownames( allTargets) <- 1:nrow(allTargets)

	# store the complete set, and the current one
	assign( "AllTargets", value=allTargets, envir=TargetEnv)
	assign( "CurrentTarget", value=allTargets[myRow, ], envir=TargetEnv)

	# force the load of these mapsets
	for ( i in 1:length(speciesVec)) {
		spec <- speciesVec[i]
		if ( spec %in% getAllSpecies()) next
		cat( "  ", spec, "..", sep="")
		mapsetName <- paste( prefixVec[i], "MapSet", sep=".")
		mapSet <- findAndLoadMapSet( mapsetName, mapset.path=mapset.path)
		if ( is.null( mapSet)) {
			cat( "\nFailed to find/load MapSet data...  Species: ", spec, "\tPrefix: ", prefixVec[i], "\n")
		}
	}

	setCurrentSpecies( speciesVec[1])
	return( targetID)
}


# reset to factory defaults
`target.defaults` <- function() {

	# load the package mapSets, and stuff them into the mapSet environment
	assign( "AllTargets", value=data.frame(), envir=TargetEnv)
	assign( "CurrentTarget", value=data.frame(), envir=TargetEnv)

	cat( "\nLoading Target Species: ")
	addTarget( "Pf3D7", "Pf3D7", "Pf")
	addTarget( "Hs_grc", "Hs_grc", "Hs")
	addTarget( "HsPf", "Pf3D7,Hs_grc", "Pf,Hs")

	addTarget( "PCO", "PCO", "Pco")
	addTarget( "PcoAH", "PcoAH", "PcoAH")
	addTarget( "PkH", "PkH", "PkH")
	addTarget( "MacMu", "MacMu", "MacMu")
	addTarget( "MacFas", "MacFas", "MacFas")
	addTarget( "Anan", "Anan", "Anan")
	addTarget( "MacMuPCO", "PCO,MacMu", "Pco,MacMu")
	addTarget( "MacMuPcoAH", "PcoAH,MacMu", "PcoAH,MacMu")
	addTarget( "MacMuPkH", "PkH,MacMu", "PkH,MacMu")
	addTarget( "MacFasPCO", "PCO,MacFas", "Pco,MacFas")
	addTarget( "MacFasPcoAH", "PcoAH,MacFas", "PcoAH,MacFas")
	addTarget( "MacFasPkH", "PkH,MacFas", "PkH,MacFas")
	addTarget( "AnanPcoAH", "PcoAH,Anan", "PcoAH,Anan")
	addTarget( "AnanPf", "Pf3D7,Anan", "Pf,Anan")

	addTarget( "PbANKA", "PbANKA", "Pb")
	addTarget( "PbMmu", "PbANKA,Mmu_grc", "Pb,Mmus")
	addTarget( "PbHs", "PbANKA,Hs_grc", "Pb,Hs")

	addTarget( "Py17X", "Py17X", "Py17X")
	addTarget( "Mmu_grc", "Mmu_grc", "Mmus")
	addTarget( "PyMmu", "Py17X,Mmu_grc", "Py17X,Mmus")
	addTarget( "PyHsMmu", "Py17X,Hs_grc,Mmu_grc", "Py17X,Hs,Mmus")

	addTarget( "Ocun", "Ocun", "Ocun")
	addTarget( "Rnor", "Rnor", "Rnor")

	addTarget( "PfHsMmu", "Pf3D7,Hs_grc,Mmu_grc", "Pf,Hs,Mmus")

	addTarget( "PyYM", "PyYM", "Pym")
	addTarget( "PymMmu", "PyYM,Mmu_grc", "Pym,Mmus")
	addTarget( "PymHsMmu", "PyYM,Hs_grc,Mmu_grc", "Pym,Hs,Mmus")

	addTarget( "PvSal1", "PvSal1", "Pv")
	addTarget( "HsPv", "Hs_grc,PvSal1", "Hs,Pv")
	addTarget( "PvHsMmu", "PvSal1,Hs_grc,Mmu_grc", "Pv,Hs,Mmus")
	addTarget( "PvP01", "PvP01", "PvP01")
	addTarget( "HsPvP01", "Hs_grc,PvP01", "Hs,PvP01")

	addTarget( "PchAS", "PchAS", "Pch")
	addTarget( "PchMmu", "PchAS,Mmu_grc", "Pch,Mmus")

	addTarget( "Pcy", "Pcy", "Pcy")
	addTarget( "MacFasPcy", "Pcy,MacFas", "Pcy,MacFas")
	addTarget( "MacFasPkH", "PkH,MacFas", "PkH,MacFas")

	addTarget( "Pvvv", "Pvvv", "Pvvv")
	addTarget( "PvvvMmu", "Pvvv,Mmu_grc", "Pvvv,Mmus")

	addTarget( "Povale", "Povale", "Povale")
	addTarget( "Pmalar", "Pmalar", "Pmalar")

	addTarget( "Agam", "Agam", "Ag")
	addTarget( "AgPf", "Pf3D7,Agam", "Pf,Ag")
	addTarget( "AgPv", "PvSal1,Agam", "Pv,Ag")
	addTarget( "AgPy", "Py17X,Agam", "Py17X,Ag")
	addTarget( "Asteph", "Asteph", "Asteph")
	addTarget( "AstephPf", "Pf3D7,Asteph", "Pf,Asteph")

	addTarget( "Ecoli", "Ecoli", "Eco")

	addTarget( "Dmel", "Dmel", "Dmel")
	# for Sanaria SPX, a blend of PF + A.steph + Drosophila
	addTarget( "AstephDmelPf", "Pf3D7,Dmel,Asteph", "Pf,Dmel,Asteph")

	addTarget( "MT_H37", "MT_H37", "MTb")
	addTarget( "MT_HN878", "MT_HN878", "MTbHN878")
	addTarget( "HsMTb", "MT_H37,Hs_grc", "MTb,Hs")
	addTarget( "HsMTbHN878", "MT_HN878,Hs_grc", "MTbHN878,Hs")
	addTarget( "MTbMmu", "MT_H37,Mmu_grc", "MTb,Mmus")
	addTarget( "MTbHN878Mmu", "MT_HN878,Mmu_grc", "MTbHN878,Mmus")
	addTarget( "MTbOcun", "MT_H37,Ocun", "MTb,Ocun")

	addTarget( "Styph", "Styph_sl1344", "Styph")
	addTarget( "StyphHs", "Styph_sl1344,Hs_grc", "Styph,Hs")
	addTarget( "StyphMmu", "Styph_sl1344,Mmu_grc", "Styph,Mmus")

	addTarget( "Msmeg_mc2_155", "Msmeg_mc2_155", "Msmeg")
	addTarget( "Mabsc", "Mabsc", "Mabsc")
	addTarget( "Mchel", "Mchel", "Mchel")
	addTarget( "Mavium", "Mavium", "Mavium")
	addTarget( "Mmar", "Mmar", "Mmar")
	addTarget( "Bbrevis", "Bbrevis", "Bbrevis")

	addTarget( "Tb927", "Tb927", "Tbr")
	addTarget( "Tb427", "Tb427", "Tbr427")
	
	addTarget( "PfDD2", "PfDD2", "PfDD2")
	addTarget( "PfIT", "PfIT", "PfIT")

	addTarget( "VarGenes", "VarGenes", "VarGenes")
	addTarget( "PF.VSA.JOS", "Pf3D7,VSA,JOS", "Pf,VSA,JOS")
	addTarget( "VSA.JOS", "VSA,JOS", "VSA,JOS")

	#addTarget( "GrammDol", "GrammDol", "GrammDol")
	addTarget( "Gsurd", "Gsurd", "Gsurd")

	addTarget( "Drerio", "Drerio", "Dr")

	addTarget( "KSHV", "KSHV", "KSHV")
	addTarget( "HsKSHV", "Hs_grc,KSHV", "Hs,KSHV")

	addTarget( "CMV", "CMV", "CMV")
	addTarget( "BCG", "BCG", "BCG")

	# Pichia for cell production systems
	addTarget( "Kphaf", "Kphaf", "Kphaf")

	setCurrentTarget( "HsPf")
	cat( "  Done.\n")
	return()
}


`exportTargets` <- function( fileout="DuffyTools.Targets.txt") {

	tmp <- getAllTargets()
	write.table( tmp, file=fileout, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	cat( "\nWrote targets file:  ", fileout, "\nN_Targets: ", nrow( tmp), "\n")
	return()
}


`importTargets` <- function( filein="DuffyTools.Targets.txt") {


	if ( ! file.exists(filein)) stop( paste( "importTargets:  cannot file targets file",
			"\nTried: ", filein))

	tmp <- read.delim( filein, header=TRUE, as.is=TRUE)

	# re-initialize to empty
	assign( "AllTargets", value=data.frame(), envir=TargetEnv)
	assign( "CurrentTarget", value=data.frame(), envir=TargetEnv)

	for ( i in 1:nrow(tmp)) {
		addTarget( tmp$TargetID[i], tmp$SpeciesSet[i], tmp$PrefixSet[i])
	}

	# use the first as default
	assign( "CurrentTarget", value=tmp$TargetID[1], envir=TargetEnv)
	# force the current species to the first in this target
	speciesSet <- getCurrentTargetSpecies()
	setCurrentSpecies( speciesID=speciesSet[1])

	return( getAllTargets()$TargetID)
}

