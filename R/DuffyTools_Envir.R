# DuffyTools_Envir.R

# set up and manage the 'global environment' for the DuffyTools package

# one shared environment for all DuffyTools objects
DuffyToolsEnv <- new.env( hash=TRUE, parent=emptyenv())


# global constants
FASTQ_COLUMNS <- c( "READ_ID", "READ_SEQ", "SCORE")


# mapSet environment
MapSetEnv <- new.env( hash=TRUE, parent=emptyenv())
MAPSET_NAMES <- c( "speciesID", "speciesFilePrefix", "speciesText", "seqMap", "geneMap", 
		"exonMap", "rrnaMap")


# Targets environment
TargetEnv <- new.env( hash=TRUE, parent=emptyenv())
TARGET_NAMES <- c( "TargetID", "SpeciesSet", "PrefixSet")

# group names for classes of species
MAMMAL_SPECIES <- c( "Hs_grc", "Mmu_grc", "Rnor", "MacMu", "MacFas", "Anan", "Drerio", "Gsurd", "Ocun")
MAMMAL_PREFIXES <- c( "Hs", "Mmus", "Rnor", "MacMu", "MacFas", "Anan", "Dr", "Gsurd", "Ocun")
PARASITE_SPECIES <- c( "Pf3D7", "Py17X", "PyYM", "PvSal1", "PvP01", "PbANKA", "PkH", "PchAS", "PCO", "Pcy", "Pvvv", 
			"PcoAH", "Pvvv", "Povale", "Pmalar")
PARASITE_PREFIXES <- c( "Pf", "Py17X", "PyYM", "Pv", "PvP01", "Pb", "Pk", "Pch", "Pco", "Pcy", "Pvvv", 
			"PcoAH", "Pvvv", "Povale", "Pmalar")
ORIGID_PARASITE_SPECIES <- c( "Pf3D7", "Py17X", "PyYM", "PbANKA")
ORIGID_PARASITE_PREFIXES <- c( "Pf", "Py17X", "PyYM", "Pb")
BACTERIA_SPECIES <- c( "MT_H37", "Msmeg_mc2_155", "Styph", "Mchel", "Mabsc", "Mavium", "MT_HN878", "BCG")
BACTERIA_PREFIXES <- c( "MTb", "Msmeg", "Styph", "Mchel", "Mabsc", "Mavium", "MTbHN878", "BCG")
INSECT_SPECIES <- c( "Agam", "Asteph", "Dmel")
INSECT_PREFIXES <- c( "Ag", "Asteph", "Dmel")


# Codon environment
STOP_CODON <- "*"
UNKNOWN_CODON <- "?"
STOP_CODON_PATTERN <- "\\*|\\?"
CodonEnv <- new.env( hash=TRUE, parent=emptyenv())


# Alias tools
AliasEnv <- new.env( hash=TRUE, parent=emptyenv())


# LifeCycle & CellType tools
LifeCycleEnv <- new.env( hash=TRUE, parent=emptyenv())
CellTypeEnv <- new.env( hash=TRUE, parent=emptyenv())


# Ortholog tools
OrthoEnv <- new.env( hash=TRUE, parent=emptyenv())


# GeneSet tools
ALL_GENE_SETS <- c( "Blood.GeneModules", "BloodGen3.GeneModules", "BTM.AntibodyResponses", 
			"GO.BiologicalProcess", "GO.CellularComponent", "GO.MolecularFunction", 
			"GOslim.BiologicalProcess", "GOslim.CellularComponent", "GOslim.MolecularFunction", 
			"KEGG.Pathways", "MetabolicPathways", "Panther.Pathways", "PBMC.GeneModules", 
			"HGNC.GeneFamily", "Reactome", "WikiPathways", "HumanImmuneSubsets", "MouseImmuneSubsets", 
			"STRING", "CellTypes", "MPMP.Pathways", "ParasiteLifeCycle")
MAMMAL_GENE_SETS <- c( "Blood.GeneModules", "BloodGen3.GeneModules", "BTM.AntibodyResponses", 
			"GO.BiologicalProcess", "GO.CellularComponent", "GO.MolecularFunction", 
			"GOslim.BiologicalProcess", "GOslim.CellularComponent", "GOslim.MolecularFunction", 
			"KEGG.Pathways", "MetabolicPathways", "Panther.Pathways", "PBMC.GeneModules", 
			"HGNC.GeneFamily", "Reactome", "WikiPathways", "HumanImmuneSubsets", "MouseImmuneSubsets", 
			"STRING", "CellTypes")
PARASITE_GENE_SETS <- c( "Blood.GeneModules", "BloodGen3.GeneModules", 
			"GO.BiologicalProcess", "GO.CellularComponent", "GO.MolecularFunction", 
			"GOslim.BiologicalProcess", "GOslim.CellularComponent", "GOslim.MolecularFunction", 
			"KEGG.Pathways", "MetabolicPathways", "Panther.Pathways", "PBMC.GeneModules", "STRING",
			"HGNC.GeneFamily", "Reactome", "WikiPathways", "MPMP.Pathways", "ParasiteLifeCycle")
BACTERIA_GENE_SETS <- c( "GO.BiologicalProcess", "GO.CellularComponent", "GO.MolecularFunction", 
			"GOslim.BiologicalProcess", "GOslim.CellularComponent", "GOslim.MolecularFunction", 
			"KEGG.Pathways", "MetabolicPathways", "HGNC.GeneFamily", "Panther.Pathways", "Reactome", 
			"WikiPathways", "STRING",
			"TFOE.Regulons", "Tuberculist.FunctionalGroups", "Tuberculist.GO.Ontology", "ISB.Corems",
			"iModulons","Regulons")

# TimeHour tools
TimeHourEnv <- new.env( hash=TRUE, parent=emptyenv())


# Base Depth constant
EMPTY_BASE_DEPTH_TABLE <- data.frame( "START"=vector( mode="numeric", length=0), 
					"STOP"=vector( mode="numeric", length=0), 
					"DEPTH"=vector( mode="numeric", length=0))

EMPTY_BASE_DEPTH_VECTOR <- vector( mode="numeric", length=0)


# .onLoad() function is called when the package get loaded at runtime
# any needed setup goes here...
`.onLoad` <- function( libname, pkgname) {
}

`.onUnload` <- function( libpath) {
}


# .onAttach() function is called when the package gets attached,
# i.e. at the time the user first has access to the package
`.onAttach`  <- function( libname, pkgname) {

	# wake-up message
	cat( "\nPackage: \t\t", pkgname)

	# save the library and package name...
	assign( "LibraryName", value=libname, envir=DuffyToolsEnv)
	assign( "PackageName", value=pkgname, envir=DuffyToolsEnv)

	# initialize...
	DuffyTools.defaults()
}


`DuffyTools.defaults` <- function() {

	#mapset.defaults()
	target.defaults()
	codon.defaults()
}

