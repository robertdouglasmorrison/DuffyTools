# checkDuffyPackages.R -- verify needed packages...

`checkDuffyPackages` <- function() {

	# see what is installed
	current <- installed.packages()
	packages <- current[ ,"Package"]

	Nadd <- 0
	Nfail <- 0
	addSet <- failSet <- vector()

	# some are regular CRAN packages
	CRAN_Set <- c( "RCurl", "RUnit", "methods", "plotrix", "gplots", "heatmap.plus", "ape", "GenSA")

	for ( pack in CRAN_Set) {
		if ( ! ( pack %in% packages)) {
			cat( "\nPackage '", pack, "' not installed yet.   Loading from CRAN..\n", sep="")
			success <- TRUE
			tryCatch( install.packages( pack), error=function(e) { success <<- FALSE})
			if ( success) {
				Nadd <- Nadd + 1
				addSet[Nadd] <- pack
			} else {
				Nfail <- Nfail + 1
				failSet[Nfail] <- pack
			}
		}
	}


	# some are Bioconductor...
	source( "http://bioconductor.org/biocLite.R")

	BIOC_Set <- c( "Biostrings", "DESeq", "DESeq2", "edgeR", "siggenes", "ROC", "sangerseqR", "qusage")

	for ( pack in BIOC_Set) {
		if ( ! ( pack %in% packages)) {
			cat( "\nPackage '", pack, "' not installed yet.   Loading from Bioconductor..\n", sep="")
			biocLite( pack, suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
			success <- require( pack, character.only=TRUE)
			if ( success) {
				Nadd <- Nadd + 1
				addSet[Nadd] <- pack
			} else {
				Nfail <- Nfail + 1
				failSet[Nfail] <- pack
			}
		}
	}


	if ( Nadd + Nfail == 0) {
		cat( "\n\nAll packages used by DuffyTools are already installed.\n")
		return()
	}
	if ( Nadd > 0) {
		cat( "\n\nInstalled ", Nadd, " packages wanted by DuffyTools.\n")
		cat( addSet, "\n")
	}
	if ( Nfail > 0) {
		cat( "\n\nFailled to Instal ", Nfail, " packages wanted by DuffyTools.\n")
		cat( failSet, "\n")
	}
	if ( Nadd > 0 && Nfail == 0) cat( "\n\nAll packages used by DuffyTools are now installed.\n")
}
