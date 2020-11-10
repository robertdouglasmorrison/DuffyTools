# checkDuffyPackages.R -- verify needed packages...

`checkDuffyPackages` <- function() {

	# see what is installed
	current <- installed.packages()
	packages <- current[ ,"Package"]
	versions <- current[ ,"Version"]

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
		} else {
			where <- which( packages == pack)
			cat( "\nPackage '", pack, "' found.  \tVersion = ", versions[where],  sep="")
		}
	}


	# some are Bioconductor...
	BIOC_Set <- c( "Biostrings", "DESeq", "DESeq2", "edgeR", "siggenes", "ROC", "sangerseqR", "qusage")

	# as of R3.5, Bioconductor uses a newer method
	useBiocLite <- ( version$major == "2" || ( version$major == "3" && as.numeric( version$minor) < 5))
	if ( useBiocLite) {
		source( "http://bioconductor.org/biocLite.R")
	} else {
		if ( ! requireNamespace( "BiocManager", quietly=T)) install.packages( "BiocManager")
	}

	for ( pack in BIOC_Set) {
		if ( ! ( pack %in% packages)) {
			cat( "\nPackage '", pack, "' not installed yet.   Loading from Bioconductor..\n", sep="")
			if ( useBiocLite) {
				biocLite( pack, suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
			} else {
				BiocManager::install( pack, suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
			}
			success <- require( pack, character.only=TRUE)
			if ( success) {
				Nadd <- Nadd + 1
				addSet[Nadd] <- pack
			} else {
				Nfail <- Nfail + 1
				failSet[Nfail] <- pack
			}
		} else {
			where <- which( packages == pack)
			cat( "\nPackage '", pack, "' found.  \tVersion = ", versions[where],  sep="")
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
