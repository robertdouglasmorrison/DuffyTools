# checkPackages.R -- verify needed packages...

`checkPackages` <- function() {

	current <- installed.packages()
	packages <- current[ ,"Package"]
	Nmiss <- 0


	# some are regular CRAN packages
	if ( ! ( "RCurl" %in% packages)) {
		cat( "\nPackage 'RCurl' not installed yet.   Loading from CRAN..\n")
		try( install.packages( "RCurl"))
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "RUnit" %in% packages)) {
		cat( "\nPackage 'RUnit' not installed yet.   Loading from CRAN..\n")
		try( install.packages( "RUnit"))
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "plotrix" %in% packages)) {
		cat( "\nPackage 'plotrix' not installed yet.   Loading from CRAN..\n")
		try( install.packages( "plotrix"))
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "gplots" %in% packages)) {
		cat( "\nPackage 'gplots' not installed yet.   Loading from CRAN..\n")
		try( install.packages( "gplots"))
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "heatmap.plus" %in% packages)) {
		cat( "\nPackage 'heatmap.plus' not installed yet.   Loading from CRAN..\n")
		try( install.packages( "heatmap.plus"))
		Nmiss <- Nmiss + 1
	}

	# some are Bioconductor...
	source( "http://bioconductor.org/biocLite.R")
	if ( ! ( "Biostrings" %in% packages)) {
		cat( "\nPackage 'Biostrings' not installed yet.   Loading from Bioconductor..\n")
		biocLite( "Biostrings", suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "DESeq" %in% packages)) {
		cat( "\nPackage 'DESeq' not installed yet.   Loading from Bioconductor..\n")
		biocLite( "DESeq", suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "edgeR" %in% packages)) {
		cat( "\nPackage 'edgeR' not installed yet.   Loading from Bioconductor..\n")
		biocLite( "edgeR", suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "ROC" %in% packages)) {
		cat( "\nPackage 'ROC' not installed yet.   Loading from Bioconductor..\n")
		biocLite( "ROC", suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
		Nmiss <- Nmiss + 1
	}

	if ( ! ( "siggenes" %in% packages)) {
		cat( "\nPackage 'siggenes' not installed yet.   Loading from Bioconductor..\n")
		biocLite( "siggenes", suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
		Nmiss <- Nmiss + 1
	}

	if ( Nmiss > 0) {
		cat( "\n\nInstalled ", Nmiss, " packages used by DuffyTools.\n")
	} else {
		cat( "\n\nAll packages used by DuffyTools are installed.\n")
	}
}
