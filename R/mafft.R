
`mafft` <- function( fastaFile, outFile="mafft.aln", mafftProgram="~/bin/mafft",
			mode=c("global", "local", "genaffine"), 
			iterations=1000, outfmt=c("clustal", "fasta"), 
			guidetree=FALSE, mafftArgs="", verbose=FALSE) {

	if ( is.null( mafftProgram) || nchar( mafftProgram) < 1) {
		cat( "\nRequired 'mafftProgram' argument missing or empty...")
		return(NULL)
	}

	mode <- match.arg( mode)
	outfmt <- match.arg( outfmt)
	modeStr <- " --globalpair "
	if ( mode == "local") modeStr <- " --localpair "
	if ( mode == "genaffine") modeStr <- " --genafpair "

	if ( is.null( iterations)) interations <- 1000
	iterStr <- paste( " --maxiterate", as.character(iterations))

	treeStr <- if (guidetree) " --treeout " else ""
	outfmtStr <- if ( outfmt == "fasta") "" else " --clustalout "

	# build the command line
	cmdLine <- paste( mafftProgram, modeStr, iterStr, treeStr, outfmtStr, mafftArgs, fastaFile, " > ", outFile)

	# clean away the expected result
	file.delete( outFile)

	# call it!
	if (verbose) {
		system( cmdLine)
	} else {
		sink( "MAFFT.log.txt")
		system( cmdLine, ignore.stdout=FALSE, ignore.stderr=TRUE)
		sink()
	}
	if ( ! file.exists( outFile)) {
		cat( "\nError:  No result made by MAFFT")
		return(NULL)
	}

	# consume the result
	ans <- NULL
	if ( outfmt == "clustal") {
		ans <- readALN( outFile, verbose=verbose)
	} else {
		ans <- loadFasta( outFile)
	}

	return( ans)
}

