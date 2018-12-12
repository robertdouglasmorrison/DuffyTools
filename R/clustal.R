`clustal` <- function( fastaFile, outFile="clustal.aln", clustalProgram=Sys.getenv("CLUSTAL"), 
			iterations=NULL, outfmt=c("clustal", "fasta", "phylip"), 
			guidetree=NULL, clustalArgs=NULL, verbose=FALSE) {

	if ( is.null( clustalProgram) || nchar( clustalProgram) < 1) {
		cat( "\nRequired 'clustalProgram' argument missing or empty...")
		return(NULL)
	}

	outfmt <- match.arg( outfmt)

	# build the command line
	cmdLine <- paste( clustalProgram, " -i ", fastaFile, " -o ", outFile,
				" --outfmt ", outfmt, "  --force")
	if ( ! is.null( iterations)) cmdLine <- paste( cmdLine, "  --iterations ", iterations)
	if ( ! is.null( guidetree)) cmdLine <- paste( cmdLine, "  --guidetree-out ", guidetree)
	if ( ! is.null( clustalArgs)) cmdLine <- paste( cmdLine, clustalArgs)
	if (verbose) cmdLine <- paste( cmdLine, "  -v")

	# clean away the expected result
	file.delete( outFile)

	# call it!
	system( cmdLine)
	if ( ! file.exists( outFile)) {
		cat( "\nError:  No result made by CLUSTAL")
		return(NULL)
	}

	# consume the result
	ans <- NULL
	if ( outfmt == "clustal") ans <- readALN( outFile, verbose=verbose)

	return( ans)
}

