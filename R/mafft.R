
`mafft` <- function( fastaFile, outFile="mafft.aln", mafftProgram="~/bin/mafft",
			mode=c("global", "local", "genaffine"), anysymbol=FALSE,
			iterations=1000, outfmt=c("clustal", "fasta"), 
			guidetree=FALSE, mafftArgs="", verbose=FALSE) {

	if ( is.null( mafftProgram) || nchar( mafftProgram) < 1) {
		cat( "\nRequired 'mafftProgram' argument missing or empty...")
		return(NULL)
	}
	if ( ! file.exists( mafftProgram)) {
		cat( "\nMAFFT executable not found.  Tried: ", mafftProgram)
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
	anysymbolStr <- if (anysymbol) " --anysymbol" else ""

	# build the command line
	cmdLine <- paste( mafftProgram, modeStr, iterStr, treeStr, outfmtStr, anysymbolStr, 
				mafftArgs, fastaFile, " > ", outFile)

	# clean away the expected result
	file.delete( outFile)

	# call it!
	if (verbose) {
		catch.system( cmdLine)
	} else {
		sink( "MAFFT.log.txt")
		on.exit( sink())
		catch.system( cmdLine, ignore.stdout=FALSE, ignore.stderr=TRUE)
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

