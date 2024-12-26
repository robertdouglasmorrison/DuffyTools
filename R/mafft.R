# mafft.R -- wrapper for doing MSA by MAFFT

`mafft` <- function( fastaFile, outFile="mafft.aln", mafftProgram=Sys.which("mafft"),
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


`fastMafft` <- function( fastaFile, outFile="mafft.aln", mafftProgram=Sys.which("mafft"),
			mode=c("global", "local", "genaffine"), anysymbol=FALSE,
			iterations=1000, outfmt=c("clustal", "fasta"), 
			guidetree=FALSE, mafftArgs="", verbose=FALSE) {

	# this is intended to speed up MAFFT for large sequence sets, when many of the sequences
	# are exact duplicates.  By sending MAFFT just the uniquely different ones and then
	# rebuilding what the full dataset would look like
	
	doFAST <- TRUE
	# step 1a: pre-read the FASTA and decide which sequences to use
	fa <- loadFasta( fastaFile, short=F, verbose=F)
	fullSeqSet <- fa$seq
	fullDescSet <- fa$desc
	nFull <- length(fullSeqSet)
	seqTbl <- table( fullSeqSet)
	nUnique <- length(seqTbl)
	if (nUnique > (nFull*0.8)) doFAST <- FALSE
	seqPcts <- as.numeric(seqTbl) * 100 / nFull
	bigPct <- max( seqPcts)
	
	# set 1b: find the number of copies of each to retain.  Because MSA is consensus based, make sure
	# we keep the copies of the small set proportional to the copies of the original
	MAX_COPIES <- 20
	nKeep <- round( seqPcts * MAX_COPIES / bigPct)
	nKeep[ nKeep < 1] <- 1
	if ( sum(nKeep) > (nFull*0.8)) doFAST <- FALSE
	
	# step 1c: make that small set
	if (doFAST) {
		smlSeqSet <- rep( names(seqTbl), times=nKeep)
		smlDescSet <- fullDescSet[ match( smlSeqSet,fullSeqSet)] 
		smlFastaFile <- paste( fastaFile, "Fast.Tmp.fasta", sep=".")
		smlOutFile <- paste( outFile, "Fast.Tmp.aln", sep=".")
		writeFasta( as.Fasta( smlDescSet, smlSeqSet), smlFastaFile, line=100)
	}
	
	# step 2: call MAFFT
	if ( doFAST) {
		smlAns <- mafft( smlFastaFile, smlOutFile, mafftProgram=mafftProgram,
				mode=mode, anysymbol=anysymbol,
				iterations=iterations, outfmt=outfmt, 
				guidetree=guidetree, mafftArgs=mafftArgs, verbose=verbose)
	} else {
		# not enough benefit, just to the usual
		return( mafft( fastaFile, outFile, mafftProgram=mafftProgram,
				mode=mode, anysymbol=anysymbol,
				iterations=iterations, outfmt=outfmt, 
				guidetree=guidetree, mafftArgs=mafftArgs, verbose=verbose))
	}
	
	# step 3: expand the result to look like it was full size
	# note that MAFFT never rearranges the order of sequences, so we know where they all ended up
	smlALNmatrix <- smlAns$alignment
	wh <- match( fullSeqSet, smlSeqSet)
	fullALNmatrix <- smlALNmatrix[ wh, ]
	rownames(fullALNmatrix) <- fullDescSet
	# write it out as the original result
	out <- smlAns
	out$alignment <- fullALNmatrix
	writeALN( out, outFile)
	# remove our temp files
	file.delete( c( smlFastaFile, smlOutFile))
	
	return( out)
}
