# blastTools.R -- wrapper functions for the most command BLAST calls


`callBlast` <- function( fastafile, outfile="blastOut.txt", program="blastn", db="nt", 
			path=Sys.getenv("BLASTINDEX_PATH"), task="", 
			wordsize=8, evalue=1, threads=4, outfmt=6, 
			filter="no", maxhits=5, blastArgs="", verbose=T) {

	# validate the arguments
	useprogram <- Sys.which( program)
	if ( useprogram == "") stop( "Can't find executable: ", program)
	cmdline <- useprogram 
	if ( task != "") cmdline <- paste( cmdline, " -task", task)

	if ( ! file.exists( path)) stop( "Blast database path not found: ", path)
	usedb <- file.path( path, db)
	chkpattern <- paste( "^", basename(db), ".*[np]sq$", sep="")
	if ( length( dir( dirname(usedb), pattern=chkpattern)) < 1) stop( "Can't find database file: ", 
			usedb, "\nRegexpr pattern = ", chkpattern)
	cmdline <- paste( cmdline, "  -db ", usedb)

	if ( ! file.exists( fastafile)) stop( "Can't find FASTA query file: ", fastafile)

	cmdline <- paste( cmdline, " -query", fastafile, " -out", outfile)
	cmdline <- paste( cmdline, " -evalue", evalue, " -word_size", wordsize)

	if( regexpr( "blastp", program) > 0) {
		cmdline <- paste( cmdline, " -seg", filter)
	} else if( regexpr( "tblastn", program) > 0) {
		cmdline <- paste( cmdline, " -seg", filter)
	} else if( regexpr( "blastn", program) > 0) {
		cmdline <- paste( cmdline, " -dust", filter)
	}

	cmdline <- paste( cmdline, " -outfmt", outfmt, " -num_threads", threads)

	cmdline <- paste( cmdline, " -culling_limit", maxhits*2, " -num_alignments", maxhits)
	# some BLAST args changed their names...
	#cmdline <- paste( cmdline, " -max_hsps", maxhits, " -culling_limit", maxhits*2, " -num_alignments", maxhits)
	#cmdline <- paste( cmdline, " -max_hsps", maxhits, " -num_alignments", maxhits)

	# call BLAST
	if (verbose) cat( "\nCalling BLAST: \nCommand line:  ", cmdline, "\n")
	catch.system( cmdline)
	if (verbose) cat( "\nDone.\nWrote file: ", outfile, "\n")
	return()
}


`callBlastp` <- function( fastafile, outfile="blastOut.txt", program="blastp", db="NR/nr",
			path=Sys.getenv( "BLASTINDEX_PATH"), task="", wordsize=3, evalue=1,
			threads=4, outfmt=6, filter="no", maxhits=5, blastArgs="", verbose=T) {

	callBlast( fastafile, outfile=outfile, program=program, db=db, path=path, task=task, 
			wordsize=wordsize, evalue=evalue, threads=threads, outfmt=outfmt, 
			filter=filter, maxhits=maxhits, blastArgs, verbose=verbose)
	return()
}


`callBlastn` <- function( fastafile, outfile="blastOut.txt", program="blastn", db="NT/nt",
			path=Sys.getenv( "BLASTINDEX_PATH"), task="", wordsize=8, evalue=1,
			threads=4, outfmt=6, filter="no", maxhits=5, blastArgs="", verbose=T) {

	callBlast( fastafile, outfile=outfile, program=program, db=db, path=path, task=task, 
			wordsize=wordsize, evalue=evalue, threads=threads, outfmt=outfmt, 
			filter=filter, maxhits=maxhits, blastArgs, verbose=verbose)
	return()
}


`callBlastx` <- function( fastafile, outfile="blastOut.txt", program="blastx", db="NR/nr",
			path=Sys.getenv( "BLASTINDEX_PATH"), task="", wordsize=3, evalue=1,
			threads=4, outfmt=6, filter="no", maxhits=5, blastArgs="", verbose=T) {

	callBlast( fastafile, outfile=outfile, program=program, db=db, path=path, task=task, 
			wordsize=wordsize, evalue=evalue, threads=threads, outfmt=outfmt, 
			filter=filter, maxhits=maxhits, blastArgs, verbose=verbose)
	return()
}


`callTBlastn` <- function( fastafile, outfile="blastOut.txt", program="tblastn", db="NT/nt",
			path=Sys.getenv( "BLASTINDEX_PATH"), task="", wordsize=3, evalue=1,
			threads=4, outfmt=6, filter="no", maxhits=5, blastArgs="", verbose=T) {

	callBlast( fastafile, outfile=outfile, program=program, db=db, path=path, task=task, 
			wordsize=wordsize, evalue=evalue, threads=threads, outfmt=outfmt, 
			filter=filter, maxhits=maxhits, blastArgs, verbose=verbose)
	return()
}


`makeBlastDB` <- function( fastafile, dbname=sub(".fasta$","",basename(fastafile)), 
			dbtype=c("prot", "nucl"), path=Sys.getenv("BLASTINDEX_PATH"), 
			blastArgs="") {

	# validate the arguments
	dbtype <- match.arg( dbtype)
	useprogram <- Sys.which( "makeblastdb")
	if ( useprogram == "") stop( "Can't find executable: ", "makeblastdb")

	if ( ! file.exists( fastafile)) {
		stop( paste( "Input FASTA file not found: ", fastafile))
	}

	if ( ! file.exists( path)) {
		stop( paste( "BLAST database destination folder not found: ", path))
	}
	outname <- file.path( path, dbname)

	cmdline <- paste( useprogram, " -in ", fastafile, " -dbtype ", dbtype, 
			" -out ", outname, " ", blastArgs)
	catch.system( cmdline)
}


`blastMatchPercentage` <- function( v) {

	# decode a vector of character "(xx/yy)" blast match scores into a 0 to 100 score
	out <- vector( length=length(v))

	vv <- base::gsub( " ","",v)
	x <- as.numeric( base::sub( "(\\()([0-9]+)(/.*)", "\\2", vv))
	y <- as.numeric( base::sub( "(.*/)([0-9]+)(\\))", "\\2", vv))
	out <- 100 * x / y
	return( base::pmin( out, rep( 100.0, times=length(out))))
}


`readBlastOutput` <- function( infile, outfmt=6, verbose=TRUE, nKeep=NULL) {

	# scan the output of Blast for the data columns we need
	# depending on the version of BLAST, the column mode number is different
	if ( outfmt %in% c(6,8)) {
		what <- list( "PROBE_ID"="character", "SEQ_ID"="character", "PCT_MATCH"="numeric",
			"LEN_MATCH"="integer", "MIS_MATCH"="integer", "GAP"="integer", 
			"P_FIRST"="integer", "P_LAST"="integer", "S_BEG"="integer", "S_END"="integer", 
			"E_VALUE"="numeric", "SCORE"="integer")
	} else {
		warning( paste( "Unsupported Blast Output view format...:  outfmt=", outfmt))
		return( NULL)
	}

	blastData <- scan( file=infile, what=what, flush=TRUE, fill=TRUE, quiet=verbose )
	if ( length( blastData[[1]]) < 1) {
		cat( "\nProblem reading file:  ", basename(infile))
		cat( "\nNo records found.")
		return( data.frame())
	}

	# extract what we need
	pSet <- blastData$PROBE_ID
	sSet <- blastData$SEQ_ID
	pctMatch <- as.numeric( blastData$PCT_MATCH)
	lenMatch <- as.integer( blastData$LEN_MATCH)
	misMatch <- as.integer( blastData$MIS_MATCH)
	gap <- as.integer( blastData$GAP)
	firstPmatch <- as.integer( blastData$P_FIRST)
	lastPmatch <- as.integer( blastData$P_LAST)
	sBeg <- as.integer( blastData$S_BEG)
	sEnd <- as.integer( blastData$S_END)
	# force begin < end
	strand <- rep( "+", times=length( sBeg))
	strand[ sBeg > sEnd] <- "-"
	sBegNew <- base::pmin( sBeg, sEnd)
	sEndNew <- base::pmax( sBeg, sEnd)
	sBeg <- sBegNew
	sEnd <- sEndNew
	evalue <- as.numeric( blastData$E_VALUE)
	score <- as.numeric( blastData$SCORE)
	rm( blastData, what)

	blastDF <- data.frame( pSet, sSet, pctMatch, lenMatch, misMatch, gap, firstPmatch, lastPmatch, 
			sBeg, sEnd, evalue, score, strand, stringsAsFactors=FALSE)
	colnames( blastDF) <- c( "PROBE_ID", "SEQ_ID", "PCT_MATCH", "LEN_MATCH", "MIS_MATCH", "GAP", 
			"P_FIRST", "P_LAST", "S_BEG", "S_END", "E_VALUE", "SCORE", "STRAND")
	rownames(blastDF) <- 1:nrow(blastDF)

	if ( ! is.null(nKeep)) {

		# find and keep just the top N for each construct
		ans <- tapply( 1:nrow(blastDF), factor( blastDF$PROBE_ID), function(x) {
					if ( length(x) <= nKeep) return( x)
					ord <- order( blastDF$SCORE[x], decreasing=TRUE)
					return( x[ ord[ 1:nKeep]])
				}, simplify=FALSE)
		keep <- sort( unlist( ans))
		blastDF <- blastDF[ keep, ]
		rownames(blastDF) <- 1:nrow(blastDF)
	}

	return( blastDF)
}
