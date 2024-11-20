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
	defaultWhat <- list( "PROBE_ID"="character", "SEQ_ID"="character", "PCT_MATCH"="numeric",
			"LEN_MATCH"="integer", "MIS_MATCH"="integer", "GAP"="integer", 
			"P_FIRST"="integer", "P_LAST"="integer", "S_BEG"="integer", "S_END"="integer", 
			"E_VALUE"="numeric", "SCORE"="integer")
	nOtherTerms <- 0
	if ( outfmt %in% c(6,8)) {
		what <- defaultWhat
	} else if ( grepl( "^6 std", outfmt)) {
		otherTerms <- toupper( strsplit( sub( "^6 std ", "", outfmt), split=" +")[[1]])
		nOtherTerms <- length( otherTerms)
		what <- defaultWhat
		if (nOtherTerms) for (k in 1:nOtherTerms) {
			what <- c( what, "character")
			names(what)[ length(what)] <- otherTerms[k]
		}
	} else {
		warning( paste( "Unsupported Blast Output view format...:  outfmt=", outfmt))
		return( NULL)
	}

	# single quote is very common in Blast descriptors, make sure that is not seen as a quote character by R
	blastData <- scan( file=infile, sep="\t", what=what, flush=TRUE, fill=TRUE, quote='"', quiet=verbose )
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
	if (nOtherTerms) {
		otherTextM <- matrix( "", nrow=length(pSet), ncol=nOtherTerms)
		for ( k in 1:nOtherTerms) otherTextM[ , k] <- as.character( blastData[[12+k]])
		colnames(otherTextM) <- otherTerms
	}
	rm( blastData)

	blastDF <- data.frame( pSet, sSet, pctMatch, lenMatch, misMatch, gap, firstPmatch, lastPmatch, 
			sBeg, sEnd, evalue, score, stringsAsFactors=FALSE)
	if ( nOtherTerms) blastDF <- cbind( blastDF, otherTextM, stringsAsFactors=F)
	colnames( blastDF) <- names(what)
	rownames(blastDF) <- 1:nrow(blastDF)
	blastDF$STRAND <- strand

	# if we got certain text string fields, there may be hypertext escape seqs
	if ("STITLE" %in% colnames(blastDF)) blastDF$STITLE <- convertHypertext( blastDF$STITLE)

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


extractBlastXMLdetails <- function( filein, IDs, IDprefix="", IDsuffix="") {

	# we want to extract the details that show the 'best hits'...
	out <- vector( mode="list")
	txt <- readLines( con=filein)
	iterDefLines <- grep( "Iteration_query-def", txt, fixed=TRUE)
	if ( length( iterDefLines) < 1) return( out)

	for( i in 1:length(iterDefLines)) {
		if ( i < length( iterDefLines)) {
			lastLine <- iterDefLines[i+1]
		} else {
			lastLine <- length( txt)
		}
		out[[i]] <- extractOneXMLiterationSet( txt[ iterDefLines[i] : lastLine] )
	}

	# some CRs may not return anything, but try to make the answer more complete
	idsOut <- sapply( out, function(x) x$id)
	finalOut <- vector( mode="list")
	for ( j in 1:length(IDs)) {
		i <- IDs[j]
		thisID <- paste( IDprefix, i, IDsuffix, sep="")
		where <- match( thisID, idsOut, nomatch=0)
		if ( where > 0) {
			finalOut[[j]] <- out[[where]]
		} else {
			finalOut[[j]] <- list( "n"=0, "id"=thisID, "accession"=NA, "definition"=NA, "evalue"=NA, 
						"score"=NA, "match"=NA, "hitFrom"=NA, "hitTo"=NA)
		}
	}

	return( finalOut)
}


extractOneXMLiterationSet <- function( txt) {

	# given a subset of text lines that contain one XML blast output alignment...
	id <- sub( "( *<Iteration_query-def>)(.+)(</Iteration_query-def>)", "\\2", txt[1])
	hitLines <- grep( "<Hit>", txt, fixed=TRUE)
	nHits <- length( hitLines)
	if ( nHits < 1) return( list( "n"=0, "id"=id))

	hitACClines <- grep( "<Hit_accession>", txt, fixed=TRUE)
	hitACCs <- sub( "( *<Hit_accession>)(.+)(</Hit_accession>)", "\\2", txt[ hitACClines])
	hitDEFlines <- grep( "<Hit_def>", txt, fixed=TRUE)
	hitDEFs <- sub( "( *<Hit_def>)(.+)(</Hit_def>)", "\\2", txt[ hitDEFlines])
	
	# use just the first HSP for each hit, static offsets from hit line...
	hitEvalue <- hitScore <- hitMatch <- hitFrom <- hitTo <- vector()
	for ( i in 1:length( hitLines)) {
		line <- hitLines[i]
		hitScore[i] <- as.numeric( sub( "( +<Hsp_bit-score>)(.+)(</Hsp_bit-score>)", "\\2", txt[line+9]))
		hitEvalue[i] <- as.numeric( sub( "( +<Hsp_evalue>)(.+)(</Hsp_evalue>)", "\\2", txt[line+11]))
		qfrom <- as.integer( sub( "( +<Hsp_query-from>)(.+)(</Hsp_query-from>)", "\\2", txt[line+12]))
		qto <- as.integer( sub( "( +<Hsp_query-to>)(.+)(</Hsp_query-to>)", "\\2", txt[line+13]))
		hfrom <- as.integer( sub( "( +<Hsp_hit-from>)(.+)(</Hsp_hit-from>)", "\\2", txt[line+14]))
		hto <- as.integer( sub( "( +<Hsp_hit-to>)(.+)(</Hsp_hit-to>)", "\\2", txt[line+15]))
		hitFrom[i] <- hfrom
		hitTo[i] <- hto

		# there may be a 'gaps' line... there may not be...
		hasGaps <- ( regexpr( "Hsp_gaps", txt[ line+20], fixed=TRUE) > 0)
		if ( hasGaps) line <- line + 1

		seq <- sub( "( +<Hsp_qseq>)(.+)(</Hsp_qseq>)", "\\2", txt[line+21])
		qtxt <- base::paste( formatC( qfrom, digits=12, format="d"), seq, formatC( qto, digits=12, 
				format="d"), sep=" ")
		seq <- sub( "( +<Hsp_hseq>)(.+)(</Hsp_hseq>)", "\\2", txt[line+22])
		htxt <- base::paste( formatC( hfrom, digits=12, format="d"), seq, formatC( hto, digits=12, 
				format="d"), sep=" ")
		seq <- sub( "( +<Hsp_midline>)(.+)(</Hsp_midline>)", "\\2", txt[line+23])
		mtxt <- base::paste( "             ", seq, "             ", sep=" ")
		hitMatch[i] <- base::paste( qtxt, mtxt, htxt, sep="\n")
	}

	out <- list( "n"=nHits, "id"=id, "accession"=hitACCs, "definition"=hitDEFs, "evalue"=hitEvalue, 
			"score"=hitScore, "match"=hitMatch, "hitFrom"=hitFrom, "hitTo"=hitTo)

	return( out)
}

