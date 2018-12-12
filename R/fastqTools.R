# fastqTools.R


`readFastq` <- function( filein, maxReads=NULL, verbose=TRUE) {

	fileToUse <- allowCompressedFileName( filein)

	if ( ! file.readable( fileToUse)) {
		warning( paste( "readFastq:  unable to read FASTQ file", fileToUse))
		return( data.frame())
	}

	# catch compressed files
	conIn <- openCompressedFile( fileToUse, open="r")
	if (verbose) cat( "\nReading file: ", fileToUse)

	# read in the raw text
	chunkSize <- 1000000
	if ( ! is.null(maxReads) && (maxReads*4) < chunkSize) chunkSize <- maxReads * 4

	readIDs <- readSeqs <- scores <- vector()

	repeat {

		fastqText <- readLines( con=conIn, n=chunkSize, warn=FALSE)
		if ( length( fastqText) < 1) break

		# get the delimited lines
		idLines <- grep( "^@", fastqText)
		scoreMarks <- grep( "^\\+", fastqText)
		if ( length( idLines) < 1 || length(scoreMarks) < 1 ) {
			warning( paste( "readFastq:  not a FASTQ format file: ", fileToUse))
			return( data.frame())
		}

		idLines <- base::sort( intersect( idLines, seq( 1, length(fastqText), by=4)))
		scoreMarks <- base::sort( intersect( scoreMarks, seq( 3, length(fastqText), by=4)))
		if ( (length( idLines) != length(scoreMarks)) || any( idLines >= scoreMarks)) {
			warning( paste( "readFastq:  bad FASTQ file: ", fileToUse, 
					"\tMismatched readIDs and/or quality scores..."))
			return( data.frame())
		}

		# now get the parts
		readIDs <- base::append( readIDs, sub( "^@", "", fastqText[ idLines]))
		readSeqs <- base::append( readSeqs, fastqText[ (idLines+1)])
		scores <- base::append( scores, fastqText[ (scoreMarks+1)])
	
		if (verbose) cat(".")
		if ( ! is.null( maxReads)) {
			if ( length(readIDs) >= maxReads) break
		}
	} # end of 'repeat' loop...
	close( conIn)
	
	if ( any( base::nchar( readSeqs) != base::nchar( scores))) 
			warning( "readFastq:  some Read & Score lengths disagree")

	out <- data.frame( readIDs, readSeqs, scores, stringsAsFactors=FALSE)
	colnames( out) <- FASTQ_COLUMNS
	rownames( out) <- 1:nrow(out)

	if ( ! is.null( maxReads)) {
		if ( nrow(out) >= maxReads) out <- out[ 1:maxReads, ]
	}

	if (verbose) cat( "\nN_Reads: ", prettyNum( nrow(out), big.mark=","), "\n")
	return( out)
}


`test.readFastq` <- function() {

	tmpFile <- build.testFastq()
	zz <- readFastq( tmpFile, verbose=FALSE)
	checkEquals( dim(zz), c(1000,3))
	checkEquals( colnames(zz), FASTQ_COLUMNS)
	remove.testFile( tmpFile)
}


# efficient writing of .fastq to a file

`writeFastq` <- function( x, fileout, compress=FALSE, verbose=TRUE) {

	if ( ! all( colnames( x) == FASTQ_COLUMNS)) {
		warning( paste( "writeFastq:  unexpected column names: ", colnames(x)))
		warning( "No file written...")
		return( NULL)
	}

	fileToUse <- fileout

	# catch compressed file
	hasGZ <- (regexpr( "\\.gz$", fileout) > 0)
	if ( compress || hasGZ) {
		if ( ! hasGZ) fileToUse <- paste( fileout, "gz", sep=".")
		conOut <- gzfile( fileToUse, open="w")
	} else {
		conOut <- file( fileToUse, open="w")
	}

	# prepend the .fastq markers to the needed spots
	idText <- base::paste( "@", x$READ_ID, sep="")
	out <- base::paste( idText, x$READ_SEQ, "+", x$SCORE, sep="\n")

	# write it out
	writeLines( out, con=conOut, sep="\n")
	if( verbose) cat( "\nWrote file: \t", fileToUse, "\nN_Reads: \t", nrow( x),"\n")
	close( conOut)
	return( fileToUse)
}


`test.writeFastq` <- function() {

	tmpFile <- build.testFastq()
	zz <- readFastq( tmpFile, verbose=FALSE)
	tmpFile2 <- writeFastq( zz, fileout="DuffyTools.test2.fastq", verbose=FALSE)
	zz2 <- readFastq( tmpFile2, verbose=FALSE)
	checkEquals( zz, zz2)
	remove.testFile( tmpFile)
	remove.testFile( tmpFile2)
}


# make shorter reads from the ends for splice discovery

`fastqToOverlapReadlets` <- function( filein, fileout, segmentLength=12, verbose=TRUE) {

	fileToUse <- allowCompressedFileName( filein)
	if ( ! file.exists( fileToUse)) stop( paste("Can't find input file: ", fileToUse))
	conIn <- openCompressedFile( fileToUse, open="r")

	if ( regexpr( ".gz$", fileout) > 0) {
		conOut <- gzfile( fileout, open="w")
	} else {
		conOut <- file( fileout, open="w")
	}

	chunkSize <- 400000
	nread <- 0
	if (verbose) cat( "\nReading file: ", fileToUse)

	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		# get the lines we want
		ids <- chunk[ seq( 1, length(chunk), by=4)]
		seqs <- chunk[ seq( 2, length(chunk), by=4)]
		scores <- chunk[ seq( 4, length(chunk), by=4)]
		N <- length( ids)

		# get the actual read lengths
		if ( nread == 0) {
			alllens <- base::nchar( seqs)
			maxlen <- max( alllens)
			segmentLength <- ceiling(segmentLength/2) * 2
			overlap <- segmentLength / 2
			nOverlaps <- ceiling( (maxlen - segmentLength) / overlap) + 1
			if (verbose) cat( "\n readLen, segLength, overlap, nSegs: ", maxlen, segmentLength, overlap, nOverlaps)
		}

		myStarts <- seq( 1, (maxlen-overlap), by=overlap)
		myStops <- seq( segmentLength, (maxlen+overlap-1), by=overlap)
		if ( length( myStarts) != nOverlaps || length( myStarts) != nOverlaps) {
			cat( "bad overlap sizes...?")
			stop("")
		}
		myStarts[nOverlaps] <- maxlen - segmentLength + 1
		myStops[nOverlaps] <- maxlen

		# get the substring chunks
		allSeqs <- allScores <- allIDs <- matrix( nrow=N, ncol=nOverlaps)
		for ( k in 1:nOverlaps) {
			allSeqs[ ,k] <- base::substr( seqs, myStarts[k], myStops[k])
			allScores[ ,k] <- base::substr( scores, myStarts[k], myStops[k])
			allIDs[ ,k] <- base::paste( ids, "/seg", k, sep="")
		}

		# interleave the output to keep the pairs together
		outIDs <- as.vector( t( allIDs))
		outSeqs <- as.vector( t( allSeqs))
		outScores <- as.vector( t( allScores))

		# format as Fastq (the ID still has the '@'...)
		outTxt <- base::paste( outIDs, "\n", outSeqs, "\n+\n", outScores, sep="")
		writeLines( outTxt, con=conOut)

		nread <- nread + length( chunk)

		if (verbose) cat( ".")
	}

	close( conIn)
	close( conOut)
	if (verbose) cat( "\nN_reads:            ", round( as.integer(nread)/4))
	if (verbose) cat( "\nWrote file:         ", fileout, "\n")

	return( list( "Segments"=nOverlaps, "Overlap"=overlap))
}


`fastqToChunks` <- function( filein, fileroot=sub( "(fq|fastq)(.gz)?$","",filein), 
				filetail=sub( fileroot, "", filein, fixed=T), 
				chunk.size=1000000) {

	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")


	blockSize <- min( chunk.size, 500000) * 4
	nread <- 0
	nchunk <- 0
	cat( "\nBreaking: ", basename(filein), " into chunks of ", chunk.size, "reads.\n")

	nThisChunk <- 0
	makeNewFile <- TRUE

	repeat {

		if ( makeNewFile) {
			nchunk <- nchunk + 1
			fileout <- paste( fileroot, "chunk", nchunk, ".", filetail, sep="")
			if ( regexpr( ".gz$", fileout) > 0) {
				conOut <- gzfile( fileout, open="w")
			} else {
				conOut <- file( fileout, open="w")
			}
			makeNewFile <- FALSE
		}

		thisBlockSize <- min( blockSize, (chunk.size-nThisChunk)*4)
		chunk <- readLines( conIn, n=thisBlockSize)
		if ( length( chunk) < 1) break

		nread <- nread + length(chunk)/4
		nThisChunk <- nThisChunk + length(chunk)/4

		# get the lines we want
		who <- seq( 1, length(chunk), by=4)
		ids <- chunk[ who]
		seqs <- chunk[ who + 1]
		scores <- chunk[ who + 3]
		writeLines( base::paste( ids, seqs, "+", scores, sep="\n"), con=conOut)
		cat( ".")

		if ( nThisChunk >= chunk.size) {
			close( conOut)
			cat( " ", nread, basename(fileout), "\n")
			makeNewFile <- TRUE
			nThisChunk <- 0
		}
	}
	
	# close the last one...
	close( conOut)
	cat( "\n", nread, basename(fileout))

	close( conIn)
	cat( "\nN_reads:            ", round( as.integer(nread)/4))
	return()
}


`fastqMerge` <- function( filesin, fileout, verbose=TRUE) {

	# use system level zcat to do the merge
	cmdText <- "zcat -f "
	infileText <- paste( filesin, collapse=" ")
	
	doZip <- (regexpr( "gz$", fileout) > 0)
	if ( doZip) {
		cmdText <- paste( cmdText, infileText, " | gzip > ", fileout, sep=" ")
	} else {
		cmdText <- paste( cmdText, infileText, " > ", fileout, sep=" ")
	}

	if (verbose) cat( "\nMerging: \n  ", cmdText)
	system( cmdText)
	if (verbose) cat( "\nDone.\n")
}


`fastqToFasta` <- function( filein, fileout, Qscores=TRUE) {

	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	if ( regexpr( ".gz$", fileout) > 0) {
		conOut <- gzfile( fileout, open="w")
	} else {
		conOut <- file( fileout, open="w")
	}

	chunkSize <- 800000
	nread <- 0

	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		nread <- nread + length( chunk)

		# get the lines we want
		ids <- chunk[ seq( 1, length(chunk), by=4)]
		seqs <- chunk[ seq( 2, length(chunk), by=4)]

		if ( Qscores) {
			scores <- chunk[ seq( 4, length(chunk), by=4)]
			avgScores <- apply( phredScoreStringToInt( scores), MARGIN=1, FUN=mean)
			newIds <- base::paste( sub( "@", ">", ids, fixed=TRUE), ":Q=", formatC( avgScores, digits=2, format="f"), sep="")
		} else {
			newIds <- sub( "@", ">", ids, fixed=TRUE)
		}
		writeLines( base::paste( newIds, seqs, sep="\n"), con=conOut)
		cat( ".")
	}

	close( conIn)
	close( conOut)
	cat( "\nN_reads:            ", round( as.integer(nread)/4))

	return()
}


`clipBases` <- function( mySEQ, myScores, laneID, clip5prime, clip3prime, scoresAsIntegers=TRUE) {

	# parse the clip instructions...
	if ( is.null( clip5prime)) clip5prime <- 0
	if ( is.null( clip3prime)) clip3prime <- 0
	clipLanes <- vector()
	if ( ! is.null( names(clip5prime))) clipLanes <- names( clip5prime)
	if ( ! is.null( names(clip3prime))) clipLanes <- unique.default( base::append( clipLanes, names( clip3prime)))
	nClipLanes <- length( clipLanes)
	doByLane <- (nClipLanes > 0)
	if ( doByLane) {
		# there were some named lanes, so expand to a complete set
		trims <- matrix( 0, nrow=max(as.numeric(clipLanes)), ncol=2)
		trims[ as.numeric( names( clip5prime)), 1] <- base::unlist( clip5prime)
		trims[ as.numeric( names( clip3prime)), 2] <- base::unlist( clip3prime)
	} 

	# get the lengths
	nlines <- length( mySEQ)
	nbases <- base::nchar( mySEQ)
	nscoreCh <- base::nchar( myScores)

	# and remove the leading spaces before the first base score
	if ( scoresAsIntegers) {
	  repeat {
		hasBlank <- which( base::substr( myScores,1,1) == " ")
		if ( length( hasBlank) < 1) break
		myScores[ hasBlank] <- sub( " ", "", myScores[ hasBlank], fixed=TRUE)
		nscoreCh[ hasBlank] <- nscoreCh[ hasBlank] - 1
	  }
	  scoreList <- strsplit( myScores, split=" ", fixed=TRUE)
	  lens <- sapply( scoreList, length)
	  if ( any( lens != nbases)) warning("clipBases:  Read -- Quality score length mis-match")
	}

	# ready to do the clipping...
	if ( ! doByLane) {
		cat( " clipping all lanes(5',3')=", clip5prime, clip3prime)
		new5 <- 1 + clip5prime
		for ( i in 1:nlines) {
			new3 <- nbases[i] - clip3prime
			mySEQ[i] <- base::substr( mySEQ[i], new5, new3)
			if ( scoresAsIntegers) {
				scoretmp <- scoreList[[i]][ new5:new3]
				myScores[i] <- base::paste( scoretmp, collapse=" ")
			} else {
				myScores[i] <- base::substr( myScores[i], new5, new3)
			}
		}

	} else {
		for ( ilane in 1:nClipLanes) {
			lane <- clipLanes[ ilane]
			this5 <- trims[ as.integer( lane),1]
			this3 <- trims[ as.integer( lane),2]
			if ( all( c( this5, this3) == 0)) next
			who <- which( laneID == lane)
			if ( length( who) < 1) next
			new5 <- 1 + this5
			cat( " clipping lane",lane,"(5',3')=", this5, this3)
			for ( i in who) {
				new3 <- nbases[i] - this3
				mySEQ[i] <- base::substr( mySEQ[i], new5, new3)
				if ( scoresAsIntegers) {
					scoretmp <- scoreList[[i]][ new5:new3]
					myScores[i] <- base::paste( scoretmp, collapse=" ")
				} else {
					myScores[i] <- base::substr( myScores[i], new5, new3)
				}
			}
		}
	}

	# ready.
	out <- list( "seqs"=mySEQ, "scores"=myScores)
	return( out)
}



clipFastq <- function( filein, fileout, clip5prime=0, clip3prime=0) {

	# clip bases off an existing fastq file

	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	if ( regexpr( ".gz$", fileout) > 0) {
		conOut <- gzfile( fileout, open="w")
	} else {
		conOut <- file( fileout, open="w")
	}

	chunkSize <- 800000
	nread <- 0

	cat( "\nClipping  ( 5', 3') = ", clip5prime, clip3prime, "\n")
	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		nread <- nread + length( chunk)

		# get the lines we want
		ids <- chunk[ seq( 1, length(chunk), by=4)]
		seqs <- chunk[ seq( 2, length(chunk), by=4)]
		scores <- chunk[ seq( 4, length(chunk), by=4)]

		ans <- clipBases( seqs, scores, ids, clip5=clip5prime, clip3=clip3prime, scoresAsIntegers=FALSE)

		newSeqs <- ans$seq
		newScores <- ans$scores

		newLen <- base::nchar( newSeqs[1])
		newIds <- sub( "=[0-9]+$", base::paste( "=",newLen, sep=""), ids) 
		newId2 <- sub( "@","+", newIds, fixed=TRUE)
		writeLines( base::paste( newIds, newSeqs, newId2, newScores, sep="\n"), con=conOut)
		cat( ".")
	}

	close( conIn)
	close( conOut)
	cat( "\nN_reads trimmed:            ", round( as.integer(nread)/4))
	cat( "\nNew Read Length:  ", base::nchar( newSeqs[1]), "\n")

	return()
}


`fastqPatternSearch` <- function( filein, patterns, max.mismatch=0, chunkSize=400000) {

	require( Biostrings)

	fileToUse <- allowCompressedFileName( filein)
	if ( ! file.exists( fileToUse)) stop( paste("Can't find input file: ", fileToUse))
	conIn <- openCompressedFile( fileToUse, open="r")

	nread <- 0
	cat( "\nReading file: ", fileToUse, "\n")

	nPatt <- length( patterns)
	findCounts <- rep( 0, times=nPatt)

	repeat {
		chunk <- readLines( conIn, n=chunkSize)
		if ( length( chunk) < 1) break

		# get the lines we want
		ids <- chunk[ seq( 1, length(chunk), by=4)]
		seqs <- chunk[ seq( 2, length(chunk), by=4)]
		N <- length( ids)

		# turn these into a serchable string
		subject <- DNAString( paste( seqs, collapse="N"))
		for ( k in 1:nPatt) {
			v <- countPattern( patterns[k], subject, max.mismatch=max.mismatch)
			findCounts[k] <- findCounts[k] + v
			if ( k %% 20 == 0) cat( "\r", k, v, sum( findCounts[1:k]))
		}
		nread <- nread + N
		cat( "\nReads: ", formatC( as.integer(nread), big.mark=","), "\tHits: ", sum( findCounts))
	}
	cat( "\n")

	close( conIn)

	out <- findCounts
	names(out) <- patterns
	return( out)
}



`bam2fastq` <- function( bamfile, outfile=sub( ".bam$", "", bamfile), paired.end=TRUE, clobber=TRUE) {


	cat( "\nConverting BAM file: ", bamfile)
	cat( "\nCreating FASTQ(s):   ", outfile, "\n")
	cmdline <- paste( "bam2fastq.pl ", " --filter '-F 0x000'   --prefix ", outfile, "  ", bamfile)
	if ( clobber) {
		cmdline <- paste( "bam2fastq.pl ", " --filter '-F 0x000'  --yes   --prefix ", outfile, "  ", bamfile)
	}

	system( cmdline)

	cat( "  Done.\n")
	return( )
}


`fastqToPeptides` <- function( filein, fileout, chunkSize=100000, maxReads=NULL, maxPeptides=NULL,
				lowComplexityFilter=FALSE, trim5=0, trim3=0, ...) {

	filein <- allowCompressedFileName( filein)
	if ( ! file.exists( filein)) stop( paste("Can't find input file: ", filein))
	conIn <- openCompressedFile( filein, open="r")

	chunkLines <- chunkSize * 4

	LAPPLY <- base::lapply
	NCHAR <- base::nchar
	SETDIFF <- base::setdiff
	WHICH <- base::which

	# local function to parallelize
	myDNAtoAAfunc <- function( x, peps, cnts) {

		peptides <- counts <- vector()
		nout <- 0

		LAPPLY( x, function(i) {
				dna <- peps[ i]
				# they have to be full length with no N's or non-AA calls
				# 'full length is a bit too restrictive...  we lose many peptides
				# at the true stop codon of the protein...  Relax a bit.
				# keep any where the STOP was in the last 25% of the peptide
				lenFullAA <- floor( NCHAR(dna)/3 * 0.75)
				AAs <- DNAtoAA( dna)
				good <- SETDIFF( WHICH( NCHAR(AAs) >= lenFullAA), grep( "?", AAs, fixed=T))
				if ( (ngood <- length(good)) > 0) {
					nnow <- (nout+1) : (nout+ngood)
					peptides[ nnow] <<- AAs[ good]
					counts[ nnow] <<- cnts[i]
					nout <<- nout + ngood
				}
				return(NULL)
			})

		# the mapping from DNA to AA may have generated duplicates
		if ( length( peptides) < 1) return( data.frame())
		tapply( 1:length(peptides), factor( peptides), FUN=function(x) {
				if ( length(x) < 2) return()
				totalCnt <- sum( counts[x])
				counts[ x[1]] <<- totalCnt
				counts[ x[2:length(x)]] <<- 0
				return()
			}, simplify=FALSE)
		keep <- which( counts > 0)

		# final ans is a table of peptides with counts
		ans <- list( "Peptide"=peptides[keep], "Count"=counts[keep])
		return(ans)
	}  # end of local function


	nread <- ntotal <- nUtotal <- 0
	repeat {
		chunk <- readLines( conIn, n=chunkLines)
		if ( length( chunk) < 1) break

		nReadsNow <- length(chunk)/4
		nread <- nread + nReadsNow
		cat( "\nReads: ", prettyNum( as.integer( nread), big.mark=","))

		# get the raw read lines we want
		seqs <- chunk[ seq.int( 2, length(chunk), 4)]

		# allow trimming
		if ( trim5 > 0) {
			cat( " trim5:", trim5)
			N <- max( nchar(seqs))
			seqs <- substr( seqs, (trim5+1), N)
		}
		if ( trim3 > 0) {
			cat( " trim3:", trim3)
			N <- max( nchar(seqs))
			seqs <- substr( seqs, 1, (N-trim3))
		}

		# but only do one of each duplicate
		seqCntsTbl <- table( seqs)
		uniqSeqs <- names( seqCntsTbl)
		uniqCnts <- as.vector( seqCntsTbl)
		nReadsNow <- N <- length(uniqCnts)
		rm( seqCntsTbl)
		nUtotal <- nUtotal + N
		cat( "  Unique: ", prettyNum( as.integer( nUtotal)))

		coreOpt <- getOption("cores")
		nCores <- if ( is.null(coreOpt)) 1 else as.integer( coreOpt)
		if ( nCores < 2 || N < 100) {
			doMULTICORE <- FALSE
			ans <- myDNAtoAAfunc( 1:N, uniqSeqs, uniqCnts)
		} else {
			doMULTICORE <- TRUE
			starts <- seq.int( 1, nReadsNow, round(nReadsNow/nCores))
			stops <- c( starts[2:length(starts)] - 1, nReadsNow)
		        cat( "  use", nCores, "cores..")
			locs <- vector( mode="list")
			for ( i in 1:length(starts)) locs[[i]] <- starts[i] : stops[i]

			ans <- multicore.lapply( locs, FUN=myDNAtoAAfunc, peps=uniqSeqs, cnts=uniqCnts, preschedule=TRUE)
		}

		# extract all the little data frame answers
		bigp <- bigc <- vector()
		if ( ! doMULTICORE) {
			bigp <- ans$Peptide
			bigc <- ans$Count
		} else {
		    cat( "  merge..")
		    for ( k in 1:length(ans)) {
			obj <- ans[[k]]
			lenNow <- length(obj$Peptide)
			if ( lenNow < 1) next
			newlocs <- (length(bigp) + 1) : (length(bigp) + lenNow)
			bigp[ newlocs] <- obj$Peptide
			bigc[ newlocs] <- obj$Count
		    }
		}
		N <- length(bigp)
		if ( ! N) next

		# the mapping from DNA to AA may have generated duplicates
		if ( N > 1) {
			cat( "  dups..")
			tapply( 1:N,  factor( bigp), FUN=function(x) {
					if ( length(x) < 2) return()
					totalCnt <- sum( bigc[x])
					bigc[ x[1]] <<- totalCnt
					bigc[ x[2:length(x)]] <<- 0
					return()
				}, simplify=FALSE)
		}
		keep <- which( bigc > 0)
		ansDF <- data.frame( "Peptide"=bigp[keep], "Count"=bigc[keep], stringsAsFactors=F)
		nout <- nrow(ansDF)

		# allow dropping of unwanted peptides
		if ( lowComplexityFilter) {
			cat( "  lowComplexity..")
			toDrop <- lowComplexityPeptides( ansDF$Peptide, ...)
			if ( length( toDrop)) {
				ansDF <- ansDF[ -toDrop, ]
				nout <- nrow(ansDF)
			}
		}
		if ( ! nout) next

		ntotal <- ntotal + nout
		cat( "  Peptides: ", prettyNum( as.integer( ntotal)), "  ", as.percent( ntotal, big.value=nUtotal*100, 
				percentSign=FALSE), "per Read")

		#write.table( ansDF, file=conOut, col.names=(ntotal == nout), sep="\t", quote=F, row.names=F)
		if ( ntotal == nout) {
			write.table( ansDF, file=fileout, append=FALSE, col.names=TRUE, sep="\t", quote=F, row.names=F)
		} else {
			write.table( ansDF, file=fileout, append=TRUE, col.names=FALSE, sep="\t", quote=F, row.names=F)
		}

		# perhaps stop early...
		if ( ! is.null( maxReads)) {
			if ( nRead >= maxReads) break
		}
		if ( ! is.null( maxPeptides)) {
			if ( ntotal >= maxPeptides) break
		}
	}

	close( conIn)
	#close( conOut)
	cat( "\nTotal_Reads:            ", prettyNum( as.integer(nread), big.mark=","))
	cat( "\nTotal_Peptides:         ", prettyNum( as.integer(ntotal), big.mark=","))
	cat( "\n")
	return( ntotal)
}


`fastqReader` <- function() {

	filename <- ""
	con <- NULL
	N <- 0

	initialize <- function( file) {

		cat( "\nInitializing FASTQ reader for: ", file)
		if ( ! file.exists( file)) {
			cat( "\nFile not found: ", file)
			return( "Error")
		}
		con <<- gzfile( file, open="rt")
		filename <<- file
		return( file)
	}

	read <- function( n=1) {
	
		nlines <- n * 4
		txt <- readLines( con, n=nlines)
		nread <- length( txt)
		if ( nread < 4) return(NULL)
		isRID <- seq( 1, nread, by=4)
		N <<- N + length(isRID)
		return( list( "rid"=base::sub( "^@","",txt[isRID]), "seq"=txt[isRID+1], "score"=txt[isRID+3]))
	}

	finalize <- function() {
		if ( !is.null(con)) base::close(con)
		cat( "\nRead", N, "records from: ", filename, "\n")
		return()
	}

	return( environment())
}


`fastqWriter` <- function() {

	filename <- ""
	con <- NULL
	N <- 0

	initialize <- function( file) {

		cat( "\nInitializing FASTQ writer for: ", file)
		con <<- gzfile( file, open="wt")
		filename <<- file
		return( file)
	}

	write <- function( rid, seq, score) {
	
		txt <- paste( "@", rid, "\n", seq, "\n+\n", score, sep="")
		writeLines( txt, con=con)
		N <<- N + 1
		return()
	}

	finalize <- function() {
		if ( !is.null(con)) base::close(con)
		cat( "\nWrote", N, "records to: ", filename, "\n")
		return()
	}

	return( environment())
}


`fastqToPairedFiles` <- function( file, secondFile=NULL, max.buf=20000) {

	options( warn=-1)
	on.exit( options( warn=0))

	fin <- fastqReader()
	if ( fin$initialize(file) != file) return()

	f1name <- sub( "fastq.gz", "1.fastq.gz", file)
	fout1 <- fastqWriter()
	fout1$initialize(f1name)
	on.exit( fout1$finalize(), add=TRUE)

	f2name <- sub( "fastq.gz", "2.fastq.gz", file)
	fout2 <- fastqWriter()
	fout2$initialize(f2name)
	on.exit( fout2$finalize(), add=TRUE)

	rids1 <- rids2 <- seqs1 <- seqs2 <- scores1 <- scores2 <- rep.int("", 100)
	n1 <- n2 <- nout <- 0

	find1 <- function( rid) {
			if ( n1 < 1) return(0)
			who <- (1:n1)[ rid == rids1[1:n1]][1]
			return( if (is.na(who)) 0 else who)
	}

	find2 <- function( rid) {
			if ( n2 < 1) return(0)
			who <- (1:n2)[ rid == rids2[1:n2]][1]
			return( if (is.na(who)) 0 else who)
	}

	store1 <- function( rid, seq, score) {
			where <- if ( n1 > 0) (1:n1)[ rids1[1:n1] == ""][1] else NA
			#where <- if ( n1 > 0) which( rids1[1:n1] == "")[1] else NA
			if ( is.na( where)) {
				n1 <<- n1 + 1
				where <- n1
			}
			rids1[where] <<- rid
			seqs1[where] <<- seq
			scores1[where] <<- score
			return()
	}

	store2 <- function( rid, seq, score) {
			where <- if ( n2 > 0) (1:n2)[ rids2[1:n2] == ""][1] else NA
			#where <- if ( n2 > 0) which( rids2[1:n2] == "")[1] else NA
			if ( is.na( where)) {
				n2 <<- n2 + 1
				where <- n2
			}
			rids2[where] <<- rid
			seqs2[where] <<- seq
			scores2[where] <<- score
			return()
	}

	squeeze1 <- function() {
			isBlank <- which( rids1 == "")
			if ( length( isBlank) > 0) {
				rids1 <<- rids1[ -isBlank]
				seqs1 <<- seqs1[ -isBlank]
				scores1 <<- scores1[ -isBlank]
				n1 <<- length(rids1)
				if ( n1 > max.buf) {
					drops <- 1 : (n1-max.buf)
					rids1 <<- rids1[ -drops]
					seqs1 <<- seqs1[ -drops]
					scores1 <<- scores1[ -drops]
					n1 <<- length(rids1)
				}
			}
	}

	squeeze2 <- function() {
			isBlank <- which( rids2 == "")
			if ( length( isBlank) > 0) {
				rids2 <<- rids2[ -isBlank]
				seqs2 <<- seqs2[ -isBlank]
				scores2 <<- scores2[ -isBlank]
				n2 <<- length(rids2)
				if ( n2 > max.buf) {
					drops <- 1 : (n2-max.buf)
					rids2 <<- rids2[ -drops]
					seqs2 <<- seqs2[ -drops]
					scores2 <<- scores2[ -drops]
					n2 <<- length(rids2)
				}
			}
	}



	# the main loop is to look for mate pairs within the file, buffering up as we go
	# so read one at a time
	cat("\n")
	repeat {
		item <- fin$read(1)
		if (is.null(item)) break
		
		rid1 <- rid2 <- item[[1]]
		nc <- nchar( rid1)
		mate <- substr( rid1, nc, nc)
		rid <- substr( rid1, 1, nc-1)
		seq <- item[[2]]
		score <- item[[3]]
		is1 <- (mate == '1')
		hit <- FALSE

		if (is1) {
			where2 <- find2( rid)
			if ( where2 == 0) {
				store1( rid, seq, score)
			} else {
				fout1$write( rid1, seq, score)
				rid2 <- paste( rid, "2", sep="")
				fout2$write( rid2, seqs2[where2], scores2[where2])
				rids2[where2] <- ""
				nout <- nout + 1
				hit <- TRUE
			}
		} else {
			where1 <- find1( rid)
			if ( where1 == 0) {
				store2( rid, seq, score)
			} else {
				fout2$write( rid2, seq, score)
				rid1 <- paste( rid, "1", sep="")
				fout1$write( rid1, seqs1[where1], scores1[where1])
				rids1[where1] <- ""
				nout <- nout + 1
				hit <- TRUE
			}
		}
		if ( hit && nout %% 100 == 0) {
			squeeze1()
			squeeze2()
			cat( "\rPairs:", nout, "  Buf1:", n1, "  Buf2:", n2)
		}
	}

	# done reading the file...
	fin$finalize()
	cat( "\rPairs:", nout, "  Buf1:", n1, "  Buf2:", n2)

	# if a second file was given, (the 'No Hits'), use it as a read only source of possible second mates
	if ( ! is.null( secondFile)) {
		if ( fin$initialize(secondFile) != secondFile) break
		cat("\n")
		repeat {
			item <- fin$read( n=10000)
			if (is.null(item)) break
			
			ridset <- item[[1]]
			seqset <- item[[2]]
			scoreset <- item[[3]]
			hit <- FALSE

			# see if any of whats in our buffers match these
			try1 <- paste( rids1, "2", sep="")
			where <- match( ridset, try1, nomatch=0)
			if ( any( where > 0)) {
				set2 <- which( where > 0)
				set1 <- where[ set2]
				for (j in 1:length(set2)) {
					j1 <- set1[j]
					j2 <- set2[j]
					fout1$write( rids1[j1], seqs1[j1], scores1[j1])
					fout2$write( ridset[j2], seqset[j2], scoreset[j2])
					rids1[j1] <- ""
					nout <- nout + 1
				}
				hit <- TRUE
			}

			try2 <- paste( rids2, "1", sep="")
			where <- match( ridset, try2, nomatch=0)
			if ( any( where > 0)) {
				set1 <- which( where > 0)
				set2 <- where[ set1]
				for (j in 1:length(set1)) {
					j1 <- set1[j]
					j2 <- set2[j]
					fout1$write( ridset[j1], seqset[j1], scoreset[j1])
					fout2$write( rids2[j2], seqs2[j2], scores2[j2])
					rids2[j2] <- ""
					nout <- nout + 1
				}
				hit <- TRUE
			}
			if (hit) {
				squeeze1()
				squeeze2()
				cat( "\rPairs:", nout, "  Buf1:", n1, "  Buf2:", n2, "   ")
				if ( n1 == 0 && n2 == 0) break
			}
		}
		fin$finalize()
	}

	# all done now
}


`fastqCompressFolder` <- function( folder, recursive=FALSE) {

	fastqPattern <- "(fq|fastq)$"
	fastqFiles <- dir( folder, recursive=recursive, pattern=fastqPattern, full.name=T)
	N <- length( fastqFiles)
	if ( ! N) return( NULL)

	out <- data.frame()
	cat( "\nNumber of .FASTQ files to compress:  ", length(fastqFiles), "\n")
	for ( i in 1:N) {
		f <- fastqFiles[i]
		cat( "\r",i, " Compressing ", basename(f))
		ans <- fastqCompressFile( f)
		if ( is.null(ans)) next
		out <- rbind( out, as.data.frame(ans))
		cat( "  done.  Pct_Reduction: ", ans$PctReduced)
	}

	rownames(out) <- 1:nrow(out)
	ttlReduction <- sum( out$Reduction)
	cat( "\nFile Size Reduction: ", prettyNum( ttlReduction, big.mark=","), "bytes.\n")
	return( out)
}


`fastqCompressFile` <- function( f) {

	info1 <- file.info( f)
	if ( info1$isdir) {
		cat( "\nSkipping folder: ", f, "\n")
		return( NULL)
	}
	f2 <- paste( f, "gz", sep=".")
	cmdline <- paste( "gzip", f, sep=" ")
	system( cmdline)
	info2 <- file.info( f2)

	size1 <- info1$size
	size2 <- info2$size
	reduction <- size1 - size2
	pct <- round( reduction * 100 / size1, digits=2)
	out <- data.frame( "File"=basename(f), "Original"=size1, "Compressed"=size2, 
			"Reduction"=reduction, "PctReduced"=pct, stringsAsFactors=F)
	return(out)
}
