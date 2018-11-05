# validateMapSet.R

# verify that the mapSet (extracted from a .gff file or any other sources) is internally consistent.


validateMapSet <- function( mapset, checkOverlaps=FALSE) {

# local functions...


validateSeqMap <- function( mydf, speciesID) {

	neededColumns <- c( "SEQ_ID", "LENGTH")
	anyBad <- FALSE
	if ( ! all( neededColumns %in% colnames( mydf))) stop( paste( "seqMap is missing required columns: ", neededColumns))
	if ( length( whoBad <- which( duplicated( mydf$SEQ_ID)))) {
		cat( "\nseqMap has non-unique SEQ_ID entries:\n")
		print( mydf[whoBad, ])
		anyBad <- TRUE
	}
	if ( any( is.na( as.integer( mydf$LENGTH)))) stop( "seqMap has non-integer LENGTH entries")
	if ( length( whoBad <- which( mydf$LENGTH < 2))) {
		cat( "\nseqMap has invalid LENGTH entries:\n")
		print( mydf[whoBad,])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( base::substr( mydf$SEQ_ID, 1, base::nchar(speciesID)) != speciesID))) {
		cat( "\nseqMap has SEQ_IDs that do not begin with 'speciesID':\n")
		print( mydf[whoBad,])
		anyBad <- TRUE
	}
	if (anyBad) stop( "Fix SeqMap issues before maps can be validated..")
	return( )
}


validateGeneMap <- function( mydf, seqMap, checkOverlaps=FALSE) {

	neededColumns <- c( "GENE_ID", "POSITION", "END", "SEQ_ID", "STRAND")
	anyBad <- FALSE
	if ( ! all( neededColumns %in% colnames( mydf))) 
		stop( paste( "geneMap is missing required columns: ", neededColumns))
	if ( length( whoBad <- which( duplicated( mydf$GENE_ID)))) {
		cat( "geneMap has non-unique GENE_ID entries:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( ! all( unique.default( mydf$SEQ_ID) %in% seqMap$SEQ_ID)) stop( "some geneMap SEQ_ID entries not in seqMap")
	if ( ! all( unique.default( mydf$STRAND) %in% c("+","-","",NA))) stop( "some geneMap STRAND entries not valid")
	if ( any( is.na( as.integer( mydf$POSITION)))) stop( "geneMap has non-integer POSITION entries")
	if ( any( is.na( as.integer( mydf$END)))) stop( "geneMap has non-integer END entries")
	if ( length( whoBad <- which( mydf$POSTION < 1))) {
		cat( "geneMap has POSITION entries less than 1:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END < 2))) {
		cat( "geneMap has END entries less than 2:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	lastBase <- getSeqLength( mydf$SEQ_ID, seqMap)
	if ( length( whoBad <- which( mydf$POSITION > lastBase))) {
		cat( "geneMap has POSITION entries past SEQ_ID's length:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END > lastBase))) {
		cat( "geneMap has END entries past SEQ_ID's length:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if (anyBad) stop( "Fix GeneMap issues before maps can be validated..")

	# force non-decreasing order within each chromosome
	checkNonDecreasingOrder( mydf)

	if ( checkOverlaps) checkGeneOverlaps( mydf)
	return( )
}


checkNonDecreasingOrder <- function( gmap) {

	seqFac <- factor( gmap$SEQ_ID)
	Good <- TRUE
	tapply( 1:nrow(gmap), seqFac, function(x) {
		pos <- gmap$POSITION[x]
		diffs <- diff( pos)
		if ( any( diffs < 0)) {
			cat( "\nGenes out of order: ", gmap$SEQ_ID[x[1]])
			cat( "\nRows: ", x[ which( diffs < 0)])
			Good <<- FALSE
		}
	})
	if (! Good) {
		cat("\n")
		stop( "Gene POSITIONs must be in non-decreasing order")
	} else {
		return(TRUE)
	}
}


checkGeneOverlaps <- function( gmap) {

	prevG <- 0
	prevSeq <- ""
	hasNonGenes <- ( ("REAL_G" %in% colnames( gmap)) && any( gmap$REAL_G == FALSE))

	# run along the genes and see if any overlap
	nOver <- 0
	for( i in 1:nrow( gmap)) {
		if ( hasNonGenes && gmap$REAL_G[i] != TRUE) next

		if ( gmap$SEQ_ID[i] != prevSeq) {
			prevSeq <- gmap$SEQ_ID[i]
			prevG <- i
			next
		}

		if ( gmap$POSITION[i] <= gmap$END[ prevG]) {
			if ( gmap$STRAND[i] == gmap$STRAND[prevG]) {
				cat( "\nOverlap:  ", gmap$GENE_ID[prevG], gmap$GENE_ID[i])
				cat( "\nExtents:  ", gmap$POSITION[prevG], "-",gmap$END[prevG], "  <->  ", gmap$POSITION[i], "-", gmap$END[i])
				cat( "\n")
				nOver <- nOver + 1
			}
		}

		prevG <- i
	}
	cat( "\n\nGene Overlap summary: \nN_Overlapping_Genes = ", nOver, "\n\n")
}


validateExonMap <- function( mydf, seqMap, geneMap, checkOverlaps=FALSE) {

	neededColumns <- c( "GENE_ID", "POSITION", "END", "SEQ_ID", "STRAND")
	anyBad <- FALSE
	if ( ! all( neededColumns %in% colnames( mydf))) stop( paste( "exonMap is missing required columns: ", neededColumns))
	if ( ! all( unique.default( mydf$SEQ_ID) %in% seqMap$SEQ_ID)) stop( "some exonMap SEQ_ID entries not in seqMap")
	if ( ! all( unique.default( mydf$GENE_ID) %in% geneMap$GENE_ID)) stop( "some exonMap GENE_ID entries not in geneMap")
	if ( any( is.na( as.integer( mydf$POSITION)))) stop( "exonMap has non-integer POSITION entries")
	if ( any( is.na( as.integer( mydf$END)))) stop( "exonMap has non-integer END entries")
	if ( length( whoBad <- which( mydf$POSTION < 1))) {
		cat( "exonMap has POSITION entries less than 1:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END < 2))) {
		cat( "exonMap has END entries less than 2:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	geneBases <- getGeneLimits( mydf$GENE_ID, geneMap)
	if ( length( whoBad <- which( mydf$POSITION < geneBases$POSITION))) {
		cat( "exonMap has POSITION entries before GENE_ID's start:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$POSITION > geneBases$END))) {
		cat( "exonMap has POSITION entries past GENE_ID's end:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END < geneBases$POSITION))) {
		cat( "exonMap has END entries before GENE_ID's start:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END > geneBases$END))) {
		cat( "exonMap has END entries past GENE_ID's end:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if (anyBad) stop( "Fix ExonMap issues before maps can be validated..")

	if ( checkOverlaps) checkExonOverlaps( mydf, geneMap)
	return( )
}


validateCdsMap <- function( mydf, seqMap, geneMap, checkOverlaps=FALSE) {

	neededColumns <- c( "GENE_ID", "POSITION", "END", "SEQ_ID", "STRAND")
	anyBad <- FALSE
	if ( ! all( neededColumns %in% colnames( mydf))) stop( paste( "cdsMap is missing required columns: ", neededColumns))
	if ( ! all( unique.default( mydf$SEQ_ID) %in% seqMap$SEQ_ID)) stop( "some cdsMap SEQ_ID entries not in seqMap")
	if ( ! all( unique.default( mydf$GENE_ID) %in% geneMap$GENE_ID)) stop( "some cdsMap GENE_ID entries not in geneMap")
	if ( any( is.na( as.integer( mydf$POSITION)))) stop( "cdsMap has non-integer POSITION entries")
	if ( any( is.na( as.integer( mydf$END)))) stop( "cdsMap has non-integer END entries")
	if ( length( whoBad <- which( mydf$POSTION < 1))) {
		cat( "cdsMap has POSITION entries less than 1:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END < 2))) {
		cat( "cdsMap has END entries less than 2:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	geneBases <- getGeneLimits( mydf$GENE_ID, geneMap)
	if ( length( whoBad <- which( mydf$POSITION < geneBases$POSITION))) {
		cat( "cdsMap has POSITION entries before GENE_ID's start:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$POSITION > geneBases$END))) {
		cat( "cdsMap has POSITION entries past GENE_ID's end:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END < geneBases$POSITION))) {
		cat( "cdsMap has END entries before GENE_ID's start:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if ( length( whoBad <- which( mydf$END > geneBases$END))) {
		cat( "cdsMap has END entries past GENE_ID's end:\n")
		print( mydf[ whoBad, ])
		anyBad <- TRUE
	}
	if (anyBad) stop( "Fix CdsMap issues before maps can be validated..")

	if ( checkOverlaps) checkExonOverlaps( mydf, geneMap)
	return( )
}


checkExonOverlaps <- function( emap, geneMap) {

	prevG <- 0
	prevSeq <- ""

	# run along the genes and see if any overlap
	nOver <- 0
	for( i in 1:nrow( emap)) {

		if ( emap$SEQ_ID[i] != prevSeq) {
			prevSeq <- emap$SEQ_ID[i]
			prevG <- i
			next
		}

		if ( emap$POSITION[i] <= emap$END[ prevG]) {
			if ( emap$STRAND[i] == emap$STRAND[prevG]) {
				cat( "\nOverlap:  ", emap$GENE_ID[prevG], emap$GENE_ID[i])
				cat( "\nExtents:  ", emap$POSITION[prevG], "-",emap$END[prevG], "  <->  ", emap$POSITION[i], "-", emap$END[i])
				cat( "\n")
				nOver <- nOver + 1
			}
		}

		# also verify that it lands "in" its gene...
		gptr <- base::match( emap$GENE_ID[i], geneMap$GENE_ID, nomatch=0)
		if ( gptr == 0) {
			cat( "\n\nMissing gene: ", emap$GENE_ID[i], "\n")
		} else {
			if ( emap$POSITION[i] < geneMap$POSITION[gptr]) {
				cat( "\nBoundary error:  (exon too early)", emap$GENE_ID[i], emap$POSITION[i], geneMap$POSITION[gptr])
			}
			if ( emap$END[i] > geneMap$END[gptr]) {
				cat( "\nBoundary error:  (exon too late)", emap$GENE_ID[i], emap$END[i], geneMap$END[gptr])
			}
		}

		prevG <- i
	}
	cat( "\n\nExon Overlap summary: \nN_Overlapping_Exons = ", nOver, "\n\n")
}


validateRRNA_Map <- function( mydf) {

	if ( is.null( mydf)) {
		cat( "\nRibosomal RNA map not found...  skipping.\n")
		return( NULL)
	}
	neededColumns <- c( "GENE_ID", "POSITION", "END", "SEQ_ID", "STRAND")
	if ( ! all( neededColumns %in% colnames( mydf))) 
		stop( paste( "rrnaMap is missing required columns: ", neededColumns))
	return( )
}


getSeqLength <- function( seqids, seqMap) {

	out <- rep( 0, times=length(seqids))
	where <- base::match( seqids, seqMap$SEQ_ID, nomatch=0)
	out[ where > 0] <- seqMap$LENGTH[ where]
	return( out)
}

getGeneLimits <- function( geneids, geneMap) {

	outB <- outE <- rep( 0, times=length(geneids))
	where <- base::match( geneids, geneMap$GENE_ID, nomatch=0)
	outB[ where > 0] <- geneMap$POSITION[ where]
	outE[ where > 0] <- geneMap$END[ where]
	return( data.frame( "POSITION"=outB, "END"=outE))
}

	# end of local functions...


	# check each map
	speciesID <- mapset$speciesID
	seqMap <- mapset$seqMap
	validateSeqMap( seqMap, speciesID)
	geneMap <- mapset$geneMap
	validateGeneMap( geneMap, seqMap, checkOverlaps)
	exonMap <- mapset$exonMap
	validateExonMap( exonMap, seqMap, geneMap, checkOverlaps)
	cdsMap <- mapset$cdsMap
	validateCdsMap( cdsMap, seqMap, geneMap, checkOverlaps)
	rrnaMap <- mapset$rrnaMap
	validateRRNA_Map( rrnaMap)
	return( "OK")
}

