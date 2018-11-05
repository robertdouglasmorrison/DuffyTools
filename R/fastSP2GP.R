# fastSP2GP.R _ next try at speedier 'SeqPos to GenePos' lookup...


`fastSP2GPsetup` <- function() {

	sp2gpData <- vector( mode="list")
	sp2gpSpecies <- sp2gpSeqID <- vector()
	nout <- 0

	curSpecies <- getCurrentSpecies()
	allSpecies <- getCurrentTargetSpecies()

	# build a 'findInterval' lookup table for every possible location
	for( spec in allSpecies) {
		setCurrentSpecies( spec)
		smap <- getCurrentSeqMap()
		gmap <- getCurrentGeneMap()
		emap <- getCurrentExonMap()

		seqIDs <- smap$SEQ_ID
		for( seqid in seqIDs) {

			rows <- which( gmap$SEQ_ID == seqid)
			if ( length(rows) < 1) next
		
			allGeneStarts <- gmap$POSITION[ rows]
	
			# some gene maps may have fake anti-sense genes
			dups <- which( duplicated( allGeneStarts) & gmap$REAL_G[rows] == FALSE)
			if ( length( dups) > 0) {
				rows <- rows[ -dups]
				allGeneStarts <- gmap$POSITION[ rows]
			}
			allGenePtrs <- rows
			allGeneNames <- gmap$GENE_ID[rows]
	
			# the Human genome has lots of genes on both strands at the same location...
			# try using Exons too, to get  a closer notion of the best geneID to use
			erows <- which( emap$SEQ_ID == seqid)
			if ( length( erows) > 0) {
				exonNames <- emap$GENE_ID[erows]
				exonStarts <- emap$POSITION[erows]
				exonPtrs <- match( exonNames, gmap$GENE_ID, nomatch=0)
				allGeneStarts <- c( allGeneStarts, exonStarts)
				allGenePtrs <- c( allGenePtrs, exonPtrs)
				allGeneNames <- c( allGeneNames, exonNames)
			}

			allGenePositions <- allGeneStarts
			allGenePositions[ allGenePtrs > 0] <- gmap$POSITION[ allGenePtrs]

			# add this seqID's lookup details to the set
			nout <- nout + 1
			sp2gpSeqID[ nout] <- seqid
			sp2gpSpecies[ nout] <- spec
			ord <- base::order( allGeneStarts)
			sml <- data.frame( "START"=allGeneStarts[ord], "NAME"=allGeneNames[ord], "PTR"=allGenePtrs[ord],
						"POSITION"=allGenePositions[ord], stringsAsFactors=FALSE)
			sp2gpData[[ nout]] <- sml
		}
	}
	fsp2gpData <<- sp2gpData
	fsp2gpSpecies <<- sp2gpSpecies
	fsp2gpSeqID <<- sp2gpSeqID

	setCurrentSpecies( curSpecies)

	return()
}


`fastSP2GPready` <- function() {

	if ( ! exists( "fsp2gpData")) return( FALSE)
	if ( ! all( getCurrentTargetSpecies() %in% fsp2gpSpecies)) return( FALSE)
	return( TRUE)
}


`fastSP2GPcleanup` <- function() {

	if ( exists( "fsp2gpData")) rm( fsp2gpData, envir=.GlobalEnv)
	if ( exists( "fsp2gpSpecies")) rm( fsp2gpSpecies, envir=.GlobalEnv)
	if ( exists( "fsp2gpSeqID")) rm( fsp2gpSeqID, envir=.GlobalEnv)
}


# turn chromosomal positions into gene positions
`fastSP2GP` <- function( seqid, seqbase) {

	if ( ! fastSP2GPready()) fastSP2GPsetup()

	Nbase <- length( seqbase)
	tapplyAns <- vector( mode="list")

	# set up storage to get filled
	geneOut <- speciesOut <- rep( NA, length=Nbase)
	gptrOut <- genePos <- rep( 0, length=Nbase)



	# local function for one SeqID at a time...
	`myFast_OneSeqFunc` <- function(x) {

		# given all the pointers to <seqid, base> tuples shared by one seqid
		thisSeq <- seqid[x[1]]
		who <- match( thisSeq, fsp2gpSeqID, nomatch=0)

		# if we can't find it...
		if ( who == 0) return()

		# OK for this seqID
		sp2gpData <- fsp2gpData[[ who]]
		thisSpecies <- fsp2gpSpecies[who]
		baseWhere <- fastFindInterval( seqbase[x], sp2gpData$START)
		baseWhere[ baseWhere == 0] <- 1
	
		# hand those names. positions back
		geneOut[ x] <<- sp2gpData$NAME[ baseWhere]
		speciesOut[ x] <<- thisSpecies
		gptrOut[ x] <<- sp2gpData$PTR[ baseWhere]
		genePos[ x] <<- seqbase[x] - sp2gpData$POSITION[ baseWhere] + 1
		return()
	}


	# factor and do it, one seqID at a time
	seqFac <- factor( seqid)
	tapply( 1:Nbase, INDEX=seqFac, FUN=myFast_OneSeqFunc)

	# send back the 'filled in' answer
	out <- list( "GENE_ID"=geneOut, "GENE_POSITION"=genePos, "SPECIES"=speciesOut, 
			"GMAP_PTR"=gptrOut)
	return( out)
}
