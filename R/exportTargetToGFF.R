# exprotTargetToGFF.R -- turn a target set of species into a single GFF file
#
#	as of July 2026, updating to use 'product=' attribute, on both mRNA and CDS lines
# 	(to better match GenBank submission requirements)

`exportTargetToGFF` <- function( target=getCurrentTarget()$TargetID, outfile=paste( target, "gff", sep="."), 
				genomicDNAfilePath=NULL) {

	# load the wanted target
	prevTarget <- curTarget <- getCurrentTarget()$TargetID
	prevSpecies <- getCurrentSpecies()
	if ( curTarget != target) curTarget <- setCurrentTarget( target)
	curSpeciesSet <- getCurrentTargetSpecies()

	if ( is.null(genomicDNAfilePath)) {
		cat( "\nNo checking for gene/pseudogene distinctions..")
		doAA <- FALSE
	} else {
		cat( "\nWill translate genes into proteins to check for pseudogenes..")
		doAA <- TRUE
	}

	# open the new file for writing, and put out the GFF header
	con <- file( outfile, open="wt")
	writeLines( "##gff-version 3", con=con)

	# put out each chromosome as a region pragma
	for ( s in curSpeciesSet) {
		setCurrentSpecies(s)
		smap <- getCurrentSeqMap()
		txt <- paste( "##sequence-region", smap$SEQ_ID, 1, smap$LENGTH, sep="\t")
		writeLines( txt, con=con)
	}


	# local functions
	write.Gene.Line <- function( seqid, pos, end, strand, attribs) {
		
		txt <- paste( seqid, ".", "gene", pos, end, ".", strand, ".", attribs, sep="\t")
		writeLines( txt, con=con)
	}


	write.mRNA.Line <- function( seqid, gid, pos, end, nam, strand, prod) {
		
		mRNAid <- paste( gid, "mRNA", sep=".")
		attribs <- paste( "ID=", mRNAid, ";Parent=", gid, ";Name=", nam, 
					";product=\"", prod, "\"", sep="")
		txt <- paste( seqid, ".", "mRNA", pos, end, ".", strand,
				".", attribs, sep="\t")
		writeLines( txt, con=con)
	}


	write.Exon.Lines <- function( seqid, gid, emap, cdsmap, prod) {
		
		Nexons <- nrow( emap)
		mRNAid <- paste( gid, "mRNA", sep=".")
		exonid <- paste( gid, "exon", 1:Nexons, sep=".")
		attribs <- paste( "ID=", exonid, ";Parent=", mRNAid, sep="")
		myStrand <- emap$STRAND[1]
		if ( Nexons > 1) {
			if (myStrand == "+") {
				ord <- order( emap$POSITION, decreasing=FALSE)
			} else {
				ord <- order( emap$END, decreasing=TRUE)
			}
			if ( any( ord != 1:Nexons)) emap <- emap[ ord, ]
		}
		pos <- as.integer( emap$POSITION)
		end <- as.integer( emap$END)
		txt <- paste( seqid, ".", "exon", pos, end, ".", emap$STRAND, ".", attribs, sep="\t")
		writeLines( txt, con=con)
		
		# calc the phasing for CDS
		Ncds <- nrow( cdsmap)
		if ( ! Ncds) return()
		cdsid <- paste( gid, "cds", 1:Ncds, sep=".")
		myStrand <- cdsmap$STRAND[1]
		if ( Ncds > 1) {
			if (myStrand == "+") {
				ord <- order( cdsmap$POSITION, decreasing=FALSE)
			} else {
				ord <- order( cdsmap$END, decreasing=TRUE)
			}
			if ( any( ord != 1:Ncds)) cdsmap <- cdsmap[ ord, ]
		}
		pos <- as.integer( cdsmap$POSITION)
		end <- as.integer( cdsmap$END)
		nbp <- end - pos + 1
		myPhase <- rep.int( 0, Ncds)
		if ( Ncds > 1) {
		    # perhaps the phasing is done backward when on reverse strand
			cumbp <- cumsum( nbp)
			isMod1 <- which( cumbp[1:(Ncds-1)] %% 3 == 1)
			myPhase[ isMod1 + 1] <- 2
			isMod2 <- which( cumbp[1:(Ncds-1)] %% 3 == 2)
			myPhase[ isMod2 + 1] <- 1
		}
		attribs <- paste( "ID=", cdsid, ";Parent=", mRNAid, 
					";product=\"", prod, "\"", sep="")
		txt <- paste( seqid, ".", "CDS", pos, end, ".", cdsmap$STRAND, myPhase, attribs, sep="\t")
		writeLines( txt, con=con)
	}


	# now visit every gene
	for ( s in curSpeciesSet) {
		setCurrentSpecies(s)
		smap <- getCurrentSeqMap()
		geneMap <- getCurrentGeneMap()
		exonMap <- getCurrentExonMap()
		cdsMap <- getCurrentCdsMap()
		for ( seqid in smap$SEQ_ID) {
			gmap <- subset.data.frame( geneMap, SEQ_ID == seqid & REAL_G == TRUE)
			if ( nrow( gmap) < 1) next
			emap <- subset.data.frame( exonMap, SEQ_ID == seqid)
			cmap <- subset.data.frame( cdsMap, SEQ_ID == seqid)
			ord <- order( gmap$POSITION)
			gmap <- gmap[ ord, ]
			cat( "\n", s, seqid, "\n")
			allGenes <- gmap$GENE_ID
			allNames <- gmap$NAME
			allPos <- as.integer( gmap$POSITION)
			allEnd <- as.integer( gmap$END)
			allStrands <- gmap$STRAND
			# clean the product strings a bit. Remove any attribute separators, and convert any commas and internal equal signs
			allProds <- gsub( ";", "", gmap$PRODUCT)
			allProds <- gsub( ",", "%2C", allProds)
			allProds <- gsub( "=", " ", allProds)
			allAttribs <- paste( "ID=", allGenes, ";", "Name=", allNames, sep="")
			for ( i in 1:nrow(gmap)) {
				thisG <- allGenes[i]
				if ( doAA) {
					# see we the protein is full length
					aaAns <- gene2Fasta( thisG, genomicDNAfilePath=genomicDNAfilePath, mode="aa")
					aa <- aaAns$seq[1]
					NAA <- nchar(aa)
					isPS <- TRUE
					startM <- ( substr(aa,1,1) == "M")
					endStop <- ( substr(aa,NAA,NAA) == "*")
					internStop <- grepl( "\\*", substr(aa, 1, (NAA-1)))
					if ( startM && endStop && ! internStop) isPS <- FALSE
					if ( isPS) {
						allAttribs[i] <- paste( allAttribs[i], "pseudo=true", sep=";")
					}
				}
				write.Gene.Line( seqid, allPos[i], allEnd[i], allStrands[i], allAttribs[i])
				eemap <- subset.data.frame( emap, GENE_ID == thisG)
				ccmap <- subset.data.frame( cmap, GENE_ID == thisG)
				if ( nrow( eemap)) {
					write.mRNA.Line( seqid, thisG, allPos[i], allEnd[i], allNames[i], 
							allStrands[i], prod=allProds[i])
					write.Exon.Lines( seqid, thisG, eemap, ccmap, prod=allProds[i])
				}
				if ( i %% 10 == 0) cat( "\r", i, thisG)
			}
		}
	}
	cat( "\nDone.\n")
	close( con)

	setCurrentTarget( prevTarget)
	setCurrentSpecies( prevSpecies)
}

