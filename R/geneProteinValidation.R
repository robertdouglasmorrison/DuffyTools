# geneProteinValidation.R -- check if the annotation for a gene matches it's protein

`geneProteinValidation` <- function( gene, protein.sequence, 
					genomicFastaFile=getOptionValue("Options.txt","genomicFastaFile",verbose=F)) {

	# given one geneID, get all we know about its annotation
	gmap <- subset.data.frame( getCurrentGeneMap(), GENE_ID == gene[1])
	if ( nrow( gmap) != 1) {
		cat( "GeneID not unique in current gene map:  ", gene)
		return( NULL)
	}
	cmap <- subset.data.frame( getCurrentCdsMap(), GENE_ID == gene[1])
	if ( nrow( cmap) < 1) {
		cat( "GeneID has no CDS fragments:  ", gene)
		return( NULL)
	}
	if ( nrow( cmap) > 1) {
		ord <- cmap$POSITION
		cmap <- cmap[ ord, ]
	}
	nCDS <- nrow( cmap)
	nAA <- nchar( protein.sequence[1])
	mySeqID <- gmap$SEQ_ID
	myStrand <- gmap$STRAND
	needRevComp <- ( myStrand == "-")
	geneStart <- gmap$POSITION
	geneStop <- gmap$END

	# get the genomic DNA
	if ( is.list( genomicFastaFile)) {
		fa <- genomicFastaFile
	} else {
		fa <- loadFasta( genomicFastaFile)
	}
	# grab this gene's chunk from the genomic DNA, and translate it in all 6 frames
	who <- match( mySeqID, fa$desc)
	if ( is.na( who)) {
		cat( "Chromosome not found in genomic DNA.  Looked for: ", mySeqID)
		return( NULL)
	}
	geneDNA <- base::substr( fa$seq[who], geneStart, geneStop)
	if (needRevComp) geneDNA <- myReverseComplement( geneDNA)
	#geneGenomicAA <- DNAtoAA( geneDNA, clipAtStop=F, readingFrame=1:6)

	# build some mappings between CDNA and Protein
	if (needRevComp) {
		cmap$GENE_POSITION <- geneStop - cmap$END + 1
		cmap$GENE_END <- geneStop - cmap$POSITION + 1
	} else {
		cmap$GENE_POSITION <- cmap$POSITION - geneStart + 1
		cmap$GENE_END <- cmap$END - geneStart + 1
	}
	cmap$LEN_DNA <- abs( cmap$GENE_END - cmap$GENE_POSITION + 1)
	cmap$PROT_POSITION <- NA
	cmap$PROT_END <- NA
	cmap$LEN_AA <- round( cmap$LEN_DNA / 3, digits=2)
	# for the protein, we need to step along
	nAA <- 0
	if (needRevComp) {
		for( i in nCDS:1) {
			cmap$PROT_POSITION[i] <- nAA + 1
			cmap$PROT_END[i] <- floor( nAA + cmap$LEN_AA[i])
			nAA <- cmap$PROT_END[i]
		}
	} else {
		for( i in 1:nCDS) {
			cmap$PROT_POSITION[i] <- nAA + 1
			cmap$PROT_END[i] <- floor( nAA + cmap$LEN_AA[i])
			nAA <- cmap$PROT_END[i]
		}
	}

	# at this point, the CDS chunks of protein 'should' match the reference
	cmap$REF_PROT <- NA
	cmap$ANNO_PROT <- NA
	cmap$AA_EDITDIST <- NA
	for ( i in 1:nCDS) {
		cmap$REF_PROT[i] <- aaRef <- base::substr( protein.sequence, cmap$PROT_POSITION[i], cmap$PROT_END[i])
		cmap$ANNO_PROT[i] <- aaAnno <- DNAtoBestPeptide( substr( geneDNA, cmap$GENE_POSITION[i], 
								cmap$GENE_END[i]), tieBreakMode="reference", 
								reference=protein.sequence)
		cmap$AA_EDITDIST[i] <- adist( aaRef, aaAnno)[1]
	}

	# for now, just pass that cmap back
	return( cmap)
}


`validateAllProteins` <- function( proteinFastaFile, 
				genomicFastaFile=getOptionValue("Options.txt","genomicFastaFile",verbose=F)) {

	# get the full set of protein coding genes
	genes <- sort( unique( getCurrentCdsMap()$GENE_ID))
	nGenes <- length( genes)
	cat( "\nN_Protein Coding Genes: ", nGenes)

	# grab all the proteins and all the chromosomes
	if ( is.list( proteinFastaFile)) {
		protFA <- proteinFastaFile
	} else {
		protFA <- loadFasta( proteinFastaFile)
	}
	if ( is.list( genomicFastaFile)) {
		chrFA <- genomicFastaFile
	} else {
		chrFA <- loadFasta( genomicFastaFile)
	}

	# visit each one
	MATCH <- base::match
	nBadExons <- editDist <- rep.int( NA, nGenes)
	for ( i in 1:nGenes) {
		g <- genes[i]
		whereP <- MATCH( g, protFA$desc)
		if ( is.na( whereP)) next
		protSeq <- protFA$seq[whereP]
		ans <- geneProteinValidation( g, protSeq, genomicFastaFile=chrFA)
		if ( is.null(ans)) next
		nBadExons[i] <- sum( ans$EDITDIST > 0)
		editDist[i] <- sum( ans$EDITDIST, na.rm=T)
		cat( "\r", i, g, nBadExons[i], editDist[i])
	}

	out <- data.frame( "GENE_ID"=genes, "N_BAD_EXONS"=nBadExons, "PROT_EDITDIST"=editDist,
			stringsAsFactors=F)
	out
}
