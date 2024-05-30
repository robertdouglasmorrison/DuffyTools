# expressionFileSetToMatrix.R


`expressionFileSetToMatrix` <- function( fnames, fids, geneColumn=c( "GENE_ID", "GeneID"), 
					intensityColumn=c("INTENSITY","RPKM_M", "RANK"),
					missingGenes=c("na", "drop", "fill"), sep="\t",
					keepIntergenics=FALSE, shortenGeneNames=NULL, verbose=FALSE) {

	missingGenes <- match.arg( missingGenes)

	# make sure we can read those files
	filesOK <- file.exists( fnames)
	if ( !all( filesOK)) {
		cat( "\nSome expression files not found:\n")
		print( fnames[ !filesOK])
		return(NULL)
	}
	nFiles <- length( fnames)
	if ( length(fids) != nFiles) stop( "length(fids) must equal length(fnames)")

	# load each file in turn
	allData <- vector( mode="list")
	for( i in 1:nFiles) {
		tmp <- read.delim( fnames[i], as.is=T, sep=sep)
		if ( all( intensityColumn == "RANK") && (! ("RANK" %in% colnames(tmp)))) {
			tmp$RANK <- 1:nrow(tmp)
		}
		# allow integers as well as characters
		if ( is.character(intensityColumn)) {
			intenC <- base::match( intensityColumn, colnames(tmp), nomatch=0)
			if ( any( intenC > 0)) intenC <- intenC[ intenC > 0][1]
		} else {
			intenC <- as.integer(intensityColumn)
		}
		if ( is.character(geneColumn)) {
			geneC <- base::match( geneColumn, colnames(tmp), nomatch=0)
			if ( any( geneC > 0)) geneC <- geneC[ geneC > 0][1]
		} else{
			geneC <- as.integer(geneColumn)
		}
		if ( any( c(intenC,geneC) == 0)) {
			cat( "\nSome needed columns not found:   file: ", fnames[i],
				"\n  Expected: ", geneColumn, intensityColumn,
				"\n  Found:    ", colnames(tmp))
			return(NULL)
		}
		thisGenes <- tmp[[ geneC]]
		if ( ! is.null( shortenGeneNames)) {
			thisGenes <- shortGeneName( thisGenes, keep=shortenGeneNames)
		}
		thisInten <- tmp[[ intenC]]
		smallDF <- data.frame( "GENE_ID"=thisGenes, "INTENSITY"=thisInten, stringsAsFactors=F)
		if ( ! keepIntergenics) {
			drops <- grep ( "(ng)", thisGenes, fixed=T)
			if ( length(drops) > 0) smallDF <- smallDF[ -drops, ]
			nonGenes <- subset( getCurrentGeneMap(), REAL_G == FALSE)$GENE_ID
			drops <- which( smallDF$GENE_ID %in% nonGenes)
			if ( length(drops) > 0) smallDF <- smallDF[ -drops, ]
		}
		allData[[i]] <- smallDF
		if ( i == 1) {
			allGenes <- smallDF$GENE_ID
		} else {
			if ( missingGenes == "drop") {
				nWas <- length(allGenes)
				allGenes <- base::intersect( allGenes, smallDF$GENE_ID)
				nNow <- length(allGenes)
				if ( nNow < nWas) cat( "\nSome geneIDs missing: ", basename(fnames[i]))
			} else {
				allGenes <- base::union( allGenes, smallDF$GENE_ID)
			}
		}
		if ( verbose) cat( "\nFile: ", i, basename(fnames[i]), "\tN_Genes: ", nrow(smallDF))
	}
	if (verbose) cat( "\n")

	allGenes <- base::sort( setdiff( allGenes, ""))
	nGenes <- length( allGenes)
	m <- matrix( NA, nrow=nGenes, ncol=nFiles)
	colnames(m) <- fids
	rownames(m) <- allGenes

	# now fill in the matrix
	smallestV <- max( thisInten, na.rm=T)
	for( i in 1:nFiles) {
		v <- rep( NA, times=nGenes)
		smallDF <- allData[[i]]
		thisGenes <- smallDF$GENE_ID
		thisInten <- smallDF$INTENSITY
		where <- base::match( allGenes, thisGenes, nomatch=0)
		v[ where > 0] <- thisInten[ where]
		m[ , i] <- v
		if ( ! any( intensityColumn == "RANK")) smallestV <- min( smallestV, thisInten, na.rm=T)
	}

	# if we kept all genes, fill in the holes with a small value
	if (missingGenes == "fill") {
		m[ is.na(m)] <- smallestV
	}

	return( m)
}


`expressionMatrixToFileSet` <- function( m, groups=colnames(m), geneColumn="GENE_ID", 
					intensityColumn=c("INTENSITY","RPKM_M","READS_M", "TPM_M"),
					path=".", sep="\t", AVG.FUN=sqrtmean, addNGScolumns=FALSE,
					dropDuplicateIDs=TRUE, verbose=FALSE) {

	intensityColumn <- match.arg( intensityColumn)
	prefix <- getCurrentSpeciesFilePrefix()
	suffix <- if ( sep == ",") "csv" else "txt"
	gmap <- getCurrentGeneMap()

	# make sure we can write those files
	if ( ! file.exists( path)) dir.create( path, recursive=T)

	NC <- ncol(m)
	NG <- nrow(m)
	genes <- rownames(m)
	sids <- colnames(m)
	prods <- gene2Product(genes)

	# see if we reduce replicates by group
	if ( ! all( colnames(m) == groups)) {
		groupFac <- factor( groups)
		groupIDs <- levels( groupFac)
		NC <- nlevels(groupFac)
		cat( "\nReducing", ncol(m), "samples down to", NC, "groups..")
		m2 <- matrix( NA, NG, NC)
		colnames(m2) <- groupIDs
		rownames(m2) <- genes
		for ( i in 1:NG) m2[ i, ] <- tapply( m[i,], groupFac, FUN=AVG.FUN)
		m <- m2
		sids <- groupIDs
	}

	# note that TPM units by definition should sum to 1 million.  Averaging can distort that,
	# so restore that by hand
	if ( grepl( "TPM", toupper(intensityColumn))) {
		for ( i in 1:NC) {
			v <- m[ , i]
			m[ , i] <- round( v * 1000000 / sum(v,na.rm=T), digits=4)
		}
	}

	for ( i in 1:NC) {
		v <- m[,i]
		v <- round( v, digits=4)
		sml <- data.frame( "GeneID"=genes, "Product"=prods, "Value"=v, stringsAsFactors=F)
		ord <- order( v, decreasing=T)
		sml <- sml[ ord, ]
		rownames(sml) <- 1:NG
		colnames(sml) <- c( geneColumn, "PRODUCT", intensityColumn)

		# if we need to make this look like NGS data, add a few more data fields
		if (addNGScolumns) {
			if ( intensityColumn != "READS_M") sml$READS_M <- sml[[3]]
			if ( intensityColumn != "RPKM_M") sml$RPKM_M <- sml[[3]]
			sml$SIGMA_M <- rosettaSigma( sml[[3]])
			sml$N_EXON_BASES <- 1000
			wh1 <- match( sml[[1]], gmap$GENE_ID, nomatch=0)
			sml$N_EXON_BASES[ wh1 > 0] <- gmap$N_EXON_BASES[wh1]
			wh2 <- match( sml[[1]], gmap$NAME, nomatch=0)
			sml$N_EXON_BASES[ wh2 > 0] <- gmap$N_EXON_BASES[wh2]
		}

		# many downstream tools expect fully unique gene IDs
		if ( dropDuplicateIDs) {
			drops <- which( duplicated( sml[[1]]))
			if ( length( drops)) sml <- sml[ -drops, ]
		}

		outfile <- file.path( path, paste( sids[i], prefix, "Transcript", suffix, sep="."))
		write.table( sml, outfile, sep=sep, quote=(suffix == "csv"), row.names=F)
		if (verbose) cat( "\r", i, sids[i])
	}
	return( NC)
}
