# pipe.GeneSetEnrichment.R -- investigate sets of DE genes like pathway and GO groups by enrichment

`pipe.GeneSetEnrichment` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  folderName="", 
				toolName=c( "MetaResults", "RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
				geneColumn=if(speciesID %in% MAMMAL_SPECIES) "GENE_NAME" else "GENE_ID", 
				groupColumn="Group", geneSets=defaultGeneSets(speciesID),
				descriptor="Enrichment", maxPvalue=0.05, wt.enrich=1, wt.pvalue=2,
				verbose=TRUE, ...) {

	toolName <- match.arg( toolName)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting pipe 'GeneSetEnrichment' on Sample Set: \n")
		print(sampleIDset)
		cat(  "\nUsing results from Species:  ", speciesID)
		cat(  "\nUsing results from DE tool:  ", toolName)
		if (is.character(geneSets)) cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}

	annT <- readAnnotationTable( annotationFile)
	allSamples <- unique( unlist( sampleIDset))
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allGroups <- unique( myAnnT[[ groupColumn]])
	Ngrps <- length( allGroups)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	nGenomeGenes <- nrow( subset.data.frame( getCurrentGeneMap(), REAL_G == TRUE))

	# make the paths to where we read gene lists and write results
	optT <- readOptionsTable( optionsFile)
	if (is.null( results.path)) results.path <- getOptionValue( optT, "results.path", notfound=".", verbose=F)
	dePath <- file.path( results.path, toolName, paste( prefix, folderName, sep="."))
	outPath <- file.path( dePath, descriptor)
	if ( ! file.exists( outPath)) dir.create( outPath, recursive=T)
	suffix <- c( "RR", "RP","SAM","EdgeR","DESeq","Meta")[ match( toolName, 
				c("RoundRobin","RankProduct","SAM","EdgeR","DESeq","MetaResults"))]

	for ( i in 1:Ngrps) {

		thisGroup <- allGroups[i]
		# get the DE file of gene ratios, that has the geneIDs in ranked order
		filein <- file.path( dePath, paste( thisGroup, prefix, suffix, "Ratio.txt", sep="."))
		# file name evolved for MetaResults
		if (suffix == "Meta") filein <- file.path( dePath, paste( thisGroup, prefix, suffix, "UP.txt", sep="."))
		if ( ! file.exists(filein) && suffix == "Meta") filein <- file.path( dePath, paste( thisGroup, prefix, suffix, "Ratio.txt", sep="."))
		tbl <- read.delim( filein, as.is=T)
		if ( ! (geneColumn %in% colnames(tbl))) {
			cat( "\nGene column not found.  Looked for: ", geneColumn, "\nFound: ", colnames(tbl))
			stop()
		}
		geneList <- tbl[[ geneColumn]]
		Ngenes <- length( geneList)

		# when given a full genome, do about 4-5 iterations, covering the first 10% of the genome.
		if ( Ngenes >= nGenomeGenes/2) {
			stopN <- round( (Ngenes*0.1) / 100) * 100
			startN <- round( (stopN*0.2) / 50) * 50
		} else {
			# given a smaller set of genes, use up to the first third
			stopN <- round( (Ngenes*0.333) / 10) * 10
			startN <- round( (stopN*0.2) / 10) * 10
		}
		stepN <- startN
		steps <- seq( startN, stopN, by=stepN)

		# do those enrichment calls
		ans <- geneSetMetaEnrichment( geneList, Ngenes=steps, geneSets=geneSets, wt.enrich=wt.enrich, 
					wt.pvalue=wt.pvalue, verbose=F)
		ans <- subset( ans, AVG_PVALUE <= maxPvalue)
		rownames(ans) <- 1:nrow(ans)

		# turn it to a HTML result
		outfile <- file.path( outPath, paste( thisGroup, "UP", "Enrichment.html", sep="."))
		otherGrps <- setdiff( allGroups, thisGroup)
		if ( length(otherGrps) > 1) otherGrps <- paste( "{", paste(otherGrps,collapse=" + "), "}", sep="")
		title <- paste( "Enrichment: &nbsp; Gene Sets most UP in &nbsp; ", thisGroup, " &nbsp; vs &nbsp; ", 
				otherGrps, " <br> Comparison Folder: &nbsp; ", folderName, " <br> Species: ", speciesID, sep="")
		metaEnrichment2html( ans, file=outfile, title=title)
		cat( "\nWrote file: ", outfile)

		# strip out the hyperlinks and make .CSV too
		ans$PathName <- cleanGeneSetModuleNames( ans$PathName, wrapParen=F)
		outfile <- file.path( outPath, paste( thisGroup, "UP", "Enrichment.csv", sep="."))
		write.table( ans, outfile, sep=",", quote=T, row.names=F)
		cat( "\nWrote file: ", outfile)

		# repeat for the down genes
		# in most cases, it's the exact same file, just invert it
		geneList <- rev( geneList)
		if (suffix == "Meta" && grepl( ".UP.",basename(filein), fixed=T)) {
			filein <- file.path( dePath, paste( thisGroup, prefix, suffix, "DOWN.txt", sep="."))
			tbl <- read.delim( filein, as.is=T)
			geneList <- tbl[[ geneColumn]]
		}
		ans <- geneSetMetaEnrichment( geneList, Ngenes=steps, geneSets=geneSets, wt.enrich=wt.enrich, 
					wt.pvalue=wt.pvalue, verbose=F)
		ans <- subset( ans, AVG_PVALUE <= maxPvalue)
		rownames(ans) <- 1:nrow(ans)
		outfile <- file.path( outPath, paste( thisGroup, "DOWN", "Enrichment.html", sep="."))
		title <- paste( "Pathway Enrichment: &nbsp; Meta Results: &nbsp; DOWN in '", thisGroup, "'", 
				" &nbsp;  Species: ", speciesID, sep="")
		metaEnrichment2html( ans, file=outfile, title=title)
		cat( "\nWrote file: ", outfile)
		ans$PathName <- cleanGeneSetModuleNames( ans$PathName, wrapParen=F)
		outfile <- file.path( outPath, paste( thisGroup, "DOWN", "Enrichment.csv", sep="."))
		write.table( ans, outfile, sep=",", quote=T, row.names=F)
		cat( "\nWrote file: ", outfile)
	}
	cat( "\nDone.\n")

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished pipe 'GeneSetEnrichment' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n")
	}

	return()
}
