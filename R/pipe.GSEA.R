# pipe.GSEA -- wrapper to the 'fgsea' gene set enrichment package

`pipe.GSEA` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, folderName=NULL,
				tool=c( "MetaResults", "DESeq", "EdgeR", "RankProduct", "RoundRobin", "SAM", "LimmaVoom"),
				m=NULL, contrast=NULL, useLog=TRUE,
				groupColumn="Group", geneSets=defaultGeneSets(speciesID), descriptor="CombinedGeneSets", 
				minGenesPerSet=if (speciesID %in% MAMMAL_SPECIES) 4 else 2, 
				min.rpkm=1, min.variance=1, offset=2,
				addCellTypes=TRUE, doGSEA= TRUE, verbose=TRUE, 
				# args specific to GSEA
				minSize=minGenesPerSet, maxSize=1000, eps=0, nPermSimple=5000, ...)
{

	require( fgsea)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# we need a subfolder for where to stuff these results
	if ( is.null( folderName)) stop( "Required argument 'folderName' not specified")

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting pipe 'GSEA' on Sample Set: \n")
		print(sampleIDset)
		cat(  "\nUsing results from Species:  ", speciesID)
		if ( is.character(geneSets)) cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}

	# get the samples and their groups
	sampleIDset <- unique( unlist( sampleIDset))
	annT <- readAnnotationTable( annotationFile)
	where <- match( sampleIDset, annT$SampleID, nomatch=0)
	myAnnT <- annT[ where, ]
	sids <- myAnnT$SampleID
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\nGroup column name not in Annotation table.")
		cat( "\nLooked for:  ", groupColumn)
		cat( "\nFound:  ", colnames(annT))
		stop()
	}
	# qusage has expectations about what the group string can contain
	grps <- myAnnT[[ groupColumn]]
	grps <- gsub( "-", "_", grps, fixed=T)
	grps <- gsub( " ", "_", grps, fixed=T)
	allGroups <- unique( grps)
	Ngroups <- length( allGroups)

	if ( is.null( results.path)) results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=F)

	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode="combined")
	bigSetOfGeneSets <- GSanswer$geneSets
	myGeneSets <- bigSetOfGeneSets[[1]]
	if ( is.character( geneSets) && length(geneSets) == 1) descriptor=geneSets[1]

	# force the check/drop of gene sets with too few genes in them, for this matrix of genes
	# this step is done by GSEA explicitly
	myGeneSetNames <- names(myGeneSets)
	myGeneSetShortNames <- names(myGeneSets) <- cleanGeneSetModuleNames( myGeneSetNames, wrapParentheses=FALSE)

	# ready to call GSEA
	cat( "\n")
	tool <- match.arg( tool)
	path <- file.path( results.path, tool, paste( prefix,folderName,sep="."), "GSEA")
	if ( ! file.exists( path)) dir.create( path, recursive=T)

		# do all groups, allow it to be multi-threaded
	multicore.gsea <- function( g) {
				group1 <- g
				group2 <- paste( "Not", g, sep="")
				myContrast <- paste( group1, group2, sep="-")
				cat( "\nContrast Grouping descriptor:  ", myContrast, "\n")
				do.GSEA( geneSets=myGeneSets, group1=group1, descriptor=descriptor, path=path, 
					prefix=prefix, makeDownHTML=(Ngroups>2), tool=tool, labels=allGroups, 
					shortNames=myGeneSetShortNames, longNames=myGeneSetNames, addCellTypes=addCellTypes,
					doGSEA=doGSEA, minSize=minSize, maxSize=maxSize, eps=eps, nPermSimple=nPermSimple, ...)
				return(NULL)
			}
	multicore.lapply( allGroups, multicore.gsea)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished pipe 'GSEA' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\nGroups: ", allGroups, "\n")
	}
	return()
}


do.GSEA <- function( geneSets, group1="Group1", descriptor="GeneSets",
			path=".", prefix=getCurrentSpeciesFilePrefix(), makeDownHTML=TRUE, 
			tool=c( "MetaResults", "DESeq", "EdgeR", "RankProduct", "RoundRobin", "SAM", "LimmaVoom"),
			labels="", shortNames=names(geneSets), longNames=names(geneSets), addCellTypes=FALSE, 
			doGSEA=TRUE, minSize=2, maxSize=1000, eps=0, nPermSimple=5000, ...) {

	# GSEA starts from a DE results file.  Make that filename and load that data table
	tool <- match.arg(tool)
	if (tool == "MetaResults") {
		deFilename <- paste( group1, prefix, "Meta.JOINED.txt", sep=".")
	} else if (tool == "RoundRobin") {
		deFilename <- paste( group1, prefix, "RR.Ratio.txt", sep=".")
	} else if (tool == "RankProduct") {
		deFilename <- paste( group1, prefix, "RP.Ratio.txt", sep=".")
	} else if (tool == "SAM") {
		deFilename <- paste( group1, prefix, "SAM.Ratio.txt", sep=".")
	} else if (tool == "DESeq") {
		deFilename <- paste( group1, prefix, "DESeq.Ratio.txt", sep=".")
	} else if (tool == "EdgeR") {
		deFilename <- paste( group1, prefix, "EdgeR.Ratio.txt", sep=".")
	} else if (tool == "LimmaVoom") {
		deFilename <- paste( group1, prefix, "LV.Ratio.txt", sep=".")
	} else { stop( paste( "Unknown tool for GSEA: ", tool))}
	# given the path for where to put the results, find the DE results up one level
	de.path <- dirname( path)
	deFilename <- file.path( de.path, deFilename)
	if ( ! file.exists( deFilename)) {
		cat( "\nError:  DE results file not found: ", deFilename)
		return(NULL)
	}
	deDF <- read.delim( deFilename, as.is=T)
	
	if (doGSEA) {

		# GSEA wants the vector of fold change data, with the genes for the names
		stats <- as.numeric( deDF$LOG2FOLD)
		names(stats) <- shortGeneName( deDF$GENE_ID, keep=1)
		
		# it seems GSEA can break for a few reasons:
		# 1) duplicate gene names: only reasonable fix is to drop duplicates, so first 
		# instance is retained
		drops <- which( duplicated( names(stats)))
		if ( length(drops)) stats <- stats[ -drops]
		# 2) duplicate numeric values:  only reasonable fix is to introduce random jitter to break ties
		while ( any( duplicated( stats))) stats <- jitter( stats)

		# the package seems to have some changing command line arguments, so call it differently as needed
		majorVersion <- as.numeric( R.version$major)

		# ok, call GSEA
		if (majorVersion >= 4) {
			SAVGSEA <<- ans <- suppressWarnings( fgsea( pathways=geneSets, stats=stats, minSize=minSize, maxSize=maxSize, 
						eps=eps, nPermSimple=nPermSimple, scoreType="pos", ...))
		} else {
			SAVGSEA <<- ans <- suppressWarnings( fgsea( pathways=geneSets, stats=stats, minSize=minSize, maxSize=maxSize, 
						nperm=nPermSimple, ...))
		}

		# first put the result into our wanted UP to DOWN order
		ans$NES <- as.numeric( ans$NES)
		ans$pval <- as.numeric( ans$pval)
		ord <- diffExpressRankOrder( ans$NES, ans$pval)
		ans <- ans[ ord, ]
		# now extract and clean up the wanted results
		pval <- ans$pval
		padj <- ans$padj
		nes <- round( ans$NES, digits=4)
		pathSize <- as.numeric( ans$size)
		pathName <- ans$pathway
		gNames <- ans$leadingEdge
		pathGenes <- sapply( gNames, function(x) paste( sort( unique( x)), collapse=";"))
		# create the pathway fold change from the input data
		pathFold <- sapply( gNames, function(x) {
							wh <- match(x,names(stats))
							myFC <- as.numeric( stats[wh])
							return( mean( myFC, na.rm=T))
					})
		pathFold <- round( pathFold, digits=3)
	
		# we gave GSEA the short names, put the longer full names back
		where <- match( pathName, shortNames)
		pathName <- longNames[ where]
		if ( addCellTypes) {
			cellType <- geneSetCellType( pathName, max.type=4)
			nKeep <- 6
			out <- data.frame( "PATHWAY"=pathName, "CellType"=cellType, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "PADJUST"=padj, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		} else {
			nKeep <- 5
			out <- data.frame( "PATHWAY"=pathName, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "PADJUST"=padj, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		}
		rownames(out) <- 1:nrow(out)
		outfile <- file.path( path, paste( group1, prefix, "UP.GSEA", "csv", sep="."))
		write.table( out, outfile, sep=",", quote=T, row.names=F)
	} else {
		cat( "\nUsing previously calculated GSEA results..")
		outfile <- file.path( path, paste( group1, prefix, "UP.GSEA", "csv", sep="."))
		out <- read.csv( outfile, as.is=T)
		nKeep <- if (addCellTypes) 6 else 5
	}

	otherGrpString <- paste( setdiff( sort( unique( labels)), group1), collapse=" + ")
	otherGrpString <- paste( "{", otherGrpString, "}", sep="")
	html1 <- subset( out, LOG2FOLD >= 0.1 & PVALUE <= 0.05)[ , 1:nKeep]
	if ( nrow(html1)) {
		htmlfile <- sub( "csv$", "html", outfile)
		htmltitle <- paste( "GSEA:  '", descriptor, "'  Pathways UP in: &nbsp; ", group1, 
					" &nbsp; vs &nbsp; ", otherGrpString, sep="")
		table2html( html1, htmlfile, title=addSpeciesToHtmlTitle(htmltitle), linkColumnNames=NULL)
	}
	# we can also make an enrichment table of cell types
	if (addCellTypes) {
		out1 <- subset( out, LOG2FOLD >= 0.1 & PVALUE <= 0.05)
		if ( nrow(out1)) {
			enrich <- cellTypeEnrichment( out1$CellType, mode="geneSets", upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( group1, prefix, "UP.GSEA.CellTypeEnrichment.csv", sep=".")
			f <- file.path( path, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
	}
	
	# for GSEA of down genes, run it all again
	if ( makeDownHTML ) {
	
		if (doGSEA) {
		
		# ok, call GSEA again
		if (majorVersion >= 4) {
			SAVGSEA <<- ans <- suppressWarnings( fgsea( pathways=geneSets, stats=stats, minSize=minSize, maxSize=maxSize, 
						eps=eps, nPermSimple=nPermSimple, scoreType="neg", ...))
		} else {
			SAVGSEA <<- ans <- suppressWarnings( fgsea( pathways=geneSets, stats=stats, minSize=minSize, maxSize=maxSize, 
						nperm=nPermSimple, ...))
		}

		# first put the result into our wanted UP to DOWN order
		ans$NES <- as.numeric( ans$NES)
		ans$pval <- as.numeric( ans$pval)
		ord <- diffExpressRankOrder( -ans$NES, ans$pval)
		ans <- ans[ ord, ]
		# now extract and clean up the wanted results
		pval <- ans$pval
		padj <- ans$padj
		nes <- round( ans$NES, digits=4)
		pathSize <- as.numeric( ans$size)
		pathName <- ans$pathway
		gNames <- ans$leadingEdge
		pathGenes <- sapply( gNames, function(x) paste( sort( unique( x)), collapse=";"))
		# create the pathway fold change from the input data
		pathFold <- sapply( gNames, function(x) {
							wh <- match(x,names(stats))
							myFC <- as.numeric( stats[wh])
							return( mean( myFC, na.rm=T))
					})
		pathFold <- round( pathFold, digits=3)
	
		# we gave GSEA the short names, put the longer full names back
		where <- match( pathName, shortNames)
		pathName <- longNames[ where]
		if ( addCellTypes) {
			cellType <- geneSetCellType( pathName, max.type=4)
			nKeep <- 6
			out <- data.frame( "PATHWAY"=pathName, "CellType"=cellType, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "PADJUST"=padj, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		} else {
			nKeep <- 5
			out <- data.frame( "PATHWAY"=pathName, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "PADJUST"=padj, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		}
		rownames(out) <- 1:nrow(out)
		outfile <- file.path( path, paste( group1, prefix, "DOWN.GSEA", "csv", sep="."))
		write.table( out, outfile, sep=",", quote=T, row.names=F)
	} else {
		cat( "\nUsing previously calculated GSEA results..")
		outfile <- file.path( path, paste( group1, prefix, "DOWN.GSEA", "csv", sep="."))
		out <- read.csv( outfile, as.is=T)
		nKeep <- if (addCellTypes) 6 else 5
	}

	otherGrpString <- paste( setdiff( sort( unique( labels)), group1), collapse=" + ")
	otherGrpString <- paste( "{", otherGrpString, "}", sep="")
	html2 <- subset( out, LOG2FOLD <= -0.1 & PVALUE <= 0.05)[ , 1:nKeep]
	if ( nrow(html2)) {
		htmlfile <- sub( "csv$", "html", outfile)
		htmltitle <- paste( "GSEA:  '", descriptor, "'  Pathways DOWN in: &nbsp; ", group1, 
					" &nbsp; vs &nbsp; ", otherGrpString, sep="")
		table2html( html2, htmlfile, title=addSpeciesToHtmlTitle(htmltitle), linkColumnNames=NULL)
	}
	# we can also make an enrichment table of cell types
	if (addCellTypes) {
		out2 <- subset( out, LOG2FOLD <= -0.1 & PVALUE <= 0.05)
		if ( nrow(out2)) {
			enrich <- cellTypeEnrichment( out2$CellType, mode="geneSets", upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( group1, prefix, "DOWN.GSEA.CellTypeEnrichment.csv", sep=".")
			f <- file.path( path, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
	}
	}  # end if MakeDOWN...
	
	return(out)
}

