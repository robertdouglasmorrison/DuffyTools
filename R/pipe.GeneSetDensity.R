# pipe.GeneSetDensity.R -- investigate sets of DE genes like pathway and GO groups, as Density curves

`pipe.GeneSetDensity` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  folderName="", 
				toolName=c( "MetaResults", "RoundRobin", "RankProduct", "SAM", "DESeq", "EdgeR"), 
				geneMapColumn=if(speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				groupColumn="Group", colorColumn="Color",
				geneSets=defaultGeneSets(speciesID), descriptor="CombinedGeneSets", 
				minGenesPerSet=if (speciesID %in% MAMMAL_SPECIES) 5 else 2, 
				mode=c("combined", "separate"), cutPvalue=0.05, cutRankShift=NULL, makePlots=TRUE,
				doFDR=TRUE, trimGenesToGeneMap=TRUE, makeGeneTables=!(speciesID %in% MAMMAL_SPECIES), 
				cutFold=0.1, cutFDR=0.05, NgeneSets=500, addCellTypes=(speciesID %in% MAMMAL_SPECIES), 
				addLifeCycle=(speciesID %in% PARASITE_SPECIES), PLOT.FUN=NULL, verbose=T, ...)
{

	toolName <- match.arg( toolName)
	mode <- match.arg( mode)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting pipe 'GeneSetDensity' on Sample Set: \n")
		print(sampleIDset)
		cat(  "\nUsing results from Species:  ", speciesID)
		cat(  "\nUsing results from DE tool:  ", toolName)
		if (is.character(geneSets)) cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}

	annT <- readAnnotationTable( annotationFile)
	allSamples <- unique( unlist( sampleIDset))
	myAnnT <- subset( annT, SampleID %in% allSamples)
	if ( ! all( c( groupColumn, colorColumn) %in% colnames(annT))) {
		cat( "\nGroup and/or Color column name not in Annotation table.")
		cat( "\nLooked for:  ", groupColumn, " ", colorColumn)
		cat( "\nFound:  ", colnames(annT))
		stop()
	}

	allGroups <- unique( myAnnT[[ groupColumn]])

	DE_list <- readDEgroupsData( groupIDset=allGroups, speciesID=speciesID, 
				optionsFile=optionsFile, results.path=results.path, folderName=folderName,
				toolName=toolName)
	
	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode=mode)
	nbig <- GSanswer$n
	bigSetOfDescriptors <- GSanswer$descriptor
	bigSetOfGeneSets <- GSanswer$geneSets

	# map from the group names to find the colors to use...
	where <- base::match( names( DE_list), myAnnT[[ groupColumn]])
	mycolors <- myAnnT[[ colorColumn]][ where]

	for ( i in 1:nbig) {
		geneSetDensity( DE_list, bigSetOfGeneSets[[i]], speciesID=speciesID, 
				colorset=mycolors, optionsFile=optionsFile, results.path=results.path, 
				folderName=folderName, toolName=toolName, geneMapColumn=geneMapColumn, 
				descriptor=bigSetOfDescriptors[i], minGenesPerSet=minGenesPerSet,
				cutPvalue=cutPvalue, cutRankShift=cutRankShift, doFDR=doFDR, 
				cutFold=cutFold, cutFDR=cutFDR, NgeneSets=NgeneSets,
				trimGenesToGeneMap=trimGenesToGeneMap, makeGeneTables=makeGeneTables,
				makePlots=makePlots, addCellTypes=addCellTypes, addLifeCycle=addLifeCycle,
				PLOT.FUN=PLOT.FUN)
	}

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished pipe 'GeneSetDensity' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n")
	}

	return()
}


`matrix.GeneSetDensity` <- function( x, groups=colnames(x), colors=1:ncol(x), speciesID=getCurrentSpecies(), 
				optionsFile="Options.txt", results.path=NULL,  folderName="", 
				geneSets=defaultGeneSets(speciesID), descriptor="CombinedGeneSets", 
				minGenesPerSet=if (speciesID %in% MAMMAL_SPECIES) 5 else 2, 
				geneMapColumn=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				mode=c("combined", "separate"), cutPvalue=0.05, 
				cutFold=0.1, cutFDR=0.05, NgeneSets=500, cutRankShift=NULL,
				doFDR=TRUE, trimGenesToGeneMap=TRUE, makeGeneTables=!(speciesID %in% MAMMAL_SPECIES), 
				makePlots=TRUE, addCellTypes=(speciesID %in% MAMMAL_SPECIES), 
				addLifeCycle=(speciesID %in% PARASITE_SPECIES), PLOT.FUN=NULL, verbose=T)
{

	mode <- match.arg( mode)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting 'GeneSetDensity' on a data matrix")
		cat(  "\nStated Species:              ", speciesID)
		cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}
	setCurrentSpecies( speciesID)

	# turn the columns of X into the DE data frames for each group
	DE_list <- expressionMatrixToDElist( x, groups=groups)

	# use the group names to get the right colors
	mygrps <- names(DE_list)
	mycolors <- colors[ match( mygrps, groups)]

	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode=mode)
	nbig <- GSanswer$n
	bigSetOfDescriptors <- GSanswer$descriptor
	bigSetOfGeneSets <- GSanswer$geneSets

	for ( i in 1:nbig) {
		geneSetDensity( DE_list, bigSetOfGeneSets[[i]], speciesID=speciesID, 
				colorset=mycolors, optionsFile=optionsFile, results.path=results.path, 
				folderName=folderName, toolName=NULL, geneMapColumn=geneMapColumn,
				descriptor=bigSetOfDescriptors[i], minGenesPerSet=minGenesPerSet,
				cutPvalue=cutPvalue, cutFold=cutFold, cutFDR=cutFDR, cutRankShift=cutRankShift, 
				doFDR=doFDR, NgeneSets=NgeneSets,
				trimGenesToGeneMap=trimGenesToGeneMap, makeGeneTables=makeGeneTables,
				makePlots=makePlots, addCellTypes=addCellTypes, addLifeCycle=addLifeCycle,
				PLOT.FUN=PLOT.FUN)
	}

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished 'GeneSetDensity' on data matrix:")
	}

	return()
}


`fileSet.GeneSetDensity` <- function( files, groups, colors=1:length(files), speciesID=getCurrentSpecies(), sep="\t",
				results.path=".",  geneColumn="GENE_ID", log2foldColumn="LOG2FOLD", pvalueColumn="PVALUE",
				geneSets=defaultGeneSets(speciesID), descriptor="CombinedGeneSets", 
				minGenesPerSet=if (speciesID %in% MAMMAL_SPECIES) 5 else 2, 
				geneMapColumn=if (speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
				mode=c("combined", "separate"), cutPvalue=0.05, cutFold=0.1, cutFDR=0.05, cutRankShift=NULL, 
				doFDR=TRUE, NgeneSets=500, trimGenesToGeneMap=TRUE, makeGeneTables=!(speciesID %in% MAMMAL_SPECIES), 
				makePlots=TRUE, addCellTypes=(speciesID %in% MAMMAL_SPECIES), 
				addLifeCycle=(speciesID %in% PARASITE_SPECIES), PLOT.FUN=NULL, verbose=T)
{

	mode <- match.arg( mode)

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting 'GeneSetDensity' on a file set")
		cat(  "\nStated Species:              ", speciesID)
		cat(  "\nGene Sets to analyze:        ", geneSets,"\n")
	}
	setCurrentSpecies( speciesID)

	NF <- length(files)
	if ( any( ! file.exists(files))) {
		cat( "\nSome files not found: \n")
		print( files[ ! file.exists(files)])
		return()
	}

	DE_list <- vector( mode='list', length=NF)
	for ( j in 1:NF) {
		tbl <- read.delim( files[j], as.is=T, sep=sep)
		if ( is.null( pvalueColumn)) {
			sml <- tbl[ , c( geneColumn, log2foldColumn)]
			colnames(sml) <- c( "GENE_ID", "LOG2FOLD")
			sml$PVALUE <- 1
		} else {
			sml <- tbl[ , c( geneColumn, log2foldColumn, pvalueColumn)]
			colnames(sml) <- c( "GENE_ID", "LOG2FOLD", "PVALUE")
		}
		DE_list[[j]] <- sml
	}
	names( DE_list) <- groups
	mycolors <- colors

	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode=mode)
	nbig <- GSanswer$n
	bigSetOfDescriptors <- GSanswer$descriptor
	bigSetOfGeneSets <- GSanswer$geneSets

	for ( i in 1:nbig) {
		geneSetDensity( DE_list, bigSetOfGeneSets[[i]], speciesID=speciesID, 
				colorset=mycolors, optionsFile="", results.path=results.path, 
				folderName="", toolName=NULL, geneMapColumn=geneMapColumn, 
				descriptor=bigSetOfDescriptors[i], minGenesPerSet=minGenesPerSet,
				cutPvalue=cutPvalue, cutFold=cutFold, cutFDR=cutFDR, cutRankShift=cutRankShift, 
				doFDR=doFDR, NgeneSets=NgeneSets,
				trimGenesToGeneMap=trimGenesToGeneMap, makeGeneTables=makeGeneTables,
				makePlots=makePlots, addCellTypes=addCellTypes, addLifeCycle=addLifeCycle,
				PLOT.FUN=PLOT.FUN)
	}

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished 'GeneSetDensity' on file set.")
	}

	return()
}

