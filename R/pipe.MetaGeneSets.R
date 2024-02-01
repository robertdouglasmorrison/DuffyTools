# pipe.MetaGeneSets.R

# run the 5 different DE tools and combine their results

`pipe.MetaGeneSets` <- function( sampleIDset, folderName="", speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  groupColumn="Group", colorColumn="Color", 
				geneSets=defaultGeneSets(speciesID), 
				NgeneSets=100, verbose=TRUE, label="", nFDRsimulations=0,
				doGeneSets=TRUE, doMissing=TRUE, baselineGroup=NULL, legend.cex=1, ...) {

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\nStarting pipe 'MetaGeneSets' on folder: \t", folderName)
		cat("\n", label, "\n\nUsing results from Species:  ", speciesID,"\n")
	}

	# set up for this species...
	annT <- readAnnotationTable( annotationFile)
	if ( ! (groupColumn %in% colnames(annT))) {
		cat( "\n\nError: Sample grouping column: '", groupColumn, "' not found in annotation file.", sep="")
		stop()
	}
	optT <- readOptionsTable( optionsFile)
	#if ( is.null( targetID)) targetID <- getOptionValue( optT, "targetID", notfound="HsPf", verbose=F)
	#setCurrentTarget( targetID)
	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optT, "results.path", notfound=".")
	}

	dev.type <- getPlotDeviceType( optT)
	dev.ext <- paste( ".", dev.type, sep="")
	
	# check that the setup of the cell type tools has occured
	CellTypeSetup()
	reference <- getCellTypeReference()
	doCellType <- FALSE
	if ( ! is.null( reference) && reference != "") {
		doCellType <- TRUE
	} else {
		reference <- "CellTypes"
	}

	# rather than force running all the GeneSet tools, just search for and use what you find.  Report missing ones..

	wt.fold <- 1
	wt.pval <- 0.05
	wt.rank <- 2

	metaPath <- file.path( results.path, "MetaResults", paste( prefix, folderName, sep="."))
	if ( ! file.exists( metaPath)) {
		cat( "\nGiven Meta Results folder not found.  Tried: ", metaPath)
		cat( "\nNo Meta GeneSets done...")
		return(NULL)
	}
	metaPathwayPath <- file.path( metaPath, "MetaGeneSets")
	if ( ! file.exists( metaPathwayPath)) dir.create( metaPathwayPath, recursive=T, showWarn=F)

	# see which tools have results here
	allFiles <- dir( metaPath, include.dir=T, full.name=T)
	subFolders <- allFiles[ file.info(allFiles)$isdir]
	radarFolder <- grep( "/RadarPlots", subFolders, value=T)[1]
	hasRadar <- ( ! is.na( radarFolder))
	if ( doGeneSets || (!hasRadar && doMissing)) {
		pipe.RadarPlots( sampleIDset, folderName=folderName, speciesID=getCurrentSpecies(), 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,  
				groupColumn=groupColumn, colorColumn=colorColumn, baselineGroup=baselineGroup,
				geneSets=geneSets, main=paste( "Comparison:  ",folderName), legend.cex=legend.cex, ...)
		allFiles <- dir( metaPath, include.dir=T, full.name=T)
		subFolders <- allFiles[ file.info(allFiles)$isdir]
		radarFolder <- grep( "/RadarPlots", subFolders, value=T)[1]
		hasRadar <- ( ! is.na( radarFolder))
	}
	# do Bubble plots if/when we do Radar Plots
	bubbleFolder <- grep( "/BubblePlots", subFolders, value=T)[1]
	hasBubble <- ( ! is.na( bubbleFolder))
	if ( doGeneSets || (!hasBubble && doMissing)) {
		pipe.GeneSetBubblePlots( sampleIDset, folderName=folderName, speciesID=getCurrentSpecies(), 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,  
				groupColumn=groupColumn, colorColumn=colorColumn, baselineGroup=baselineGroup,
				geneSets=geneSets, main=paste( "Comparison:  ",folderName), legend.cex=legend.cex, ...)
		# don't add to meta results, just make them
	}
	# allow old and new naming of this tool...
	densityFolder <- grep( "/Density", subFolders, value=T)[1]
	if ( is.na( densityFolder)) densityFolder <- grep( "/CombinedGeneSets", subFolders, value=T)[1]
	hasDensity <- ( ! is.na( densityFolder))
	if ( doGeneSets || (!hasDensity && doMissing)) {
		pipe.GeneSetDensity( sampleIDset, folderName=folderName, speciesID=getCurrentSpecies(), 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,  
				groupColumn=groupColumn, colorColumn=colorColumn, doFDR=F, NgeneSets=NgeneSets, 
				geneSets=geneSets, legend.cex=legend.cex, ...)
		allFiles <- dir( metaPath, include.dir=T, full.name=T)
		subFolders <- allFiles[ file.info(allFiles)$isdir]
		densityFolder <- grep( "/Density", subFolders, value=T)[1]
		if ( is.na( densityFolder)) densityFolder <- grep( "/CombinedGeneSets", subFolders, value=T)[1]
		hasDensity <- ( ! is.na( densityFolder))
	}
	enrichFolder <- grep( "/Enrichment", subFolders, value=T)[1]
	hasEnrich <- ( ! is.na( enrichFolder))
	if ( doGeneSets || (!hasEnrich && doMissing)) {
		pipe.GeneSetEnrichment( sampleIDset, folderName=folderName, speciesID=getCurrentSpecies(), 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,  
				geneSets=geneSets, groupColumn=groupColumn, ...)
		allFiles <- dir( metaPath, include.dir=T, full.name=T)
		subFolders <- allFiles[ file.info(allFiles)$isdir]
		enrichFolder <- grep( "/Enrichment", subFolders, value=T)[1]
		hasEnrich <- ( ! is.na( enrichFolder))
	}
	qusageFolder <- grep( "/QuSage", subFolders, value=T)[1]
	hasQuSage <- ( ! is.na( qusageFolder))
	if ( doGeneSets || (!hasQuSage && doMissing)) {
		pipe.QuSage( sampleIDset, folderName=folderName, speciesID=getCurrentSpecies(), 
				annotationFile=annotationFile, optionsFile=optionsFile, results.path=results.path,  
				geneSets=geneSets, groupColumn=groupColumn, ...)
		allFiles <- dir( metaPath, include.dir=T, full.name=T)
		subFolders <- allFiles[ file.info(allFiles)$isdir]
		qusageFolder <- grep( "/QuSage", subFolders, value=T)[1]
		hasQuSage <- ( ! is.na( qusageFolder))
	}

	flatSamples <- sort( unique( unlist( sampleIDset)))
	annT2 <- subset( annT, SampleID %in% flatSamples)
	myGrps <- sort( unique( annT2[[ groupColumn]]))

	for ( grp in myGrps) {

		cat( "\n\nDoing MetaGeneSets on:    ", grp)
		for (direction in c( "UP", "DOWN")) {
		cat( "\n\tDirection:  ", direction)

		# keep the full and short pathway names
		bigFullName <- bigShortName <- bigGeneCount <- vector()
		dfList <- vector( mode='list')
		nDF <- 0

		# set up to try to add extra links to info from Density/Radar/etc
		densityAns <- radarAns <- NULL

		if ( hasDensity) {
			# Combined Gene Sets has what was Significant and UP for each Group
			f <- file.path( densityFolder, paste( grp, prefix, direction, "CombinedGeneSets.txt", sep="."))
			if ( ! file.exists(f)) {
				cat( "\nCombinedGeneSets results not found:\n", f)
			} else {
				tmp <- read.delim(f, as.is=T)
				densityAns <- tmp
				myColumns <- if( "CellType" %in% colnames(tmp)) c(2,6,7,4) else c(2,5,6,3)
				sml <- tmp[ , myColumns]
				colnames(sml) <- c( "PathName","LOG2FOLD","PVALUE","N_GENES")
				bigFullName <- c( bigFullName, sml$PathName)
				bigGeneCount <- c( bigGeneCount, sml$N_GENES)
				sml$PathName <- cleanGeneSetModuleNames( sml$PathName, wrapParen=F)  # already cleaned...
				bigShortName <- c( bigShortName, sml$PathName)
				nDF <- nDF + 1
				dfList[[nDF]] <- sml
				names(dfList)[nDF] <- "Density GS"
			}
		}
		if ( hasQuSage) {
			# QuSage results are always UP only
			f <- file.path( qusageFolder, paste( grp, prefix, "UP", "QuSage.CombinedGeneSets.csv", sep="."))
			if ( ! file.exists(f)) {
				cat( "\nQuSage results not found:\n", f)
			} else {
				tmp <- read.csv(f, as.is=T)
				sml <- tmp[ ,c("PATHWAY","LOG2FOLD","PVALUE", "N_GENES")]
				colnames(sml)[1] <- "PathName"
				bigFullName <- c( bigFullName, sml$PathName)
				bigGeneCount <- c( bigGeneCount, sml$N_GENES)
				sml$PathName <- cleanGeneSetName( cleanGeneSetModuleNames( sml$PathName, wrapParen=F))
				bigShortName <- c( bigShortName, sml$PathName)
				# since it has UP only, invert the signs for DOWN
				if ( direction == "DOWN") {
					 sml$LOG2FOLD <- -(sml$LOG2FOLD)
					# the sign change causes us to rerank
					ord <- diffExpressRankOrder( sml$LOG2FOLD, sml$PVALUE, wt.fold=1, wt.pvalue=2)
					sml <- sml[ ord, ]
					if (nrow(sml)) rownames(sml) <- 1:nrow(sml)
				}
				nDF <- nDF + 1
				dfList[[nDF]] <- sml
				names(dfList)[nDF] <- "QuSage GS"
			}
		}
		if ( hasRadar) {
			# Radar Plots only have a single combined file, with all groups and all directions
			f <- file.path( radarFolder, "Radar.All.CombinedGeneSets.csv")
			if ( ! file.exists(f)) {
				if ( ! file.exists( radarFolder)) {
					cat( "\nNo RadarPlots results found: \n", radarFolder)
				} else {
					combineRadarPlotResults( radarFolder)
				}
			}
			if ( ! file.exists(f)) {
				cat( "\nRadarPlot results not found:\n", f)
			} else {
				tmp <- read.csv(f, as.is=T)
				radarAns <- tmp
				grpPattern <- paste( "_", grp, "$", sep="")
				isGrp <- grep( grpPattern, colnames(tmp))
				isFold <- grep( "^Fold_", colnames(tmp))
				isPval <- grep( "^Pvalue_", colnames(tmp))
				myColumns <- c( 1, intersect(isGrp,isFold)[1], intersect(isGrp,isPval)[1], 2)
				sml <- tmp[ , myColumns]
				colnames(sml) <- c("PathName","LOG2FOLD","PVALUE", "N_GENES")
				bigFullName <- c( bigFullName, sml$PathName)
				bigGeneCount <- c( bigGeneCount, sml$N_GENES)
				sml$PathName <- cleanGeneSetName( cleanGeneSetModuleNames( sml$PathName, wrapParen=F))
				bigShortName <- c( bigShortName, sml$PathName)
				# no separate "DOWN' data, so just invert
				if ( direction == "DOWN") {
					sml$LOG2FOLD <- -(sml$LOG2FOLD)
				}
				# sort on just this group
				ord <- diffExpressRankOrder( sml$LOG2FOLD, sml$PVALUE, wt.fold=1, wt.pvalue=2)
				sml <- sml[ ord, ]
				if (nrow(sml)) rownames(sml) <- 1:nrow(sml)
				nDF <- nDF + 1
				dfList[[nDF]] <- sml
				names(dfList)[nDF] <- "Radar GS"
			}
		}
		if ( hasEnrich) {
			f <- file.path( enrichFolder, paste( grp, direction, "Enrichment.csv", sep="."))
			if ( ! file.exists(f)) {
				cat( "\nEnrichment results not found:\n", f)
			} else {
				tmp <- read.csv(f, as.is=T)
				sml <- tmp[ , c( 1,grep("Enrichment",colnames(tmp)),grep("PVALUE",colnames(tmp)))]
				colnames(sml) <- c( "PathName", "LOG2FOLD", "PVALUE")
				# turn the enrichment into something more like Log2Fold units, and flip it's sign for DOWN
				sml$LOG2FOLD <- log2( sml$LOG2FOLD)
				if ( direction == "DOWN") {
					sml$LOG2FOLD <- -(sml$LOG2FOLD)
				}
				# remove the N= suffix  (already cleaned of any hyperlinks...)
				sml$PathName <- sub( "  \\(N=[0-9]+\\)$", "", sml$PathName)
				sml$PathName <- cleanGeneSetName( sml$PathName)
				nDF <- nDF + 1
				dfList[[nDF]] <- sml
				names(dfList)[nDF] <- "Enrich GS"
			}
		}
		
		if ( length(dfList) < 2) {
			cat( "\nNot enough Gene Set analyses found.")
			cat( "\nSkipping run of MetaGeneSets...")
			return()
		}

		ans <- metaRank.data.frames( dfList, geneColumn="PathName", valueColumn="LOG2FOLD",
					pvalueColumn="PVALUE", productColumn=NA, missingGenes="fill", 
					missingValue=0, direction=direction, nFDRsimulations=nFDRsimulations)
	
		# since these are pathways, drop the "PRODUCT" column that MetaRanks gave us
		if ( "PRODUCT" %in% colnames(ans)) ans <- ans[ ,-grep("PRODUCT",colnames(ans))]

		# put back in  the fully hyperlinked path names if we can, and the gene counts
		# and try to cache up the various pathway naming styles for adding extra links
		ans$N_Genes <- NA
		where <- match( ans$PathName, bigShortName, nomatch=0)
		ansShortPathNames <- ans$PathName
		ans$PathName[ where > 0] <- bigFullName[where]
		ansLongPathNames <- ans$PathName
		ans$N_Genes[ where > 0] <- bigGeneCount[where]

		# append extra links where we can
		radarRanks <- as.numeric( ans[["Radar GS"]])
		namesWithExtraLinks <- addExtraPathwayLinks( ansShortPathNames, ansLongPathNames, 
							densityAns=densityAns, radarAns=radarAns, 
							radarRanks=radarRanks, metaPath=metaPath,
							dev.ext=dev.ext, densityFolder=densityFolder)
		ans$PathName <- namesWithExtraLinks

		# put in the column order we expect
		out <- ans[ ,c( 1, ncol(ans), 2,4,3, 5:(ncol(ans)-1))]

		if ( doCellType) {
			out$CellType <- geneSetCellType( ansLongPathNames, max.type=4)
			out <- out[ ,c( 1, ncol(out), 2:(ncol(out)-1))]
		}

		# lastly, we may have a better ordering by using all 3 features
		myFold <- out$LOG2FOLD
		if ( direction == "DOWN") myFold <- -myFold
		myPval <- out$AVG_PVALUE
		myRank <- out$AVG_RANK
		ord <- diffExpressMetaResultOrder( myFold, myPval, myRank, wt.fold=wt.fold, wt.pvalue=wt.pval, wt.rank=wt.rank)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)

		# make the text file version of all rows
		fileout <- paste( grp, prefix, "MetaGeneSets", direction, "txt", sep=".")
		fileout <- file.path( metaPathwayPath, fileout)
		write.table( out, fileout, sep="\t", quote=F, row.names=F)

		# final answer is cropped
		N <- min( nrow(out), NgeneSets)
		out <- out[ 1:N, ]

		htmlout <- sub( ".txt$", ".html", fileout)
		otherGrpTxt <- setdiff( myGrps, grp)
		if ( length(otherGrpTxt) > 1) otherGrpTxt <- paste( "{", paste( otherGrpTxt, collapse=" + "), "}")
		title <- paste( "MetaGeneSets: &nbsp; GeneSets ", direction, " Regulated in group: &nbsp; ", grp, 
				" &nbsp; vs &nbsp; ", otherGrpTxt, "<br> Comparison Folder: &nbsp; ", folderName)
		metaResultsToHTML( out, htmlout, addSpeciesToHtmlTitle(title), maxRows=N, linkColumnName=NA)

		# write the cell type enrichment on just these top N gene sets...
		if (doCellType) {
			enrich <- cellTypeEnrichment( out$CellType, mode="geneSets", upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( grp, prefix, "MetaGeneSets", direction, reference, "Enrichment.csv", sep=".")
			f <- file.path( metaPathwayPath, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
		rm( out)
		cat("\n")
	}}  # for direction and group

	if (verbose) {
		cat( verboseOutputDivider)
		cat( "\n\nFinished pipe 'MetaGeneSets' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\n", label, "\n")
	}
	return()
}


`addExtraPathwayLinks` <- function( ansShortPathNames, ansLongPathNames, densityAns=NULL, 
					radarAns=NULL, radarRanks=NULL, metaPath=".", dev.ext=".png",
					densityFolder="./Density") {

	# given the pathway names, with and without the hyperlinks to the source URL
	N <- length( ansShortPathNames)
	out <- ansLongPathNames

	# add others if we can
	outDensity <- outGenes <- outRadar <- rep.int( "", N)
	if ( ! is.null( densityAns)) {

		# grab the fields from the Density data
		gsID <- densityAns$GenesPerGroup
		longnames <- densityAns$Name
		shortnames <- cleanGeneSetModuleNames(longnames)
		gsNumber <- sub( "GeneSet_", "", gsID)
		
		# see if we hit any of them
		myNumber <- rep.int( "", N)
		whoShort <- match( ansShortPathNames, shortnames, nomatch=0)
		myNumber[ whoShort > 0] <- gsNumber[ whoShort]
		whoLong <- match( ansLongPathNames, longnames, nomatch=0)
		myNumber[ whoLong > 0] <- gsNumber[ whoLong]

		# for the Density results, the files should always be there if
		# but the folder name has changed
		relativeDensityFolder <- "../Density"
		localDensityFolder <- "Density"
		localDensityPlotsFolder <- "Density.pngPlots"
		if ( grepl( "CombinedGeneSets", densityFolder)) {
			relativeDensityFolder <- "../CombinedGeneSets"
			localDensityFolder <- "CombinedGeneSets"
			localDensityPlotsFolder <- "CombinedGeneSets.pngPlots"
		}

		hits <- which( myNumber != "")
		if ( length( hits)) {
			fDensity <- paste( "CombinedGeneSets_", myNumber[hits], dev.ext, sep="")
			fLocalDensity <- file.path( relativeDensityFolder, localDensityPlotsFolder, fDensity)
			outDensity[hits] <- fLocalDensity
			fGlobalDensity <- file.path( metaPath, localDensityFolder, localDensityPlotsFolder, fDensity)
			notFile <- ( ! file.exists( fGlobalDensity))
			outDensity[ hits[ notFile]] <- ""
			fGenes <- paste( "GeneSet_", myNumber[hits], ".html", sep="")
			fLocalGenes <- file.path( relativeDensityFolder, localDensityPlotsFolder, fGenes)
			outGenes[hits] <- fLocalGenes
			fGlobalGenes <- file.path( metaPath, localDensityFolder, localDensityPlotsFolder, fGenes)
			notFile <- ( ! file.exists( fGlobalGenes))
			outGenes[ hits[ notFile]] <- ""
		}
	}

	if ( ! is.null( radarAns)) {

		# grab the fields from the Radar data
		myNames <- radarAns$PathName
		whoRadar <- match( ansShortPathNames, myNames, nomatch=0)

		# the image name is from the very front
		gsID <- sub( ":.+", "", myNames)
		fRadar <- paste( "Radar.", gsID, dev.ext, sep="")
		fRadar <- file.path( "../RadarPlots", fRadar)
		outRadar[ whoRadar > 0] <- fRadar[whoRadar]

		# radars only show the top K, so try to not add links that would not be seen
		if ( ! is.null( radarRanks)) {
			drops <- which( radarRanks > 24)
			if (length(drops)) outRadar[ drops] <- ""
		}
	}

	# now we know what links we can add
	use <- which( outDensity != "")
	if ( length(use)) {
		url <- outDensity[use]
		newText <- as.link( url, "(density)")
		out[use] <- paste( out[use], newText, sep=" &nbsp; ")
	}
	use <- which( outRadar != "")
	if ( length(use)) {
		url <- outRadar[use]
		newText <- as.link( url, "(radar)")
		out[use] <- paste( out[use], newText, sep=" &nbsp; ")
	}
	use <- which( outGenes != "")
	if ( length(use)) {
		url <- outGenes[use]
		newText <- as.link( url, "(genes)")
		out[use] <- paste( out[use], newText, sep=" &nbsp; ")
	}

	# done!
	return( out)
}
