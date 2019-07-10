# geneSetDensity.R

# various tools to mine and display the DE of gene sets from Pathways or GO, etc


`geneSetDensity` <- function( deList, geneSets, speciesID=getCurrentSpecies(), colorset=c( 2:(length(deList)+1)),
					optionsFile="Options.txt", results.path=NULL, folderName="", 
					toolName=c("MetaResults","RoundRobin","RankProduct","SAM","DESeq","EdgeR"),
					descriptor="GeneSet", minGenesPerSet=2, trimGenesToGeneMap=TRUE,
					geneMapColumn=if(speciesID %in% MAMMAL_SPECIES) "NAME" else "GENE_ID", 
					cutPvalue=0.05, cutRankShift=NULL, cutFold=0.1, cutFDR=0.05, 
					makePlots=TRUE, PLOT.FUN=NULL, makeGeneTables=TRUE, 
					addCellTypes=(speciesID %in% MAMMAL_SPECIES), 
					addLifeCycle=(speciesID %in% PARASITE_SPECIES), 
					doFDR=TRUE, NgeneSets=500) {

	setCurrentSpecies( speciesID)
	geneMap <- getCurrentGeneMap()
	geneMap <- subset.data.frame( geneMap, REAL_G == TRUE)
	rownames(geneMap) <- 1:nrow(geneMap)

	if ( is.null(cutRankShift)) cutRankShift <- sum( geneMap$REAL_G) * 0.05

	subsetNames <- names( geneSets)

	# try to catch and drop any missing results
	drops <- which( sapply( deList, nrow) < 1)
	if ( length( drops)) {
		cat( "Dropping NotFound/Missing Groups:  ", names(deList)[drops])
		deList <- deList[ -drops]
	}

	groupIDs <- names( deList)
	Ngrps <- length(groupIDs)
	cat( "\n\nGene Sets for:   ", descriptor, "\n\nAnalyzing ", length(subsetNames), 
		" gene sets for significance among ", Ngrps, " groups:  ", groupIDs, "\n")


	# see how many genes per set,... compare needs at least 2
	# use actual number of genes in the geneMap, not just the length of the given list
	if ( trimGenesToGeneMap) {
		validNames <- geneMap[[ geneMapColumn]]
		geneSets <- lapply( geneSets, function(x) return( x[ x %in% validNames]))
	}
	nEach <- sapply( geneSets, length)
	canBeEvaluated <- which( nEach >= minGenesPerSet)
	cat( "  Skipping ", sum( nEach < minGenesPerSet), " gene sets with less than", minGenesPerSet, " genes.\n")

	# set up directories to read from or write to...
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	}
	prefix <- getCurrentSpeciesFilePrefix()

	if ( ! is.null( toolName)) {
		toolName <- match.arg( toolName)
		toolPrefix <- c( "RR", "RP", "SAM", "Meta", "DESeq", "EdgeR")[ match( toolName, 
				c("RoundRobin","RankProduct","SAM","MetaResults","DEseq", "EdgeR"))]
		if ( base::nchar( folderName) < 1) stop( "'geneSetDensity' needs an explicit 'folderName' argument...")
		DE_path <- file.path( results.path, toolName, paste( prefix, folderName, sep="."))
	} else {
		DE_path <- file.path( results.path, folderName)
	}
	if ( ! file.exists( DE_path)) dir.create( DE_path, recursive=T, showWarnings=F)
	DE_colors <- colorset
	GS_path <- file.path( DE_path, descriptor)
	if ( ! file.exists( GS_path)) dir.create( GS_path, recursive=T, showWarnings=F)

	# clean up any global storage... 
	if ( exists( "DE_List")) try( rm( DE_List, envir=.GlobalEnv))

	out <- data.frame()
	goodSets <- vector() 
	goodProfiles <- vector( mode="list")

	if (doFDR) {
		randomGeneSetPvalues <- calcRandomGeneSetPvalues( deList, geneSets, who=canBeEvaluated)
	} else {
		randomGeneSetPvalues <- NULL
	}

	cat( "\nMeasuring gene set densities...\n")

	# we previously tossed out bad P-value sets right away...
	# but some downstream tools would like to know exactly how each set fared.
	# so keep all for the .TXT file versions, and trim down for the .HTML versions
	for( i in canBeEvaluated) {

		densityProfiles <- calcGeneSetDensity( deList, geneSets[[i]], geneMapColumn)
		ans <- measureDensityDifferences( densityProfiles, randomPvalues=randomGeneSetPvalues)
		# skip if not worthwhile
		#if ( ans$bestPvalue > cutPvalue) next

		# good, add to table
		if ( ans$bestPvalue <= cutPvalue) {
			goodSets <- c( goodSets, i)
			goodProfiles[[i]] <- densityProfiles
		}

		oneDF <- data.frame( "GeneSet"=paste( descriptor, i, sep="_"), 
				"Name"=cleanGeneSetName( names( geneSets)[i]), 
				"N Genes"=length( geneSets[[i]]), "GenesPerGroup"=paste( "GeneSet", i, sep="_"), 
				t(ans$measurements), stringsAsFactors=FALSE)
		colnames( oneDF)[1] <- descriptor
		out <- rbind( out, oneDF)
		if ( i %% 10 == 0) cat( "\r", i, trimGeneSetNameLink( names(geneSets)[i]), "   ", length( geneSets[[i]]))
	}

	# now with them all done, adjust for multiple comparisons
	bestPvalues <- rep.int( 1, nrow(out))
	pvalColumns <- grep( "P.value", colnames(out))
	for ( k in pvalColumns) {
		pvals <- as.numeric( out[[k]])
		pvals <- p.adjust( pvals, "BH")
		out[[k]] <- pvals
		bestPvalues <- pmin( bestPvalues, pvals)
	}

	# let's gather who will get tossed out, but don't delete those rows until we get to HTML steps.
	drops <- which( bestPvalues > cutPvalue)
	cat( "\nDropped by 'cutPvalue': ", length(drops))
	bestFolds <- rep.int( 0, nrow(out))
	foldColumns <- grep( "FoldChange", colnames(out))
	for ( k in foldColumns) {
		folds <- as.numeric( out[[k]])
		bestFolds <- pmax( bestFolds, abs(folds))
	}
	dropsFold <- which( bestFolds < cutFold)
	cat( "\nDropped by 'cutFold': ", length(dropsFold))
	drops <- sort( union( drops, dropsFold))
	if ( doFDR) {
		bestFDR <- rep.int( 1, nrow(out))
		fdrColumns <- grep( "FDR", colnames(out))
		for ( k in fdrColumns) {
			fdrs <- as.numeric( out[[k]])
			bestFDR <- pmin( bestFDR, fdrs)
		}
		dropsFDR <- which( bestFDR > cutFDR)
		cat( "\nDropped by 'cutFDR': ", length(dropsFDR))
		drops <- sort( union( drops, dropsFDR))
	}

	# the result is ready now.   Do we want to add in a column of "Cell Types"?
	if ( addCellTypes) {
		cellType <- getGeneSetCellType( out$Name)
		out <- cbind( out[,1:2], "CellType"=cellType, out[,3:ncol(out)], stringsAsFactors=F)
	} else if ( addLifeCycle) {
		lifeCycle <- getGeneSetLifeCycle( out$Name)
		out <- cbind( out[,1:2], "LifeCycle"=lifeCycle, out[,3:ncol(out)], stringsAsFactors=F)
	}

	# remove the entries with not enough change
	# we are switching over to just keeping one "all sets" result, and turning it into HTML on the fly
	outForHTML <- out

	# be smarter about getting the true IDs of who is getting dropped.  Simple indexes are not valid!...
	if ( length(drops)) {
		smlDrop <- out[ drops, ]
		IDsToDrop <- as.integer( sub( "^.+_", "", smlDrop[[1]]))
		#cat( "\n\nDEBUG:  ", head( IDsToDrop, 100), "\n\n")
		# the data frame rows are by index
		outForHTML <- outForHTML[ -drops, ]

		# but the set of IDs is by set membership
		#goodSets <- goodSets[ -drops]
		goodSets <- setdiff( goodSets, IDsToDrop)
	}
	cat( "\nN_GeneSets with multiple comparison adjusted Pvalue < ", cutPvalue, "and abs(Fold) > ", 
		cutFold, " = ", length(goodSets))
	if ( length(goodSets) < 1) {
		cat( "\nNo significant results.  Perhaps loosen cutoffs...")
		return()
	}

	# make a folder to hold all the linked plots and geneset files
	localPlotPath <- paste( descriptor, "pngPlots", sep=".")
	globalPlotPath <- file.path( GS_path, localPlotPath)
	if ( ! file.exists( globalPlotPath)) dir.create( globalPlotPath)
	# plots to the genes are up relative to the HTML
	localGenePlotPath <- "../../pngPlots"

	# write the result as a text file
	fileout=paste( prefix, descriptor, "txt", sep=".")
	fileout <- file.path( GS_path, fileout)
	write.table( out, file=fileout, sep="\t", quote=FALSE, row.names=FALSE)

	# write the result as a set of web pages and text for each group
	for( k in grep( "FoldChange", colnames(outForHTML))) {
		outForHTML[[k]] <- formatC( outForHTML[[k]], format="f", digits=3, flag="+")
	}
	for( k in grep( "P.value", colnames(outForHTML))) {
		outForHTML[[k]] <- formatC( outForHTML[[k]], format="e", digits=2)
	}
	for( k in grep( "FDR", colnames(outForHTML))) {
		outForHTML[[k]] <- formatC( outForHTML[[k]], format="f", digits=3)
	}
	for( k in grep( "RankShift", colnames(outForHTML))) {
		outForHTML[[k]] <- formatC( round(outForHTML[[k]]), format="d", flag="+", big.mark=",")
	}
	colnames(outForHTML) <- sub( "FoldChange", "Fold Change", colnames(outForHTML))
	colnames(outForHTML) <- sub( "RankShift", "Rank Shift", colnames(outForHTML))
	colnames(outForHTML) <- sub( "GenesPerGroup", "Genes Per Group", colnames(outForHTML))
	if (! doFDR) {
		k <- grep( "FDR", colnames(outForHTML))
		if ( length(k)) outForHTML <- outForHTML[ , -k]
	}
	if (makeGeneTables) {
		table2html( outForHTML, fileout=sub( "txt$", "html", fileout), title=addSpeciesToHtmlTitle(descriptor), 
			linkColumnNames=c( descriptor, "Genes Per Group"), 
			linkPaths=rep( localPlotPath, times=2), linkExtensions=c( ".png", ".html"))
	} else {
		table2html( outForHTML, fileout=sub( "txt$", "html", fileout), title=addSpeciesToHtmlTitle(descriptor), 
			linkColumnNames=descriptor, linkPaths=localPlotPath, linkExtensions=".png")
	}

	# also make smaller pages of each one sample, with only the relavent rows...
	geneSetsToPlot <- vector()
	# always start from the main result
	nColumnPerSample <- nColumnPerSample.html <- 4
	for ( j in 1:Ngrps) {
		thisSample <- groupIDs[j]
		if (Ngrps == 2) {
			otherGroup <- groupIDs[3-j]
		} else {
			otherGroup <- paste( "{", paste(sort(setdiff(groupIDs,thisSample)),collapse=" + "), "}", sep="")
		}
		GENEGROUP <- 4;  FOLDCHANGE <- 5;  PVALUE <- 6;  FDR <- 7;  RANKSHIFT <- 8;  
		if (addCellTypes || addLifeCycle) {
			GENEGROUP <- 5;  FOLDCHANGE <- 6;  PVALUE <- 7;  FDR <- 8;  RANKSHIFT <- 9;  
		}
		theseColumns <- ((j-1) * nColumnPerSample) + c(FOLDCHANGE:RANKSHIFT)
		theseColumns.html <- theseColumns
		RANKSHIFT.HTML <- RANKSHIFT
		if ( !doFDR) {
			nColumnPerSample.html <- 3
			theseColumns.html <- ((j-1) * nColumnPerSample.html) + c(FOLDCHANGE:FDR) # no 3rd column of 4 in this case
			RANKSHIFT.HTML <- RANKSHIFT.HTML - 1
		}

		# up rows
		sml <- out[ , c( 1:GENEGROUP, theseColumns)]
		ord <- diffExpressDistanceRankOrder( sml[ , FOLDCHANGE], sml[ , PVALUE], 
				sml[ , RANKSHIFT], wt.folds=1, wt.pvalues=0.5, wt.dist=1)
		sml <- sml[ ord, ]
		rownames(sml) <- 1:nrow(sml)
		# a copy as text
		f <- paste( thisSample, prefix, "UP", descriptor, "txt", sep=".")
		f <- file.path( GS_path, f)
		write.table( sml, f, sep="\t", quote=F, row.names=F)

		# only keep what we would draw
		keep <- which( sml[,RANKSHIFT] > cutRankShift & sml[,PVALUE] < cutPvalue & sml[,FOLDCHANGE] > cutFold)
		if ( length( keep) > 0) {
			smlForHTML <- sml[ keep, ]
			if ( nrow(smlForHTML) > NgeneSets) smlForHTML <- smlForHTML[ 1:NgeneSets, ]
			geneSetsToPlot <- c( geneSetsToPlot, smlForHTML[[1]])
			for( k in grep( "FoldChange", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( smlForHTML[[k]], format="f", digits=3, flag="+")
			}
			for( k in grep( "P.value", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( smlForHTML[[k]], format="e", digits=2)
			}
			for( k in grep( "FDR", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( smlForHTML[[k]], format="f", digits=3)
			}
			for( k in grep( "RankShift", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( round(smlForHTML[[k]]), format="d", flag="+", big.mark=",")
			}
			colnames(smlForHTML) <- sub( "FoldChange", "Fold Change", colnames(smlForHTML))
			colnames(smlForHTML) <- sub( "RankShift", "Rank Shift", colnames(smlForHTML))
			colnames(smlForHTML) <- sub( "GenesPerGroup", "Genes Per Group", colnames(smlForHTML))
			if (! doFDR) {
				k <- grep( "FDR", colnames(smlForHTML))
				if ( length(k)) smlForHTML <- smlForHTML[ , -k]
			}
			f <- paste( thisSample, prefix, "UP", descriptor, "html", sep=".")
			f <- file.path( GS_path, f)
			mytitle <- paste( descriptor, ": &nbsp; Gene sets Up-Regulated in group: &nbsp; ", 
						thisSample, " &nbsp; vs &nbsp; ", otherGroup, sep="") 
			table2html( smlForHTML, fileout=f, title=addSpeciesToHtmlTitle(mytitle), maxRows=NgeneSets,
				linkColumnNames=c( descriptor, "Genes Per Group"), 
				linkPaths=rep( localPlotPath, times=2), linkExtensions=c( ".png", ".html"))

			# we can also make an enrichment table of cell types
			if (addCellTypes) {
				enrich <- cellTypeEnrichment( smlForHTML$CellType, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
				f <- paste( thisSample, prefix, "UP.GeneSetCellTypeEnrichment.csv", sep=".")
				f <- file.path( GS_path, f)
				write.table( enrich, f, sep=",", quote=T, row.names=F)
			} else if (addLifeCycle) {
				enrich <- lifeCycleEnrichment( smlForHTML$LifeCycle, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
				f <- paste( thisSample, prefix, "UP.GeneSetLifeCycleEnrichment.csv", sep=".")
				f <- file.path( GS_path, f)
				write.table( enrich, f, sep=",", quote=T, row.names=F)
			}
		}
		# down rows
		sml <- out[ , c( 1:GENEGROUP, theseColumns)]
		ord <- diffExpressDistanceRankOrder( -(sml[ , FOLDCHANGE]), sml[ , PVALUE], 
				-(sml[ , RANKSHIFT]), wt.folds=1, wt.pvalues=0.5, wt.dist=1)
		sml <- sml[ ord, ]
		rownames(sml) <- 1:nrow(sml)
		# a copy as text too
		f <- paste( thisSample, prefix, "DOWN", descriptor, "txt", sep=".")
		f <- file.path( GS_path, f)
		write.table( sml, f, sep="\t", quote=F, row.names=F)

		keep <- which( sml[,RANKSHIFT] < (-cutRankShift) & sml[,PVALUE] < cutPvalue & sml[,FOLDCHANGE] < (-cutFold))
		if ( length( keep) > 0) {
			smlForHTML <- sml[ keep, ]
			if ( nrow(smlForHTML) > NgeneSets) smlForHTML <- smlForHTML[ 1:NgeneSets, ]
			geneSetsToPlot <- c( geneSetsToPlot, smlForHTML[[1]])
			for( k in grep( "FoldChange", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( smlForHTML[[k]], format="f", digits=3, flag="+")
			}
			for( k in grep( "P.value", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( smlForHTML[[k]], format="e", digits=2)
			}
			for( k in grep( "FDR", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( smlForHTML[[k]], format="f", digits=3)
			}
			for( k in grep( "RankShift", colnames(smlForHTML))) {
				smlForHTML[[k]] <- formatC( round(smlForHTML[[k]]), format="d", flag="+", big.mark=",")
			}
			colnames(smlForHTML) <- sub( "FoldChange", "Fold Change", colnames(smlForHTML))
			colnames(smlForHTML) <- sub( "RankShift", "Rank Shift", colnames(smlForHTML))
			colnames(smlForHTML) <- sub( "GenesPerGroup", "Genes Per Group", colnames(smlForHTML))
			if (! doFDR) {
				k <- grep( "FDR", colnames(smlForHTML))
				if ( length(k)) smlForHTML <- smlForHTML[ , -k]
			}
			f <- paste( thisSample, prefix, "DOWN", descriptor, "html", sep=".")
			f <- file.path( GS_path, f)
			mytitle <- paste( descriptor, ": &nbsp; Gene sets Down-Regulated in group: &nbsp; ", 
						thisSample, " &nbsp; vs &nbsp; ", otherGroup, sep="") 
			table2html( smlForHTML, fileout=f, title=addSpeciesToHtmlTitle(mytitle), maxRows=NgeneSets,
				linkColumnNames=c( descriptor, "Genes Per Group"), 
				linkPaths=rep( localPlotPath, times=2), linkExtensions=c( ".png", ".html"))

			# we can also make an enrichment table of cell types
			if (addCellTypes) {
				enrich <- cellTypeEnrichment( smlForHTML$CellType, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
				f <- paste( thisSample, prefix, "DOWN.GeneSetCellTypeEnrichment.csv", sep=".")
				f <- file.path( GS_path, f)
				write.table( enrich, f, sep=",", quote=T, row.names=F)
			} else if (addLifeCycle) {
				enrich <- lifeCycleEnrichment( smlForHTML$LifeCycle, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
				f <- paste( thisSample, prefix, "DOWN.GeneSetLifeCycleEnrichment.csv", sep=".")
				f <- file.path( GS_path, f)
				write.table( enrich, f, sep=",", quote=T, row.names=F)
			}
		}
	}
	# track just the gene set numbers that show up in HTML files
	geneSetsToPlot <- sort( unique( sub( ".+_", "", geneSetsToPlot)))

	# make all the little files of genes names and plots for each set
	if ( makeGeneTables) {

	   # files of genes in each pathway set
	   genesPathsSeen <- vector()
	   cat( "\n\nBuilding gene set tables...\n")

	   # we need to reset a few items for these tables
	   givenGeneMapColumn <- geneMapColumn

	   for( j in goodSets) {
		gset <- geneSets[[j]]

		# we were using the geneMap to grab these Product terms etc.
		# try to be a bit less restrictive...
		gmap <- makeSmallGeneProductDF( gset, geneMap=geneMap, geneMapColumn=givenGeneMapColumn)
		if ( nrow( gmap) < minGenesPerSet) next
		geneMapColumn <- "GENE_ID"

		# get the gene locations and folds for this set
		densityProfiles <- goodProfiles[[j]]
		for (k in 1:length(deList)) {
			tinyAt <- densityProfiles$at[[k]]
			tinyFold <- densityProfiles$fold[[k]]
			tinyGenes <- densityProfiles$gene[[k]]
			wh <- match( gmap[,1], tinyGenes)
			tiny <- data.frame( formatC( tinyFold[wh], format="f", digits=3, flag="+"),
						format( tinyAt[wh], format="d", big.mark=","),
						stringsAsFactors=F)
			colnames( tiny) <- gsub( "_", " ", paste( names(deList)[k], c("Fold", "Rank")))
			extra <- if ( k == 1) tiny else cbind( extra, tiny, stringsAsFactors=F)
		}
		gmap <- cbind( gmap, extra)
		f <- file.path( globalPlotPath, paste( "GeneSet_", j, ".html", sep=""))
		table2html( gmap, fileout=f, title=addSpeciesToHtmlTitle(cleanGeneSetName( names( geneSets)[j])), 
				linkColumnNames=c( geneMapColumn, "GroupsPerGene"),
				linkPaths=c( "../../pngPlots", "."),
				linkExtension=c(".png",".html"))

		thisSet <- rep( j, times=nrow(gmap))
		names( thisSet) <- gmap[[ geneMapColumn]]
		genesPathsSeen <- c( genesPathsSeen, thisSet)
		if ( j %% 10 == 0) cat( "\r", j, trimGeneSetNameLink( names(geneSets)[j]), "   ", length( geneSets[[j]]))
	   }

	   # now make the tables of all pathways for each gene
	   gpFac <- factor( names( genesPathsSeen))
	   cat( "\n\nBuilding gene membership tables...\n")
	   tapply( 1:length(genesPathsSeen), gpFac, function(x) {
	   		mygene <- names( genesPathsSeen)[ x[1]]
			mypathNums <- genesPathsSeen[x]
			# to save on tiny useless files, skip if too few groups for this gene
			if ( length( mypathNums) < minGenesPerSet) return()
			mypathNames <- names( geneSets)[ mypathNums]
			mypathNames <- sapply( mypathNames, cleanGeneSetName)
			mypathID <- paste( descriptor, mypathNums, sep="_")
			mygrpID <- paste( "GeneSet", mypathNums, sep="_")
			oneDF <- data.frame( "GeneSet"=mypathID, "Name"=mypathNames, 
					"GenesPerGroup"=mygrpID, stringsAsFactors=F)
			colnames(oneDF)[1] <- descriptor
			# order on the text of the pathway
			ord <- order( oneDF[ ,2])
			oneDF <- oneDF[ ord, ]
			rownames(oneDF) <- 1:nrow(oneDF)
			# grab these rows from the big table, and append the sample details
			who <- match( oneDF[ ,1], outForHTML[ , 1])
			extra <- outForHTML[ who, 5:ncol(outForHTML)]
			if (!doFDR) extra <- extra[ , -grep("FDR",colnames(extra))]
			colnames(extra) <- gsub( "_", " ", colnames(extra))
			colnames(extra) <- gsub( "\\.", " ", colnames(extra))
			oneDF <- cbind( oneDF, extra, stringsAsFactors=F)
			f <- paste( "GroupSet_", mygene, ".html", sep="")
			f <- file.cleanSpecialCharactersFromFileName( f)
			f <- file.path( globalPlotPath, f)
			table2html( oneDF, fileout=f, title=paste( mygene, gene2Product(mygene), sep=":  "),
				linkColumnNames=c(descriptor,"GenesPerGroup"), 
				linkPaths=rep.int(".",2), linkExtensions=c(".png",".html"))
			if ( x[1] %% 10 == 0) cat( "\r", x[1], mygene, nrow(oneDF))
		})

	   # we now have gene links that were never turned to plots so make those too.
	   #extraGenesToHTMLandPlots( genesToPlot=levels(gpFac))

	   # put back anything we overrode to true settings
	   geneMapColumn <- givenGeneMapColumn

	}

	if ( !makePlots) {
	   cat( "\n\nSkipping gene set density plots...\n")
	} else {
	   cat( "\n\nMaking gene set density plots...\n")
	   useYmin <- 1 / nrow( getCurrentGeneMap()) * 5
	   makeAllDensityPlots( geneSets, groupIDset=groupIDs, speciesID=speciesID, colorset=DE_colors,
	   		results.path=results.path, folderName=folderName, toolName=toolName,
	   		pngPath=globalPlotPath, pngName=descriptor, geneMapColumn=geneMapColumn, 
			whoToPlot=goodSets, deList=deList, yMin=useYmin, PLOT.FUN=PLOT.FUN,
			subsetInHTML=geneSetsToPlot)
	}

	# clean up any global storage... 
	if ( exists( "DE_List")) try( rm( DE_List, envir=.GlobalEnv))
	cat( "\n\nDone with Gene Sets for:  ", descriptor, "\n")

	return()
}


`sampleDensityPlot` <- function( Nrand=50, yScale=1.1, nLines=3 ){

	# measure all the density curves
	N <- 5600
	atBG <- seq(1,N,by=(N/Nrand))
	denBG <- density( atBG, n=1024,adjust=1, cut=3)
	randGbroad <- base::sort( round( rnorm( N, 2800, 1600))) [ seq( 5, N-5, by=N/Nrand) ]
	randGnarrow <- base::sort( round( rnorm( N, 2800, 700))) [ seq( 5, N-5, by=N/Nrand) ]
	atUP <- randGbroad - 1600 - 15
	atUP <- atUP[ atUP > 0]
	denUP <- density( atUP, n=1024, adjust=1, cut=3)
	atDN <- randGbroad + 1600 + 15
	atDN <- atDN[ atDN <= N]
	denDN <- density( atDN, n=1024, adjust=1, cut=3)
	atNO <- randGnarrow
	atNO <- atNO[ atNO <= N]
	atNO <- atNO[ atNO > 0]
	denNO <- density( atNO, n=1024, adjust=1, cut=3)

	maxY <- max( c( denUP$y, denNO$y, denDN$y)) * yScale
	if ( maxY < 0.0005) maxY <- 0.0005
	negY <- maxY * -0.1

	tickWidth <- 2
	if ( length( atBG) > 200) tickWidth <- 1
	if ( length( atBG) <= 60) tickWidth <- 3

	plot( denBG, xlim=c(1,N+1), ylim=c(negY, maxY), 
		main=paste( " N = ", Nrand, " genes"),
		xlab="Rank Position",
		type="l", col=1, lty=2, lwd=2, xaxs="r",
		cex.main=1.3, cex.axis=1.2, cex.lab=1.3, font.main=2, font.lab=2)
	lines( denNO, type="l", col=3, lwd=3)
	axis( 1, at=atNO, labels=NA, tcl=0.9, col=3, lwd=tickWidth)
	if (nLines > 1) {
	lines( denUP, type="l", col=2, lwd=3)
	axis( 1, at=atUP, labels=NA, tcl=0.9, col=2, lwd=tickWidth)
	}
	if (nLines > 2) {
	lines( denDN, type="l", col=4, lwd=3)
	axis( 1, at=atDN, labels=NA, tcl=0.9, col=4, lwd=tickWidth)
	}

	# lets re-draw the ticks in left to right order to better show the coloring
	if (nLines > 2) {
	allATs <- c( atUP, atDN, atNO)
	allCols <- c( rep.int(2,length(atUP)), rep.int(4,length(atDN)), rep.int(3,length(atNO))) 
	ord <- base::order( allATs)
	for( j in ord) axis( 1, at=allATs[j], labels=NA, tcl=0.9, col=allCols[j], lwd=tickWidth)
	}
	# black line at bottom
	axis( 1, at=c(0,N), labels=NA, tcl=0, col=1, lwd=tickWidth)

	legend( "topright", legend=c( "'No Net Change' Gene Set", "Up Regulated Gene Set", 
		"Down Regulated Gene Set", "uniform rank density")[c(1:nLines,4)],
		col=c(3,2,4,1)[c(1:nLines,4)], lty=c(1,1,1,2)[c(1:nLines,4)], lwd=3, bg="white", cex=1.4)
	
	text( x=c(500,2800,5000), y=negY, label=c( "Up-Regulated", "No Net Change", "Down-Regulated"), 
		cex=1.4, pos=3, font=2)

	return()
}


`readDEgroupsData` <- function( groupIDset, speciesID=getCurrentSpecies(), optionsFile="Options.txt",
					results.path=NULL, folderName="",
					toolName=c("MetaResults","RoundRobin","RankProduct","SAM","DESeq", "EdgeR")) {

	deList <- vector( mode="list")
	toolName <- match.arg( toolName)
	toolPrefix <- c( "RR", "RP", "SAM", "Meta", "DESeq", "EdgeR")[ match( toolName, 
				c("RoundRobin","RankProduct","SAM","MetaResults", "DESeq", "EdgeR"))]

	if ( getCurrentSpecies() != speciesID) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# set up directories to read from or write to...
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", notfound=".", verbose=FALSE)
	}
	if ( base::nchar(folderName) < 1) stop( "GeneSetDensity needs an explicit 'folderName' argument...")
	DE_path <- file.path( results.path, toolName, paste( prefix, folderName, sep="."))

	for( i in 1:length( groupIDset)) {

		# file name is made from sample / species / tool
		f <- paste( groupIDset[i], prefix, toolPrefix, "Ratio.txt", sep=".")
		if (toolPrefix == "Meta") f <- paste( groupIDset[i], prefix, toolPrefix, "UP.txt", sep=".")
		f <- file.path( DE_path, f)
		if ( ! file.exists(f) && toolPrefix == "Meta") {
			f <- paste( groupIDset[i], prefix, toolPrefix, "Ratio.txt", sep=".")
			f <- file.path( DE_path, f)
		}
		if ( ! file.exists(f)) {
			cat( "\n'readDEgroupsData()':  file not found: ", f)
			# allow older naming too
			if (toolPrefix == "Meta") f <- paste( groupIDset[i], prefix, toolPrefix, "Ratio.txt", sep=".")
			deList[[i]] <- data.frame()
			names( deList)[i] <- groupIDset[i]
		} else {
			deList[[i]] <- tbl <- read.delim( f, as.is=TRUE)
			names( deList)[i] <- groupIDset[i]
			expectedColumns <- c( "GENE_ID", "LOG2FOLD")
			if ( ! all( expectedColumns %in% colnames(tbl))) {
				cat( "\nRequired column names not found! \nFile: ",f, "\nExpected: ", 
						expectedColumns, "\n")
				stop()
			}
		}
	}
	return( deList)
}


`calcGeneSetDensity` <- function( deList, gSet, geneMapColumn="GENE_ID") {

	NS <- length( deList)
	atList <- densityList <- foldList <- geneList <- vector( mode="list", length=NS)
	bigYdensity <- 0
	geneMap <- getCurrentGeneMap()

	# inside the "DE_List", the gene name column is always called "GENE_ID"
	for( i in 1:NS) {

		# we expect perfect matches to GeneIDs..., but if the 'geneMapColumn' is
		# different, we need to search on those instead
		targetNames <- deList[[i]]$GENE_ID
		if ( geneMapColumn != "GENE_ID") {
			gmapRow <- base::match( targetNames, geneMap$GENE_ID, nomatch=0)
			targetNames[ gmapRow > 0] <- geneMap[[ geneMapColumn]][ gmapRow]
		}

		# for human genes, there may be non-unique names, so get all that are in 'geneSet'
		#where <- base::match( gSet, targetNames, nomatch=0)
		#atList[[i]] <- myAt <- where[ where > 0]
		where <- base::which( targetNames %in% gSet)
		atList[[i]] <- myAt <- where
		if ( length( myAt) > 1) {
			densityList[[i]] <- myDensity <- density( myAt, n=1024, adjust=1, cut=3)
			bigYdensity <- max( c( bigYdensity, myDensity$y))
		}
		
		folds <- deList[[i]]$LOG2FOLD
		foldList[[i]] <- folds[where]

		geneList[[i]] <- targetNames[where]
	}
	
	return( list( "who"=names(deList), "at"=atList, "density"=densityList, "bigY"=bigYdensity,
			"fold"=foldList, "gene"=geneList))
}


`measureDensityDifferences` <- function( densitySet, randomPvalues=NULL) {

	# turn the density curves from the DE comparisons of a gene subset into some 
	# numerical values we can use
	samples <- densitySet$who
	NS <- length( samples)
	at <- densitySet$at
	fold <- densitySet$fold
	
	out <- vector()
	nout <- 0
	bestPvalue <- 1
	bestRankShift <- 0

	nColumnPerSample <- 4
	for ( i in 1:NS) {
		us <- at[[i]]
		usFold <- fold[[i]]
		them <- vector()
		themFold <- vector()
		for( j in 1:NS) {
			if ( i == j) next
			them <- base::append( them, at[[j]])
			themFold <- base::append( themFold, fold[[j]])
		}

                if ( length( us) < 1 || length( them) < 1) {
                        deltaRank <- 0
                        pval <- 1
			netFold <- 0
                } else if ( length( us) == 1 || length(them) == 1) {
                        deltaRank <- mean.default(us) - mean.default( them)
                        pval <- 0.5
                        netFold <- mean.default(usFold) - mean.default(themFold)
                } else {
                        ans <- t.test( us, them)
                        deltaRank <- ans$estimate[2] - ans$estimate[1]
                        pval <- ans$p.value
                        netFold <- mean.default(usFold) - mean.default(themFold)
                }
		if ( ! is.null(randomPvalues)) {
			myFDR <- getFDRoneGeneSet( length(us), pval, randomPvalues)
		} else {
			myFDR <- NA
		}
		out[ nout + 1] <- netFold
		out[ nout + 2] <- pval
		out[ nout + 3] <- myFDR
		out[ nout + 4] <- deltaRank
		nout <- nout + nColumnPerSample
		if( pval < bestPvalue) bestPvalue <- pval
		if( abs(deltaRank) > bestRankShift) bestRankShift <- abs(deltaRank)
	}
	names( out) <- paste( rep( samples, each=nColumnPerSample), 
				c( "FoldChange", "P.value", "FDR", "RankShift"), sep=" ")
	return( list( "bestPvalue"=bestPvalue, "bestRankShift"=bestRankShift, "measurements"=out))
}



`densityMetrics` <- function( d) {

	pointMasses <- d$x * d$y
	totalMass <- sum( pointMasses)
	totalDensity <- sum( d$y)
	centerOfMass <- totalMass / totalDensity
	pointDev <- d$y * abs( d$x - centerOfMass)
	absDev <- sum( pointDev)
	stdDev <- absDev / totalDensity

	return( list( "center"=centerOfMass, "sd"=stdDev))
}


`geneSetsDensityPlot` <- function( gSet, groupIDset, label="", reload=TRUE, speciesID=getCurrentSpecies(), 
				results.path=NULL, folderName="", colorset=c( 2:(length(groupIDset)+1)), 
				toolName=c("RoundRobin","RankProduct","SAM","MetaResults"),
				geneMapColumn="GENE_ID", deList=NULL, yScale=1.1, yMin=0.0005) {

	NS <- length( groupIDset)

	if (reload) {
	    if ( is.null( deList)) {
		toolName <- match.arg( toolName)
		DE_List <<- readDEgroupsData( groupIDset, speciesID=speciesID, results.path=results.path, 
					folderName=folderName, toolName=toolName)
		if ( length( DE_List) != length( groupIDset)) stop( "failed reading GeneSet DE data")
	    } else {
	    	DE_List <<- deList
	    }
	}
	NG <- nrow( DE_List[[1]] )

	# measure all the density curves for this gene set, and note the biggest Y value
	ans <- calcGeneSetDensity( DE_List, gSet, geneMapColumn)
	atList <- ans$at
	densityList <- ans$density
	bigYdensity <- ans$bigY

	# scaling and display setup
	maxY <- bigYdensity * yScale
	if ( maxY < yMin) maxY <- yMin
	negY <- maxY * -0.1
	tickWidth <- 2
	if ( length( gSet) > 200) tickWidth <- 1
	if ( length( gSet) <= 60) tickWidth <- 3
	colorList <- colorset

	# draw a background density
	denBG <- density( seq(1,NG,by=10), n=1024, adjust=1, cut=3)
	plot( denBG, xlim=c(1,NG+1), ylim=c(negY, maxY), 
		main=paste( label, "\n N = ",length(gSet), " genes"),
		xlab="Rank Position",
		type="l", col=1, lty=2, lwd=2, xaxs="r",
		cex.main=1.3, cex.axis=1.2, cex.lab=1.3, font.main=2, font.lab=2, font.axis=2)
	
	# add lines for each sample
	for( i in 1:NS) {
		lines( densityList[[i]], type="l", col=colorList[i], lwd=3)
		axis( 1, at=atList[[i]], labels=NA, tcl=0.9, col=colorList[i], lwd=tickWidth)
	}

	# lets re-draw the ticks in left to right order to better show the coloring
	allATs <- allCols <- vector()
	for( i in 1:NS) {
		allATs <- base::append( allATs, atList[[i]])
		allCols <- base::append( allCols, rep.int( colorList[i], times=length(atList[[i]])))
	}
	ord <- base::order( allATs)
	for( j in ord) axis( 1, at=allATs[j], labels=NA, tcl=0.9, col=allCols[j], lwd=tickWidth)
	# black line at bottom
	axis( 1, at=c(0,NG), labels=NA, tcl=0, col=1, lwd=tickWidth)

	legend( "topright", legend=c( groupIDset, "uniform density"),
		col=c( colorList, 1), lty=c( rep(1,NS), 2), lwd=3, bg="white", cex=1.2)
	
	text( x=c( (NG*0.1), (NG*0.5), (NG*0.9)), y=negY, label=c( "Up-Regulated", "No Net Change", 
		"Down-Regulated"), cex=1.2, font=2, pos=3)

	return()
}


`makeAllDensityPlots` <- function(  allGeneSets, groupIDset, speciesID=getCurrentSpecies(), 
				colorset=c(2:(length(groupIDset)+1)), results.path=results.path, 
				folderName=folderName, toolName=c("RoundRobin","RankProduct","SAM","MetaResults"),
				pngPath="densityPlots", pngName="Pathway", geneMapColumn="GENE_ID", 
				whoToPlot=1:length(allGeneSets), deList=NULL, yMin=0.0005, PLOT.FUN=NULL,
				subsetInHTML=NULL) {

	pathnames <- names( allGeneSets)
	ngenes <- sapply( allGeneSets, length)
	firstGood <- which( ngenes > 1)[1]

	# plot once to get window set up and load the data
	toolName <- match.arg( toolName)
	pngFile <- file.path( pngPath, paste( pngName, "_", firstGood, ".png", sep=""))
	png( filename=pngFile, width=1000, height=700, bg="white")

	geneSetsDensityPlot( gSet=allGeneSets[[ firstGood]], groupIDset, 
				label=paste( pngName, ":   ", trimGeneSetNameLink( pathnames[ firstGood])), 
				results.path=results.path, folderName=folderName, toolName=toolName,
				reload=TRUE, speciesID=speciesID, colorset=colorset, 
				geneMapColumn=geneMapColumn, deList=deList, yMin=yMin)

	dev.off()

	# now do all that have genes...

	# we may have a subset of gene set numbers to be drawn.. that happen to be character format
	if ( ! is.null( subsetInHTML)) {
		subsetInHTML <- as.integer( subsetInHTML)
	}

	# for( j in 1:length( pathnames)) {
	for( j in whoToPlot) {
		if ( ngenes[j] < 2) next

		# do we skip this one for not being in any HTML file?
		if ( (! is.null( subsetInHTML)) && ( !(j %in% subsetInHTML))) next

		pngFile <- file.path( pngPath, paste( pngName, "_", j, ".png", sep=""))
		png( filename=pngFile, width=1000, height=700, bg="white")

		if ( is.null( PLOT.FUN)) {
			geneSetsDensityPlot( gSet=allGeneSets[[j]], groupIDset, 
				label=paste( pngName, ":   ", trimGeneSetNameLink( cleanGeneSetName( pathnames[j]))), 
				reload=FALSE,
				results.path=results.path, folderName=folderName, toolName=toolName,
				speciesID=speciesID, colorset=colorset, geneMapColumn=geneMapColumn, 
				yMin=yMin)
		} else {
			label <- paste( "GeneSet:           N_Genes=", length(allGeneSets[[j]]), "\n",
					trimGeneSetNameLink( cleanGeneSetName( pathnames[j])))
			PLOT.FUN( allGeneSets[[j]], label)
		}

		dev.off()

		if ( j %% 10 == 0) cat( "\r", j, trimGeneSetNameLink( pathnames[j]), "  N_genes:", ngenes[j], "    ")
	}

	return()
}


calcRandomGeneSetPvalues <- function( deList, geneSets, who=1:length(geneSets), nSimulations=1000) {

	# given an actual RR_List of DE gene rankings, let's estimate some
	# P-values from multiple random pertubations...
	cat( "\nPre-estimating FDR for all gene sets...\n")
	NS <- length( deList)
	allLists <- 1:NS
	atList <- vector( mode="list", length=NS)
	targetNames <- deList[[1]]$GENE_ID
	setLengths <- sort( unique( sapply( geneSets, function(x) {
				if (is.null(x)) return(0)
				return( length(x))
			})))
	setLengths <- setLengths[ setLengths > 1]
	NSETS <- length( setLengths)
	allSets <- 1:NSETS
	pvalueList <- vector( mode="list", length=NSETS)

	# repeat till every possible number of genes per set has enough simulations
	nSimulations <- ceiling( nSimulations / NS)
	lapply( 1:nSimulations, FUN=function(iSimu) {
		lapply( allSets, function(iSet) {
			thisSize <- setLengths[iSet]
			thisGenes <- sample( targetNames, size=thisSize, replace=(thisSize > length(targetNames)))
			lapply( allLists, function(k) { atList[[k]] <<- base::which( deList[[k]]$GENE_ID %in% thisGenes)})
			newPs <- vector( length=NS)
			lapply( allLists, function(k) {
				us <- atList[[k]]
				them <- unlist( atList[ setdiff( allLists, k)])
                        	ans <- t.test( us, them)
                        	newPs[k] <<- ans$p.value
			})
                        pvalueList[[iSet]] <<- c( pvalueList[[iSet]], newPs)
		})
		cat( "\rIter:", iSimu)
	})

	# OK, now every length of gene set has at least 'nSimu' P-values
	for ( iSet in 1:length(setLengths)) pvalueList[[iSet]] <- sort( pvalueList[[iSet]])
	names( pvalueList) <- setLengths

	return( pvalueList)
}


getFDRoneGeneSet <- function( nGenes, observedPvalue, randomPvalues) {

	# given a calculated P-value for a given size of geneSet, how likely
	if (is.null(randomPvalues)) return(NA)

	# was getting that good a P by chance
	randomSetSizes <- as.integer( names( randomPvalues))
	myPtr <- findInterval( nGenes, randomSetSizes)
	if ( myPtr < 1) myPtr <- 1

	thesePs <- randomPvalues[[myPtr]]
	nBetter <- sum( thesePs < observedPvalue)
	fdr <- nBetter / length(thesePs)
	return( fdr)
}


`trimGeneSetNameLink` <- function( text) {

	# strip any hyper link out of the the name term
	return( sub( " *<a .+", "", text))
}


`makeSmallGeneProductDF` <- function( gset, geneMap, geneMapColumn) {

	NG <- length(gset)
	gprods <- rep.int( "", NG)

	geneList <- geneMap[[ geneMapColumn]]
	where <- match( gset, geneList, nomatch=0)
	gprods[ where > 0] <- geneMap$PRODUCT[ where]

	needProd <- which( gprods == "")
	if ( length( needProd)) {
		gprods[ needProd] <- gene2ProductAllSpecies( gset[ needProd])
	}

	gmap <- data.frame( "GENE_ID"=gset, "PRODUCT"=gprods, stringsAsFactors=F)
	gmap$GroupsPerGene <- paste( "GroupSet", gset, sep="_")

	ord <- order( gmap$GENE_ID)
	gmap <- gmap[ ord, ]
	rownames(gmap) <- 1:nrow(gmap)
	return( gmap)
}
