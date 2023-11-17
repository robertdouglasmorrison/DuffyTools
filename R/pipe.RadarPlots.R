# pipe.RadarPlots.R -- visualize gene DE as a family of radar (spider) plots

`pipe.RadarPlots` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, folderName=NULL,
				tool=c( "MetaResults", "DESeq", "EdgeR", "RankProduct", "RoundRobin", "SAM"),
				groupColumn="Group", colorColumn="Color", average.FUN=median,
				geneSets=defaultGeneSets(speciesID), restrictionSets=NULL, baselineGroup=NULL,
				legend.prefix=NULL, legend.order=NULL, legend.cex=1, Nshow=24, cex=0.80,
				start=pi/4, radial.labels=FALSE, radial.margin=c( 2,2,6,2),
				radial.lim=NULL, min.radial.lim=NULL, boxed.radial=F, label.prop=1, lwd=5, 
				main=paste( "Comparison:  ",folderName), ...)
{

	require( plotrix)
	checkX11( bg='white', width=9, height=7)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# get the SampleIDs that we know are valid
	annT <- readAnnotationTable( annotationFile)
	allSamples <- sampleIDset
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allSamples <- myAnnT$SampleID
	if ( ! (groupColumn %in% colnames(myAnnT))) stop( paste( "Given grouping column not in annotation table: ", groupColumn))
	if ( ! (colorColumn %in% colnames(myAnnT))) stop( paste( "Given coloring column not in annotation table: ", colorColumn))

	# build the path to the DE results we will use, and put our Radar results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", verbose=F)
	}
	if ( is.null( folderName)) stop( "Explicit Folder Name of DE results is required")
	folderName <- paste( prefix, folderName, sep=".")
	tool <- match.arg( tool)
	dePath <- file.path( results.path, tool, folderName)
	if ( ! file.exists( dePath)) stop( paste( "DE results folder not found: ", dePath))
	radarPath <- file.path( dePath, "RadarPlots")
	if ( ! file.exists( radarPath)) dir.create( radarPath, recursive=T)

	# we have to read the transcriptomes in... to pre-load the data
	needLoad <- TRUE
	if (needLoad) {
		# still to do... get 'multi' vs 'unique' call?...
		intensityColumn <- getExpressionUnitsColumn( optionsFile, useMultiHits=TRUE)
		cat( "\nLoading ", length( allSamples), "transcriptomes..")
		files <- paste( allSamples, prefix, "Transcript.txt", sep=".")
		files <- file.path( results.path, "transcript", files)
		radarM <<- expressionFileSetToMatrix( files, allSamples, intensityColumn=intensityColumn, verbose=T)
		cat( "\nConverting Expression Abundance to M-values..")
		radarMA <<- expressionMatrixToMvalue( radarM, average.FUN=average.FUN)
		cat( "  Done.\n")
	}


	# now we are read to make those radar plots
	if ( is.list( geneSets)) {
		ans <- pipe.OneRadarPlot( sampleIDset=allSamples, speciesID=speciesID, annotationFile=annotationFile,
				optionsFile=optionsFile, results.path=results.path, groupColumn=groupColumn,
				colorColumn=colorColumn, geneSet=geneSets, restrictionSets=restrictionSets,
				baselineGroup=baselineGroup, reload=FALSE, Nshow=Nshow, start=start, 
				radial.labels=radial.labels, radial.margin=radial.margin, radial.lim=radial.lim, min.radial.lim=min.radial.lim,
				boxed.radial=boxed.radial, label.prop=label.prop, lwd=lwd, main=main, 
				legend.order=legend.order, legend.prefix=legend.prefix, legend.cex=legend.cex, cex=cex, ...)
		geneSetName <- "Radar"
		plotFile <- file.path( radarPath, geneSetName)
		printPlot( plotFile, optT=optionsFile)
		csvFile <- file.path( radarPath, paste( geneSetName, "csv", sep="."))
		write.table( ans, csvFile, sep=",", quote=T, row.names=F)

		nOut <- min( Nshow, nrow(ans))
		# if the full set is just a bit larger, go ahead and keep them all?...
		if ( Nshow * 1.25 > nrow(ans)) nOut <- nrow(ans)
		out <- ans[ 1:nOut, ]

	} else {
		out <- vector( mode="list")
		for (k in 1:length(geneSets)) {
			gs <- geneSets[k]
			cat( "\nDoing Radar Plot for:  ", gs)

			# a few geneSets was "()" for several things...
			wrapParentheses <- ( length( grep( "BTM", gs)) == 0)

			thisRestrictionSet <- restrictionSets
			if ( ! is.null( restrictionSets)) {
				where <- match( gs, names(restrictionSets), nomatch=0)
				if ( where > 0) {
					thisRestrictionSet <- restrictionSets[[ where]]
					cat( "\n  Using a restriction set of pathways..")
				}
			}

			ans <- pipe.OneRadarPlot( sampleIDset=allSamples, speciesID=speciesID, annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path, groupColumn=groupColumn,
					colorColumn=colorColumn, geneSet=gs, restrictionSets=thisRestrictionSet,
					baselineGroup=baselineGroup, reload=FALSE, Nshow=Nshow, start=start, 
					radial.labels=radial.labels, radial.margin=radial.margin, radial.lim=radial.lim, min.radial.lim=min.radial.lim,
					boxed.radial=boxed.radial, label.prop=label.prop, lwd=lwd, main=main, 
					legend.order=legend.order, legend.prefix=legend.prefix, legend.cex=legend.cex, cex=cex, 
					wrapParentheses=wrapParentheses, ...)
			if ( ! nrow( ans)) {
				out[[k]] <- ans
				next
			}

			plotFile <- file.path( radarPath, paste( "Radar", gs, sep="."))
			printPlot( plotFile, optT=optionsFile)
			csvFile <- file.path( radarPath, paste( "Radar", gs, "csv", sep="."))
			write.table( ans, csvFile, sep=",", quote=T, row.names=F)

			nOut <- min( Nshow, nrow(ans))
			# if the full set is just a bit larger, go ahead and keep them all?...
			if ( Nshow * 1.25 > nrow(ans)) nOut <- nrow(ans)
			smallOut <- ans[ 1:nOut, ]
			out[[k]] <- smallOut
		}
		names( out) <- geneSets
		cat( "\nDone Radar Plots for", length(geneSets), "GeneSets.\n")
	}

	# make a single combined .CSV file of all the gene sets together
	combineRadarPlotResults( radarPath)

	return( invisible( out))
}


`pipe.OneRadarPlot` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  
				groupColumn="Group", colorColumn="Color", average.FUN=median,
				legend.prefix=NULL, legend.order=NULL, legend.cex=1.0, cex=par("cex.axis"), 
				geneSet=defaultGeneSets(speciesID), restrictionSets=NULL, baselineGroup=NULL,
				reload=FALSE, Nshow=24, 
				start=pi/4, radial.labels=FALSE, radial.margin=c( 2,2,6,2),
				radial.lim=NULL, min.radial.lim=NULL, boxed.radial=F, label.prop=1, lwd=5, main=NULL, 
				wrapParentheses=TRUE, ...)
{

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	mainText <- "Radar Plot:  "
	if ( is.character( geneSet)) {
		if ( missing(geneSet)) geneSet <- match.arg( geneSet)
		geneSetFile <- paste( prefix, geneSet, sep=".")
		allGeneSets <- NULL
		data( list=geneSetFile, envir=environment())
		if ( is.null( allGeneSets)) {
			cat( "\nFailed to find dataset: ", geneSetFile, "   Skipping..")
			return(data.frame())
		}
		mainText <- paste( mainText, "    ", geneSet)
		if ( ! is.null( main)) mainText <- paste( mainText, main, sep="\n")
	} else {
		if ( typeof( geneSet) != "list") stop( "'geneSet' must be a character string or a list")
		allGeneSets <- geneSet
		if ( ! is.null( main)) mainText <- paste( mainText, main, sep="    ")
	}

	# use the samples and grouping column to know the colors and files to load
	annT <- readAnnotationTable( annotationFile)
	allSamples <- sampleIDset
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allSamples <- myAnnT$SampleID

	if ( ! (groupColumn %in% colnames(myAnnT))) stop( paste( "Given grouping column not in annotation table: ", groupColumn))
	grpFac <- factor( myAnnT[[ groupColumn]])
	allGroups <- tapply( myAnnT$SampleID, grpFac, FUN=c, simplify=FALSE)
	NperGroup <- sapply( allGroups, length)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", verbose=F)
	}

	# map from the group names to find the colors to use...
	if ( ! (colorColumn %in% colnames(myAnnT))) stop( paste( "Given coloring column not in annotation table: ", colorColumn))
	where <- base::match( names(allGroups), myAnnT[[ groupColumn]])
	mycolors <- myAnnT[[ colorColumn]][ where]

	# see if we have to read the transcriptomes in...
	needLoad <- TRUE
	if ( exists( "radarM") && all( colnames(radarM) == allSamples) && nrow(radarM) > 100 && !reload) needLoad <- FALSE
	if (needLoad) {
		files <- paste( allSamples, prefix, "Transcript.txt", sep=".")
		files <- file.path( results.path, "transcript", files)
		radarM <<- expressionFileSetToMatrix( files, allSamples, verbose=T)
		radarMA <<- expressionMatrixToMvalue( radarM, average.FUN=average.FUN)
	}

	# try to standarize the names of the gene sets to keep them easy to view
	fullPathNames <- names(allGeneSets)
	names(allGeneSets) <- cleanGeneSetModuleNames( names(allGeneSets), wrapParentheses=wrapParentheses)

	# we may have been given a restriction set, to only consider gene sets from some previous run
	if ( ! is.null( restrictionSets)) {
		pathNames <- restrictionSets$PathName
		if ( ! is.null( pathNames)) {
			pathNames <- cleanGeneSetModuleNames( pathNames, wrapParentheses=wrapParentheses)
			keepers <- which( names( allGeneSets) %in% pathNames)
			#cat( "\nDebug: restrictionSet: ", head( pathNames))
			#cat( "\nDebug: Keepers: ", head( keepers))
			if ( length( keepers)) {
				allGeneSets <- allGeneSets[ keepers]
				fullPathNames <- fullPathNames[ keepers]
			}
		}
	}

	# do the reduction & grouping.  This call must use mean average, never pass down the method used for the MA step.
	radarAns <- reduceMatrixToModules( radarMA, geneModules=allGeneSets, sampleTraits=allGroups,
				gene.names=shortGeneName( rownames( radarMA), keep=1), average.FUN=mean,
				sample.names=colnames(radarMA), baselineTrait=baselineGroup)
	mShow <- radarMOD <- radarAns$matrix
	pShow <- radarPvalue <- radarAns$p.value
	piShow <- radarPIvalue <- radarAns$pi.value
	validModuleNames <- radarAns$moduleNames
	nGenes <- radarAns$geneCounts

	# the reduction may have thrown out some modules
	keep <- which( names(allGeneSets) %in% validModuleNames)
	fullPathNames <- fullPathNames[keep]

	# trim down to less groups if we need to
	Nshow <- min( Nshow, nrow( radarMOD))
	# if the full set is just a bit larger, go ahead and keep them all?...
	if ( Nshow * 1.25 > nrow(radarMOD)) Nshow <- nrow(radarMOD)

	# use both Pvalue and magnitudes to decide who to draw
	# but give magnitudes more weight?...
	magnitudes <- diff( apply( radarMOD, 1, range))
	bestPs <- apply( radarPvalue, 1, min)
	ord <- diffExpressRankOrder( magnitudes, bestPs, wt.fold=2, wt.pvalue=1)

	# try to not let one single group have all the UP spokes...
	nUpSpokes <- apply( radarMOD[ ord[1:Nshow], ], MARGIN=2, function(x) sum(x > 0, na.rm=T))
	whoMostUp <- which.max( nUpSpokes)
	idealMaxUp <- round( Nshow * 0.8)
	if ( nUpSpokes[whoMostUp] > idealMaxUp) {
		# juggle a few elements of 'ord', to put some 'down' rows into the top N
		nDownWanted <- nUpSpokes[whoMostUp] - idealMaxUp
		downSpokeRows <- which( radarMOD[ ord, whoMostUp] < 0)
		# in case there are some top N in the downs, don't let them come in a second time
		downSpokeRows <- setdiff( downSpokeRows, 1:Nshow)
		# keep exchanging the least Up for the most Down, until we have enough
		tryRowNow <- Nshow
		nDownUsed <- 0
		while (nDownWanted > 0) {
			if (tryRowNow <= 1) break
			# don't swap out one that is already down
			if (radarMOD[ ord[tryRowNow], whoMostUp] < 0) {
				tryRowNow <- tryRowNow - 1
				next
			}
			# at the lowest up row, swap it out
			nDownUsed <- nDownUsed + 1
			if (nDownUsed > length(downSpokeRows)) break
			bestNegRow <- downSpokeRows[ nDownUsed]
			tmpRow <- ord[ tryRowNow]
			ord[ tryRowNow] <- ord[ bestNegRow]
			ord[ bestNegRow] <- tmpRow
			# successful, so update
			nDownWanted <- nDownWanted - 1
			tryRowNow <- tryRowNow - 1
		}
	}

	# also adding in piValues too, as a third score/rank metric 
	bestPIvals <- piValue( magnitudes, bestPs)
	mShow <- radarMOD[ ord[1:Nshow], ]
	pShow <- radarPvalue[ ord[1:Nshow], ]
	piShow <- radarPIvalue[ ord[1:Nshow], ]
	mOut <- radarMOD[ ord, ]
	pOut <- radarPvalue[ ord, ]
	piOut <- radarPIvalue[ ord, ]
	magniOut <- magnitudes[ ord]
	bestPout <- bestPs[ ord]
	bestPIout <- bestPIvals[ ord]
	validModuleNames <- validModuleNames[ ord[ 1:Nshow]]
	fullPathNames <- fullPathNames[ ord]
	nGenes <- nGenes[ ord]

	# now alphabetize on any part after the Module names
	ord <- order( sub( "M.+: +", "", validModuleNames), sub( "(^M[0-9]+\\.[0-9]+:)(.+)", "\\1", validModuleNames))
	mShow <- mShow[ ord, ]
	pShow <- pShow[ ord, ]
	piShow <- piShow[ ord, ]
	validModuleNames <- validModuleNames[ ord]

	# plot it now
	require( plotrix)
	# either set the limits from the data..
	if ( is.null( radial.lim)) {
		radial.lim <- range( as.vector( mShow)) * 1.5
		if ( radial.lim[2] < 0.1) radial.lim[2] <- 0.1
		# make sure we span the zero point, so any 'baseline' group is sure to be drawn
		negThreshold <- max( as.vector( mShow), na.rm=T) * -0.3
		if ( radial.lim[1] > negThreshold) radial.lim[1] <- negThreshold
		if ( radial.lim[1] > -0.1) radial.lim[1] <- -0.1
	} else {
		# or clip the data to those limits
		mShow <- pmax( mShow, (radial.lim[1]*0.95))
		mShow <- pmin( mShow, (radial.lim[2]*0.95))
	}
	# perhaps allow a minimum limits
	if ( ! is.null( min.radial.lim)) {
		if ( min.radial.lim[1] < radial.lim[1]) radial.lim[1] <- min.radial.lim[1]
		if ( min.radial.lim[2] > radial.lim[2]) radial.lim[2] <- min.radial.lim[2]
	}
	DuffyTools::radial.plot( t(mShow), labels=validModuleNames, radlab=radial.labels, rp.type="p", line.col=mycolors,
			start=start, clockwise=T, mar=radial.margin, radial.lim=radial.lim, label.prop=label.prop,
			show.grid.labels=3, lwd=lwd, main=mainText, cex=cex, ...)

	# take more control of the legend location
	usr <- par( "usr")
	par( 'xpd'=NA)
	legendText <- names(allGroups)
	if ( any( NperGroup > 1)) legendText <- paste( names(allGroups), "  (N=",NperGroup, ")", sep="")
	if ( ! is.null(legend.prefix)) legendText <- paste( legend.prefix, legendText, sep=":  ")
	if ( ! is.null(legend.order)) {
		legendText <- legendText[ legend.order]
		mycolors <- mycolors[ legend.order]
	}
	legend.title <- NULL
	if ( ! is.null( baselineGroup)) {
		legend.title <- paste( "Baseline = '", baselineGroup, "'", sep="")
		legend.cex <- legend.cex * 0.9
	}
	legend( x=usr[1]*1.7, y=usr[4]*1.125, legendText, lwd=lwd, col=mycolors, bg="white", cex=legend.cex, 
		title=legend.title)
	par( 'xpd'=FALSE)
	dev.flush()

	# return the full table, what we drew is at the very top.  Try to smartly truncate the values
	if ( max(magniOut,na.rm=T) > 1) {
		magniOut <- round( magniOut, digits=3)
		mOut <- round( mOut, digits=3) 
	} else if ( max(magniOut,na.rm=T) > 0.01) {
		magniOut <- round( magniOut, digits=5)
		mOut <- round( mOut, digits=5) 
	}

	# let's not keep the HTML links...
	fullPathNames <- cleanGeneSetModuleNames( fullPathNames, wrap=F)
	out <- data.frame( "PathName"=fullPathNames, "N_Genes"=nGenes, "DeltaFold"=magniOut,
				"BestPvalue"=formatC( bestPout, format="e", digits=2), 
				"BestPIvalue"=round( bestPIout, digits=3), mOut,
				formatC( pOut, format="e", digits=2), 
				round( piOut, digits=3), stringsAsFactors=F)
	colnames(out) <- c( "PathName", "N_Genes", "DeltaFold", "BestPvalue", "BestPIvalue",
				paste( "Fold", colnames(mOut), sep="_"), paste( "Pvalue", colnames(mOut), sep="_"),
				paste( "PIvalue", colnames(mOut), sep="_"))
	rownames(out) <- 1:nrow(out)

	# reorder to get the 2 (now 3!) values per group together
	#outOrd <- c( 1:4)
	#for ( k in 1:ncol(mOut)) outOrd <- c( outOrd, (k+4), (ncol(mOut)+k+4))
	outOrd <- c( 1:5)
	for ( k in 1:ncol(mOut)) outOrd <- c( outOrd, (k+5), (ncol(mOut)+k+5), (ncol(mOut)*2+k+5))
	out <- out[ , outOrd]

	return(out)
}


`radarPlot` <- function( m, row.names=rownames(m), colGroups=colnames(m), colors=1:ncol(m),
				legend.prefix=NULL, legend.order=NULL, foldChangeTransform=TRUE,
				average.FUN=median,
				geneSet=defaultGeneSets(), restrictionSets=NULL, baselineGroup=NULL,
				Nshow=24, start=pi/4, radial.labels=FALSE, radial.margin=c( 2,2,6,2),
				radial.lim=NULL, min.radial.lim=NULL, boxed.radial=F, label.prop=1, lwd=5, main=NULL, 
				legend.cex=1.1, cex=par("cex.axis"), ...)
{

	prefix <- getCurrentSpeciesFilePrefix()

	if ( ! is.null( geneSet)) {
		if ( is.character( geneSet)) {
			geneSet <- match.arg( geneSet)
			geneSetFile <- paste( prefix, geneSet, sep=".")
			allGeneSets <- NULL
			data( list=geneSetFile, envir=environment())
			if ( is.null( allGeneSets)) stop( paste( "Failed to load dataset: ", geneSetFile))
			if (is.null( main)) main <- paste( "Radar Plot:  ", geneSet)
		} else {
			if ( typeof( geneSet) != "list") stop( "'geneSet' must be a character string or a list")
			allGeneSets <- geneSet
			if (is.null( main)) main <- "Radar Plot:  "
		}
	}

	allGroups <- tapply( colnames(m), factor(colGroups), c)
	NperGroup <- sapply( allGroups, length)
	mycolors <- colors[1:length(allGroups)]

	if ( foldChangeTransform) {
		radarM <<- m
		radarMA <<- expressionMatrixToMvalue( radarM, average.FUN=average.FUN)
	} else {
		radarMA <<- m
	}

	# try to standarize the names of the gene sets to keep them easy to view
	if ( ! is.null( geneSet)) {
		fullPathNames <- names(allGeneSets)
		names(allGeneSets) <- cleanGeneSetModuleNames( names(allGeneSets))
	}

	# we may have been given a restriction set, to only consider gene sets from some previous run
	if ( ! is.null( restrictionSets)) {
		pathNames <- restrictionSets$PathName
		if ( ! is.null( pathNames)) {
			pathNames <- cleanGeneSetModuleNames( pathNames)
			keepers <- which( names( allGeneSets) %in% pathNames)
			if ( length( keepers)) {
				allGeneSets <- allGeneSets[ keepers]
				fullPathNames <- fullPathNames[ keepers]
			}
		}
	}

	# do the reduction & grouping
	if ( ! is.null( geneSet)) {
		# do the reduction & grouping.  This call must use mean average, never pass down the method used for the MA step.
		radarAns <- reduceMatrixToModules( radarMA, geneModules=allGeneSets, sampleTraits=allGroups,
						gene.names=shortGeneName(row.names,keep=1), baselineTrait=baselineGroup,
						average.FUN=mean)
		mShow <- radarMOD <- radarAns$matrix
		pShow <- radarPvalue <- radarAns$p.value
		validModuleNames <- radarAns$moduleNames
		nGenes <- radarAns$geneCounts

		# the reduction may have thrown out some modules
		keep <- which( names(allGeneSets) %in% validModuleNames)
		fullPathNames <- fullPathNames[keep]

		# trim down to less groups if we need to
		Nshow <- min( Nshow, nrow( radarMOD))
		# if the full set is just a bit larger, go ahead and keep them all?...
		if ( Nshow * 1.25 > nrow(radarMOD)) Nshow <- nrow(radarMOD)

		# use both Pvalue and magnitudes
		magnitudes <- diff( apply( radarMOD, 1, range))
		bestPs <- apply( radarPvalue, 1, min)
		ord <- diffExpressRankOrder( magnitudes, bestPs, wt.fold=1, wt.pvalue=1)
		mShow <- radarMOD[ ord[1:Nshow], ]
		pShow <- radarPvalue[ ord[1:Nshow], ]
		mOut <- radarMOD[ ord, ]
		pOut <- radarPvalue[ ord, ]
		nGenes <- nGenes[ ord]
		validModuleNames <- validModuleNames[ ord[ 1:Nshow]]
		fullPathNames <- fullPathNames[ ord]

		# alphabetize on any part after the Module names
		ord <- order( sub( "M.+: +", "", validModuleNames), sub( "(^M[0-9]+\\.[0-9]+:)(.+)", "\\1", validModuleNames))
		mShow <- mShow[ ord, ]
		pShow <- pShow[ ord, ]
		validModuleNames <- validModuleNames[ ord]
	} else {
		if ( Nshow < nrow(radarMA)) {
			mRange <- t( apply( radarMA, 1, range))
			ord <- order( apply( mRange, 1, diff), decreasing=T)
			keep <- sort( ord[ 1:Nshow])
			radarMA <- radarMA[ keep, ]
		}
		mShow <- mOut <- radarMA
		fullPathNames <- validModuleNames <- rownames(radarMA)
		allGroups <- colnames(radarMA)
		names(allGroups) <- allGroups
		NperGroup <- rep.int( 1, length(allGroups))
		pOut <- matrix( 1, nrow=nrow(mOut), ncol=ncol(mOut))
		nGenes <- rep.int( NA, nrow(mOut))
	}

	# plot it now
	require( plotrix)
	# either set the limits from the data..
	if ( is.null( radial.lim)) {
		radial.lim <- range( as.vector( mShow)) * 1.5
		# make sure we span the zero point, so any 'baseline' group is sure to be drawn
		negThreshold <- max( as.vector( mShow), na.rm=T) * -0.3
		if ( radial.lim[1] > negThreshold) radial.lim[1] <- negThreshold
	} else {
		# or clip the data to those limits
		mShow <- pmax( mShow, (radial.lim[1]*0.95))
		mShow <- pmin( mShow, (radial.lim[2]*0.95))
	}
	# perhaps allow a minimum limits
	if ( ! is.null( min.radial.lim)) {
		if ( min.radial.lim[1] < radial.lim[1]) radial.lim[1] <- min.radial.lim[1]
		if ( min.radial.lim[2] > radial.lim[2]) radial.lim[2] <- min.radial.lim[2]
	}
	DuffyTools::radial.plot( t(mShow), labels=validModuleNames, radlab=radial.labels, rp.type="p", line.col=mycolors,
			start=start, clockwise=T, mar=radial.margin, radial.lim=radial.lim, label.prop=label.prop,
			show.grid.labels=3, lwd=lwd, main=main, cex=cex, ...)

	# take more control of the legend location
	usr <- par( "usr")
	par( 'xpd'=NA)
	legendText <- names(allGroups)
	if ( any( NperGroup > 1)) legendText <- paste( names(allGroups), "  (N=",NperGroup, ")", sep="")
	if ( ! is.null(legend.prefix)) legendText <- paste( legend.prefix, legendText, sep=":  ")
	if ( ! is.null(legend.order)) {
		legendText <- legendText[ legend.order]
		mycolors <- mycolors[ legend.order]
	}
	legend.title <- NULL
	if ( ! is.null( baselineGroup)) {
		legend.title <- paste( "Baseline = '", baselineGroup, "'", sep="")
		legend.cex <- legend.cex * 0.9
	}
	legend( x=usr[1]*1.5, y=usr[4]*1.1, legendText, lwd=lwd, col=mycolors, bg="white", cex=legend.cex, 
		title=legend.title)
	par( 'xpd'=FALSE)
	dev.flush()

	# return the full table, what we drew is at the very top
	out <- data.frame( "PathName"=fullPathNames, "N_Genes"=nGenes, round(mOut,digits=4), 
				formatC( pOut,format="e",digits=2), stringsAsFactors=F)
	colnames(out) <- c( "PathName", "N_Genes", paste( "Fold", colnames(mOut), sep="_"), paste( "Pvalue", colnames(mOut), sep="_"))
	rownames(out) <- 1:nrow(out)

	# reorder to get the 2 values per group together
	outOrd <- c( 1:2, seq( 3, ncol(out), by=2), seq( 4, ncol(out), by=2))
	out <- out[ , outOrd]

	return(out)
}


`combineRadarPlotResults` <- function( radarPath) {

	# given a folder of radar plots and .CSV results, combine them into one final table
	fset <- dir( radarPath, pattern="^Radar\\..+\\.csv$", full.name=T, recursive=F)

	# don't include any previous combined results
	drops <- grep( "Radar.All.CombinedGeneSets", fset)
	if ( length(drops)) fset <- fset[ -drops]

	if ( ! length(fset)) return()

	out <- data.frame()
	for ( f in fset) {

		# grab each file, and make sure they are the same shape/content
		tbl <- read.csv( f, as.is=T)
		if ( !nrow(out)) {
			nCol <- ncol(tbl)
			colNames <- colnames(tbl)
		} else {
			if ( ncol(tbl) != nCol) next
			if ( any( colnames(tbl) != colNames)) next
		}

		# grap the gene set name from the file, to prepend to the joined names
		geneSetPrefix <- sub( "(^Radar\\.)(.+)(\\.csv$)", "\\2", basename(f))
		tbl$PathName <- paste( geneSetPrefix, tbl$PathName, sep=": ")
		tbl$PathName <- cleanGeneSetModuleNames( tbl$PathName, wrapParen=FALSE)

		out <- rbind( out, tbl)
	}

	# final ordering
	if ( "DeltaFold" %in% colnames(out)) {
		ord <- diffExpressRankOrder( out$DeltaFold, out$BestPvalue, wt.fold=1, wt.pvalue=1)
	} else {
		ord <- diffExpressRankOrder( out[[3]], out[[4]], wt.fold=1, wt.pvalue=1)
	}
	out <- out[ ord, ]
	if ( nrow(out)) rownames(out) <- 1:nrow(out)

	# Note about Pi Values:   their scaling is relative to each set's distribution of P values
	# to be most fair after combining, we should recalculate them
	if ( "BestPIvalue" %in% colnames(out)) {
		out$BestPIvalue <- piValue( out$DeltaFold, out$BestPvalue)
		myFoldCols <- grep( "^Fold_", colnames(out))
		for (k in myFoldCols) {
			# triplets of Fold, Pval, PIval
			out[[k+2]] <- piValue( out[[k]], out[[k+1]])
		}
	}

	outfile <- file.path( radarPath, "Radar.All.CombinedGeneSets.csv")
	write.table( out, outfile, sep=",", quote=T, row.names=F)

	return( out)
}

