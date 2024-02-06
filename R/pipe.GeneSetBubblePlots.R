# pipe.GeneSetBubblePlots.R -- visualize gene set DE as a family of bubble plots

`pipe.GeneSetBubblePlots` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, folderName=NULL,
				tool=c( "MetaResults", "DESeq", "EdgeR", "RankProduct", "RoundRobin", "SAM"),
				groupColumn="Group", groupLevels=NULL, average.FUN=median,
				geneSets=defaultGeneSets(speciesID), restrictionSets=NULL, baselineGroup=NULL,
				label.cex=1, Nshow=24, cex=1, main=paste( "Comparison:  ",folderName),
				max.label.length=80, xBubbleLim=c(0.6,0.9))
{

	checkX11( bg='white', width=9, height=7)
	require( gplots)
	require( RColorBrewer)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# get the SampleIDs that we know are valid
	annT <- readAnnotationTable( annotationFile)
	allSamples <- sampleIDset
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allSamples <- myAnnT$SampleID
	if ( ! (groupColumn %in% colnames(myAnnT))) stop( paste( "Given grouping column not in annotation table: ", groupColumn))
	allGroups <- myAnnT[[ groupColumn]]

	# build the path to the DE results we will use, and where we will put our Bubble plot results
	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", verbose=F)
	}
	if ( is.null( folderName)) stop( "Explicit Folder Name of DE results is required")
	folderName <- paste( prefix, folderName, sep=".")
	tool <- match.arg( tool)
	dePath <- file.path( results.path, tool, folderName)
	if ( ! file.exists( dePath)) stop( paste( "DE results folder not found: ", dePath))
	bubblePath <- file.path( dePath, "BubblePlots")
	if ( ! file.exists( bubblePath)) dir.create( bubblePath, recursive=T)

	# we have to read the transcriptomes in... to pre-load the data
	needLoad <- TRUE
	if (needLoad) {
		# still to do... get 'multi' vs 'unique' call?...
		intensityColumn <- getExpressionUnitsColumn( optionsFile, useMultiHits=TRUE)
		cat( "\nLoading ", length( allSamples), "transcriptomes..")
		files <- paste( allSamples, prefix, "Transcript.txt", sep=".")
		files <- file.path( results.path, "transcript", files)
		bubbleM <<- expressionFileSetToMatrix( files, allSamples, intensityColumn=intensityColumn, verbose=T)
		cat( "\nConverting Expression Abundance to M-values..")
		# if we were given a baseline group, assert that at the Mvalue step, not at the ReductionToModules step
		# also, let's try using the grouping info at the Mvalue step, to get an average that is 'between the groups'
		baselineColumns <- NULL
		if ( ! is.null(baselineGroup)) baselineColumns <- which( allGroups == baselineGroup)
		bubbleMA <<- expressionMatrixToMvalue( bubbleM, average.FUN=average.FUN, groupNames=allGroups, 
						baselineColumns=baselineColumns)
		cat( "  Done.\n")
	}


	# now we are read to make those bubble plots.   We were either given one gene set (list) or
	# a vector of character gene set names
	if ( is.list( geneSets)) {
		ans <- pipe.OneBubblePlot( sampleIDset=allSamples, speciesID=speciesID, annotationFile=annotationFile,
				optionsFile=optionsFile, results.path=results.path, groupColumn=groupColumn, groupLevels=groupLevels, 
				geneSet=geneSets, restrictionSets=restrictionSets,
				baselineGroup=baselineGroup, reload=FALSE, Nshow=Nshow, 
				label.cex=label.cex, cex=cex, main=main, max.label.length=max.label.length,
				xBubbleLim=xBubbleLim)
		geneSetName <- "Bubble"
		plotFile <- file.path( bubblePath, geneSetName)
		printPlot( plotFile, optT=optionsFile)
		csvFile <- file.path( bubblePath, paste( geneSetName, "csv", sep="."))
		write.table( ans, csvFile, sep=",", quote=T, row.names=F)
		nOut <- min( Nshow, nrow(ans))
		# if the full set is just a bit larger, go ahead and keep them all?...
		if ( Nshow * 1.25 > nrow(ans)) nOut <- nrow(ans)
		out <- ans[ 1:nOut, ]
	} else {
		# given a vector of names of gene sets, do each one
		out <- vector( mode="list")
		for (k in 1:length(geneSets)) {
			gs <- geneSets[k]
			cat( "\nDoing Bubble Plot for:  ", gs)
			thisRestrictionSet <- restrictionSets
			if ( ! is.null( restrictionSets)) {
				where <- match( gs, names(restrictionSets), nomatch=0)
				if ( where > 0) {
					thisRestrictionSet <- restrictionSets[[ where]]
					cat( "\n  Using a restriction set of pathways..")
				}
			}

			ans <- pipe.OneBubblePlot( sampleIDset=allSamples, speciesID=speciesID, annotationFile=annotationFile,
					optionsFile=optionsFile, results.path=results.path, groupColumn=groupColumn, groupLevels=groupLevels, 
					geneSet=gs, restrictionSets=thisRestrictionSet,
					baselineGroup=baselineGroup, reload=FALSE, Nshow=Nshow, 
					label.cex=label.cex, cex=cex, main=main, max.label.length=max.label.length, 
					xBubbleLim=xBubbleLim)
			if ( ! nrow( ans)) {
				out[[k]] <- ans
				next
			}
			plotFile <- file.path( bubblePath, paste( "Bubble", gs, sep="."))
			printPlot( plotFile, optT=optionsFile)
			csvFile <- file.path( bubblePath, paste( "Bubble", gs, "csv", sep="."))
			write.table( ans, csvFile, sep=",", quote=T, row.names=F)
			nOut <- min( Nshow, nrow(ans))
			# if the full set is just a bit larger, go ahead and keep them all?...
			if ( Nshow * 1.25 > nrow(ans)) nOut <- nrow(ans)
			smallOut <- ans[ 1:nOut, ]
			out[[k]] <- smallOut
		}
		names( out) <- geneSets
		cat( "\nMade Bubble Plots for", length(geneSets), "GeneSets.\n")
	}

	return( invisible( out))
}


`pipe.OneBubblePlot` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL,  
				groupColumn="Group", groupLevels=NULL, average.FUN=median, 
				geneSet=defaultGeneSets(speciesID), restrictionSets=NULL, baselineGroup=NULL,
				reload=FALSE, Nshow=24, label.cex=1, cex=par("cex.axis"), main=NULL, 
				max.label.length=80, crop.min.pvalue=1e-30, las=3, xBubbleLim=c(0.6,0.9))
{

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	mainText <- "Bubble Plot:  "
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

	# use the samples and grouping column to know the files to load
	annT <- readAnnotationTable( annotationFile)
	allSamples <- sampleIDset
	myAnnT <- subset( annT, SampleID %in% allSamples)
	allSamples <- myAnnT$SampleID
	if ( ! (groupColumn %in% colnames(myAnnT))) stop( paste( "Given grouping column not in annotation table: ", groupColumn))
	allGroups <- myAnnT[[ groupColumn]]
	if ( is.null( groupLevels)) groupLevels <- sort( unique( allGroups))
	grpFac <- factor( allGroups, levels=groupLevels)
	nGroups <- nlevels( grpFac)
	groupNames <- levels( grpFac)
	allGroupsList <- tapply( myAnnT$SampleID, grpFac, FUN=c, simplify=FALSE)
	NperGroup <- sapply( allGroupsList, length)

	if ( is.null( results.path)) {
		results.path <- getOptionValue( optionsFile, "results.path", verbose=F)
	}

	# see if we need to read the transcriptomes in...
	needLoad <- TRUE
	if ( exists( "bubbleM") && all( colnames(bubbleM) == allSamples) && nrow(bubbleM) > 100 && !reload) needLoad <- FALSE
	if (needLoad) {
		intensityColumn <- getExpressionUnitsColumn( optionsFile, useMultiHits=TRUE)
		files <- paste( allSamples, prefix, "Transcript.txt", sep=".")
		files <- file.path( results.path, "transcript", files)
		bubbleM <<- expressionFileSetToMatrix( files, allSamples, intensityColumn=intensityColumn, verbose=T)
		# if we were given a baseline group, assert that at the Mvalue step, not at the ReductionToModules step
		baselineColumns <- NULL
		if ( ! is.null(baselineGroup)) baselineColumns <- which( allGroups == baselineGroup)
		bubbleMA <<- expressionMatrixToMvalue( bubbleM, average.FUN=average.FUN, baselineColumns=baselineColumns)
	}

	# try to standarize the names of the gene sets to keep them easy to view
	fullPathNames <- names(allGeneSets)
	names(allGeneSets) <- cleanGeneSetModuleNames( names(allGeneSets), wrapParentheses=F)

	# we may have been given a restriction set, to only consider gene sets from some previous run
	if ( ! is.null( restrictionSets)) {
		pathNames <- restrictionSets$PathName
		if ( ! is.null( pathNames)) {
			pathNames <- cleanGeneSetModuleNames( pathNames, wrapParentheses=F)
			keepers <- which( names( allGeneSets) %in% pathNames)
			if ( length( keepers)) {
				allGeneSets <- allGeneSets[ keepers]
				fullPathNames <- fullPathNames[ keepers]
			}
		}
	}

	# do the reduction & grouping.  This call must use mean average, never pass down the method used for the MA step.
	bubbleAns <- reduceMatrixToModules( bubbleMA, geneModules=allGeneSets, sampleTraits=allGroupsList,
				gene.names=shortGeneName( rownames( bubbleMA), keep=1), average.FUN=mean,
				sample.names=colnames(bubbleMA))
	mShow <- bubbleMOD <- bubbleAns$matrix
	pShow <- bubblePvalue <- bubbleAns$p.value
	validModuleNames <- bubbleAns$moduleNames
	nGenes <- bubbleAns$geneCounts

	# the reduction may have thrown out some modules
	keep <- which( names(allGeneSets) %in% validModuleNames)
	fullPathNames <- fullPathNames[keep]

	# trim down to less groups if we need to
	Nshow <- min( Nshow, nrow( bubbleMOD))
	# if the full set is just a bit larger, go ahead and keep them all?...
	if ( Nshow * 1.25 > nrow(bubbleMOD)) Nshow <- nrow(bubbleMOD)

	# use both Pvalue and magnitudes to decide which subset to draw.
	# if all the Module signs point in th same direction (especially when we used a baseline group),
	# down weight those elements
	magnitudes <- apply( bubbleMOD, 1, function(x) {
			dx <- diff( range( x, na.rm=T))
			xSign <- sign(x)
			if ( all(sign(xSign) == 1,na.rm=T)) dx <- dx/2
			if ( all(sign(xSign) == -1,na.rm=T)) dx <- dx/2
			return(dx)
		})
	bestPs <- apply( bubblePvalue, 1, min)
	# but give magnitudes more weight?...
	ord <- diffExpressRankOrder( magnitudes, bestPs, wt.fold=5, wt.pvalue=1)
	mShow <- bubbleMOD[ ord[1:Nshow], ]
	pShow <- bubblePvalue[ ord[1:Nshow], ]
	mOut <- bubbleMOD[ ord, ]
	pOut <- bubblePvalue[ ord, ]
	magniOut <- magnitudes[ ord]
	bestPout <- bestPs[ ord]
	validModuleNames <- validModuleNames[ ord[ 1:Nshow]]
	fullPathNames <- fullPathNames[ ord]
	nGenes <- nGenes[ ord]

	# now alphabetize on any part after the Module names.  We will draw from bottom up, so reverse them.
	ord <- rev( order( sub( "M.+: +", "", validModuleNames), sub( "(^M[0-9]+\\.[0-9]+:)(.+)", "\\1", validModuleNames)))
	mShow <- mShow[ ord, ]
	pShow <- pShow[ ord, ]
	validModuleNames <- validModuleNames[ ord]
	
	# when given a baseline group that all other groups were compared against, remove
	# that group from what we show
	if ( ! is.null(baselineGroup)) {
		dropCol <- which( groupNames == baselineGroup)[1]
		mShow <- mShow[ , -dropCol]
		pShow <- pShow[ , -dropCol]
		groupNames <- groupNames[ -dropCol]
		nGroups <- nGroups - 1
	}

	# plot it now, as a box of colored dots, with labels on the left
	# use the passed in x limits for where the bubbles go
	xBubbleBoxLeft <- xBubbleLim[1]
	xBubbleBoxRight <- xBubbleLim[2]
	xBubbleBoxWidth <- xBubbleBoxRight - xBubbleBoxLeft
	xBubbleHalfGap <- xBubbleBoxWidth / (nGroups*1.75)
	xGroupLocs <- seq( xBubbleBoxLeft+xBubbleHalfGap, xBubbleBoxRight-xBubbleHalfGap, length.out=nGroups) 
	plot( 1, 1, type="n", main=mainText, xlab=NA, ylab=NA, xaxt="n", yaxt="n", xlim=c(0,1),
			ylim=c(0,Nshow+2), yaxs="i", frame.plot=F)
	rect( xBubbleBoxLeft, 0, xBubbleBoxRight, Nshow+1, border=1, col='grey99')
	axis( side=1, at=xGroupLocs, label=groupNames, line=F, tick=T, las=las)
	for ( xx in 1:nGroups) lines( rep.int(xGroupLocs[xx],2), c(0,Nshow+1), col='grey70', lwd=0.5, lty=1)
	for ( yy in 1:Nshow) lines( c(xBubbleBoxLeft,xBubbleBoxRight), c(yy,yy), col='grey70', lwd=0.5, lty=1)
	text( rep.int(xBubbleBoxLeft,Nshow), 1:Nshow, clipLongString(validModuleNames,max.length=max.label.length,pct.front=0.8), 
			cex=label.cex, pos=2, offset=1.0, col=1)
	
	# set up to draw the bubbles, using a fixed palette for now
	nColors <- 51
	palette <- rev( brewer.pal( 9, "RdBu"))
	getPalette <- colorRampPalette(palette)
	# get the fold range, and force the fold range to be symmetric about zero
	maxPosMagni <- max( mShow, 1, na.rm=T) * 1.01
	maxNegMagni <- min( mShow, -1, na.rm=T) * 1.01
	if ( abs(maxNegMagni) > maxPosMagni) maxPosMagni <- abs(maxNegMagni)
	if ( maxNegMagni > (-maxPosMagni)) maxNegMagni <- (-maxPosMagni)
	cropPshow <- pmax( pShow, crop.min.pvalue)
	minPval <- min( cropPshow, 0.0001, na.rm=T)
	# visit each point on the bubble plot rectange, and draw that significance.
	for ( j in 1:nGroups) {
		# inspect the size of the bubbles, to draw them in an order that best shows overlaps
		ord <- order( cropPshow[ , j])
		for ( i in ord) {
			myMagni <- mShow[ i, j]
			myPval <- cropPshow[ i, j]
			if ( myPval > 0.95) next
			if (myMagni > 0) {
				myColorPtr <- nColors - round( (1-(myMagni/maxPosMagni))*nColors/2)
			} else {
				myColorPtr <- round( (1-(myMagni/maxNegMagni))*nColors/2) + 1
			}
			myCol <- getPalette(nColors)[myColorPtr]
			myCEX <-  (-log10(myPval)) / (-log10(minPval)) * 5
			points( xGroupLocs[j], i, pch=21, cex=myCEX*cex, col='grey40', lwd=0.5, bg=myCol)
		}
	}

	# make the legends
	legLeft <- xBubbleBoxRight + (1-xBubbleBoxRight)*0.2
	legRight <- xBubbleBoxRight + (1-xBubbleBoxRight)*0.6
	foldCorners <- c( legLeft, Nshow*0.6, legRight, Nshow*0.8)
	foldColors <- getPalette(nColors)
	yStep <- diff(foldCorners[c(2,4)]) / nColors
	for (k in 1:nColors) rect( foldCorners[1], foldCorners[2]+(yStep*(k-1)), foldCorners[3], foldCorners[2]+(yStep*k), border=foldColors[k], col=foldColors[k])
	rect( foldCorners[1], foldCorners[2], foldCorners[3], foldCorners[4], border=1, lwd=1, col=NA)
	text( foldCorners[3], foldCorners[4], "Log2 Fold", pos=3, offset=0.65, cex=0.85)
	prettyVals <- pretty( c( maxNegMagni, 0, maxPosMagni), 4)
	magSteps <- ((prettyVals - maxNegMagni) / (maxPosMagni - maxNegMagni))
	ySteps <- diff(range(foldCorners[c(2,4)])) * magSteps
	yTicks <- foldCorners[2] + ySteps
	show <- which( prettyVals >= maxNegMagni & prettyVals <= maxPosMagni)
	#for (yy in yTicks) lines( foldCorners[3]+c(0,0.1), rep.int(yy,2), col=1, lwd=1)
	text( rep.int(foldCorners[3],length(show)), yTicks[show], prettyVals[show], pos=4, offset=0.35, cex=0.85)

	pvalCorners <- c( legLeft, Nshow*0.18, legRight, Nshow*0.45)
	text( pvalCorners[3], pvalCorners[4], "-Log10(P) ", pos=3, offset=0.25, cex=0.85)
	prettyVals <- unique( c( 1, pretty( -log10( c( 0.05, minPval)), 7)))
	use <- which( prettyVals >= 1 & prettyVals <= -log10(minPval))
	# have the step size grow as the circles grow
	pvalStep <- diff( range( pvalCorners[c(2,4)])) / (length(use)*3.5)
	ynow <- pvalCorners[4]
	for( kk in 1:length(use)) {
		jj <- use[kk]
		myCEX <-  ( prettyVals[jj] / (-log10(minPval))) * 5
		ynow <- ynow - (pvalStep * kk)
		points( mean(pvalCorners[c(1,3)]), ynow, cex=myCEX*cex, pch=1, col=1, lwd=1)
		text( pvalCorners[3], ynow, prettyVals[jj], col=1, pos=4, offset=0.35, cex=0.85)
	}
	
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
				mOut, formatC( pOut, format="e", digits=2),  stringsAsFactors=F)
	colnames(out) <- c( "PathName", "N_Genes", "DeltaFold", "BestPvalue", 
				paste( "Fold", colnames(mOut), sep="_"), paste( "Pvalue", colnames(mOut), sep="_"))
	rownames(out) <- 1:nrow(out)

	# reorder to get the 2 values per group together
	outOrd <- c( 1:4)
	for ( k in 1:ncol(mOut)) outOrd <- c( outOrd, (k+4), (ncol(mOut)+k+4))
	out <- out[ , outOrd]

	return(out)
}
