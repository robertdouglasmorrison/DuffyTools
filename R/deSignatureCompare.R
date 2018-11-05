# deSignatureCompare.R -- compare the Diff Expression signatures of 2 sets of analyses
#				to build a way of calling how much does a set of 
#				DE results "Look Like" a specific other DE result


`deSignatureCompare` <- function( refPath, refGroup="", refLabel="Reference", 
				targetPath, targetGroup=NULL, targetLabel="Target", 
				speciesID=getCurrentSpecies(),
				min.fold=0.1, max.pvalue=0.1, max.members=2000, ...) {

	if ( refGroup == "") stop( "Explicit DE 'referenceGroup' is required.")
	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)

	# get the reference signature, and verify it was found
	cat( "\nGathering Reference Signatures..")
	refSig <- gatherAllSignatures( refPath, refGroup, min.fold=min.fold, max.pvalue=max.pvalue,
					max.members=max.members)
	if ( is.null( refSig$Gene)) stop( "Reference DE Gene results not found.")
	if ( length( refSig$Gene) > 1) stop( "Reference DE Gene results must be a single 'Group'.")

	# get the target signature, and verify something was found
	cat( "\n\nGathering Target Signatures..")
	targetSig <- gatherAllSignatures( targetPath, targetGroup, min.fold=min.fold, max.pvalue=max.pvalue,
					max.members=max.members)
	if ( is.null( targetSig$Gene)) stop( "No Target DE Gene results found.")
	if ( (N <- length(targetSig$Gene)) > 0) cat( "\n\nTarget Groups found: ", N, "\n", names(targetSig$Gene))

	# do the comparisons
	ans <- compareSignatures( refSig, targetSig, refLabel=refLabel, targetLabel=targetLabel, ...)
	cat( "\nDone.\n")

	return( ans)
}


`gatherAllSignatures` <- function( path, groupSet=NULL, min.fold=0.1, max.pvalue=0.1, max.members=2000) {

	geneFiles <- selectSignatureFiles( path, groupSet, mode="gene")
	densityFiles <- selectSignatureFiles( path, groupSet, mode="density")
	qusageFiles <- selectSignatureFiles( path, groupSet, mode="qusage")

	# grab the desired content from those files
	geneSignature <- NULL
	if ( nFiles <- length(geneFiles)) {
		geneSignature <- vector( mode="list")
		nSig <- 0
		for ( i in 1:nFiles) {
			smlDF <- gatherGeneSignature( geneFiles[i], min.fold=min.fold, max.pvalue=max.pvalue, 
						max.members=max.members)
			if ( ! nrow(smlDF)) next
			nSig <- nSig + 1
			geneSignature[[nSig]] <- smlDF
			names(geneSignature)[nSig] <- names(geneFiles)[i]
		}
	}

	# do the same for the GeneSet results
	densitySignature <- NULL
	if ( nFiles <- length(densityFiles)) {
		densitySignature <- vector( mode="list")
		nSig <- 0
		for ( i in 1:nFiles) {
			smlDF <- gatherDensitySignature( densityFiles[i], min.fold=min.fold, max.pvalue=max.pvalue, 
						max.members=max.members)
			if ( ! nrow(smlDF)) next
			nSig <- nSig + 1
			densitySignature[[nSig]] <- smlDF
			names(densitySignature)[nSig] <- names(densityFiles)[i]
		}
	}

	# do the same for the QuSage results
	qusageSignature <- NULL
	if ( nFiles <- length(qusageFiles)) {
		qusageSignature <- vector( mode="list")
		nSig <- 0
		for ( i in 1:nFiles) {
			smlDF <- gatherQusageSignature( qusageFiles[i], min.fold=min.fold, max.pvalue=max.pvalue*5, 
						max.members=max.members)
			if ( ! nrow(smlDF)) next
			nSig <- nSig + 1
			qusageSignature[[nSig]] <- smlDF
			names(qusageSignature)[nSig] <- names(qusageFiles)[i]
		}
	}

	# package it all up to send back
	out <- list( "Gene"=geneSignature, "Density"=densitySignature, "QuSage"=qusageSignature)
	return( out)
}


`selectSignatureFiles` <- function( path, groupSet=NULL, mode=c("gene","density","qusage")) {

	# create the filename pattern from the tool/mode
	mode <- match.arg( mode)
	filePattern <- "Meta.UP.txt$"
	if ( mode == "density") {
		path <- file.path( path, "CombinedGeneSets")
		filePattern <- "UP.CombinedGeneSets.txt$"
	} else if ( mode == "qusage") {
		path <- file.path( path, "QuSage")
		filePattern <- "UP.QuSage.CombinedGeneSets.csv$"
	}

	# add on the prefix to get just the right organism
	prefix <- getCurrentSpeciesFilePrefix()
	filePattern <- paste( prefix, filePattern, sep=".")

	# grab those file names
	files <- dir( path, pattern=filePattern, full.name=T)

	# give them names of their group name
	groupNames <- sub( filePattern, "", basename(files))
	groupNames <- sub( "\\.$", "", groupNames)
	names(files) <- groupNames

	# if given a group specifier, only send those back
	if ( ! is.null( groupSet)) {
		keep <- vector()
		for (grp in groupSet) {
			grpPattern <- paste( "^", grp, sep="")
			hit <- grep( grpPattern, basename(files))
			if ( length(hit)) keep <- c( keep, hit)
		}
		if ( length(keep)) files <- files[ sort( unique( keep))]
	}

	# done
	return( files)
}


`gatherGeneSignature` <- function( geneFile, min.fold=0.1, max.pvalue=0.1, max.members=2000) {

	# grab the genes that are DE UP, and 
	# make a small DF of the important values
	tbl <- read.delim( geneFile, sep="\t", as.is=T)
	tbl$RANK <- 1:nrow(tbl)
	tblUP <- tbl[ tbl$LOG2FOLD >= min.fold & tbl$AVG_PVALUE <= max.pvalue, ]
	N <- min( max.members, nrow(tblUP))
	if ( N < nrow(tblUP)) tblUP <- tblUP[ 1:N, ]

	# the META DOWN results are in a separate file
	geneFile2 <- sub( ".UP.", ".DOWN.", geneFile, fixed=T)
	tbl <- read.delim( geneFile2, sep="\t", as.is=T)
	tbl <- tbl[ rev( 1:nrow(tbl)), ]
	tbl$RANK <- 1:nrow(tbl)
	tblDOWN <- tbl[ tbl$LOG2FOLD <= -min.fold & tbl$AVG_PVALUE <= max.pvalue, ]
	N <- min( max.members, nrow(tblDOWN))
	if ( N < nrow(tblDOWN)) tblDOWN <- tblDOWN[ (nrow(tblDOWN)-N+1):nrow(tblDOWN), ]
	tbl <- rbind( tblUP, tblDOWN)
	N <- nrow(tbl)
	geneName <- shortGeneName( tbl$GENE_ID, keep=1)
	out <- data.frame( "GENE_ID"=geneName, "PRODUCT"=tbl$PRODUCT, "LOG2FOLD"=tbl$LOG2FOLD, 
				"RANK"=tbl$RANK, row.names=1:N, stringsAsFactors=F)
	cat( "\nGenes:    fold:",min.fold," pval:",max.pvalue," file:",basename(geneFile)," N=", nrow(out))
	return( out)
}


`gatherDensitySignature` <- function( densityFile, min.fold=0.1, max.pvalue=0.1, max.members=2000) {

	# grab the gene sets that are DE UP, and 
	# make a small DF of the important values
	tbl <- read.delim( densityFile, sep="\t", as.is=T)
	tbl$RANK <- 1:nrow(tbl)
	# make the column names
	groupName <- names(densityFile)[1]
	foldColumn <- grep( "Fold",colnames(tbl))[1]
	pvalColumn <- grep( "P.value",colnames(tbl))[1]
	tblUP <- tbl[ tbl[[foldColumn]] >= min.fold & tbl[[pvalColumn]] <= max.pvalue, ]
	N <- min( max.members, nrow(tblUP))
	if ( N < nrow(tblUP)) tblUP <- tblUP[ 1:N, ]

	# the META DOWN results are in a separate file
	densityFile2 <- sub( ".UP.", ".DOWN.", densityFile, fixed=T)
	tbl <- read.delim( densityFile2, sep="\t", as.is=T)
	tbl <- tbl[ rev( 1:nrow(tbl)), ]
	tbl$RANK <- 1:nrow(tbl)
	tblDOWN <- tbl[ tbl[[foldColumn]] <= -min.fold & tbl[[pvalColumn]] <= max.pvalue, ]
	N <- min( max.members, nrow(tblDOWN))
	if ( N < nrow(tblDOWN)) tblDOWN <- tblDOWN[ (nrow(tblDOWN)-N+1):nrow(tblDOWN), ]
	tbl <- rbind( tblUP, tblDOWN)
	N <- nrow(tbl)
	topGeneSets <- tbl$Name
	topGeneSets <- cleanGeneSetModuleNames( topGeneSets, wrapParentheses=FALSE)
	out <- data.frame( "GENESET_ID"=topGeneSets, "N_Genes"=tbl$N.Genes,
				"LOG2FOLD"=tbl[[ foldColumn]], "RANK"=tbl$RANK,
				row.names=1:N, stringsAsFactors=F)
	cat( "\nDensity:  fold:",min.fold," pval:",max.pvalue," file:",basename(densityFile)," N=", nrow(out))
	return(out)
}


`gatherQusageSignature` <- function( qusageFile, min.fold=0.1, max.pvalue=0.1, max.members=2000) {

	tbl <- read.csv( qusageFile, sep=",", as.is=T)
	tbl$RANK <- 1:nrow(tbl)
	tblUP <- subset( tbl, LOG2FOLD >= min.fold & PVALUE <= max.pvalue)
	N <- min( max.members, nrow(tblUP))
	if ( N < nrow(tblUP)) tblUP <- tblUP[ 1:N, ]

	tblDOWN <- subset( tbl, LOG2FOLD <= -min.fold & PVALUE <= max.pvalue)
	N <- min( max.members, nrow(tblDOWN))
	if ( N < nrow(tblDOWN)) tblDOWN <- tblDOWN[ (nrow(tblDOWN)-N+1):nrow(tblDOWN), ]
	tbl <- rbind( tblUP, tblDOWN)
	N <- nrow(tbl)
	topQusageSets <- tbl$PATHWAY
	topQusageSets <- cleanGeneSetModuleNames( topQusageSets, wrapParentheses=FALSE)
	out <- data.frame( "GENESET_ID"=topQusageSets, "N_Genes"=tbl$N_GENES,
				"LOG2FOLD"=tbl$LOG2FOLD, "RANK"=tbl$RANK,
					row.names=1:N, stringsAsFactors=F)
	cat( "\nQusage:   fold:",min.fold," pval:",max.pvalue," file:",basename(qusageFile)," N=", nrow(out))
	return(out)
}


`compareSignatures` <- function( refSig, targetSig, refLabel="", targetLabel="", targetColor=NULL, 
				targetOrder=NULL, legend.cex=1, ...) {

	comparePlotSetup( refLabel, targetLabel, ...)

	# grab the parts of the signatures
	refGene <- refSig$Gene[[1]]
	refDensity <- refSig$Density[[1]]
	refQusage <- refSig$QuSage[[1]]
	targetGene <- targetSig$Gene
	targetDensity <- targetSig$Density
	targetQusage <- targetSig$QuSage
	groupNames <- names( targetGene)
	nGroups <- length( groupNames)

	# allow control of colors and legend ordering
	#if ( is.null( targetOrder)) targetOrder <- 1:nGroups
	if ( is.null( targetColor)) {
		myOrder <- if ( is.null( targetOrder)) 1:nGroups else targetOrder
		targetColor <- vector( length=nGroups)
		targetColor[ myOrder] <- rainbow( nGroups, end=0.7)
	}

	# we will compare our one reference against all the target groups
	outCor <- outPval <- outN <- outNgene <- outNdensity <- outNqusage <- vector( length=nGroups)
	outDetails <- vector( mode="list")

	for ( i in 1:nGroups) {
		thisGrp <- groupNames[i]
		thisGeneSig <- targetGene[[i]]
		f1 <- f2 <- vector()
		nGene <- nDensity <- nQusage <- 0
		thisDetails <- data.frame()

		# we always have genes
		ans <- compareOneSignature( refGene, thisGeneSig, groupName=thisGrp, 
					idColumn="GENE_ID", col=targetColor[i], pch=21, ...)
		if ( ! is.null(ans)) {
			f1 <- c( f1, ans$Fold.Ref)
			f2 <- c( f2, ans$Fold.Target)
			nGene <- length( ans$Fold.Ref)
			cat( "\r  Add Gene:   ", thisGrp, nGene)
			# the gene name terms in the answer are not as useful as the Density/QuSage terms
			# so fill them out
			ans$Name <- paste( ans$Name, gene2Product( ans$Name), sep=": ")
			thisDetails <- rbind( thisDetails, data.frame( "Class"="Gene", ans, stringsAsFactors=F))
		}
		
		# we may have gene set density
		if ( length( targetDensity)) {
			who <- match( thisGrp, names( targetDensity), nomatch=0)
			if (who) {
				thisSig <- targetDensity[[who]]
				ans <- compareOneSignature( refDensity, thisSig, groupName=thisGrp, 
							idColumn="GENESET_ID", col=targetColor[i], pch=23, ...)
				if ( ! is.null(ans)) {
					f1 <- c( f1, ans$Fold.Ref)
					f2 <- c( f2, ans$Fold.Target)
					nDensity <- length( ans$Fold.Ref)
					cat( "\r  Add Density: ", thisGrp, nDensity)
					thisDetails <- rbind( thisDetails, data.frame( "Class"="Density", ans, stringsAsFactors=F))
				}
			}
		}
		# we may have gene set QuSage
		if ( length( targetQusage)) {
			who <- match( thisGrp, names( targetQusage), nomatch=0)
			if (who) {
				thisSig <- targetQusage[[who]]
				ans <- compareOneSignature( refQusage, thisSig, groupName=thisGrp, 
							idColumn="GENESET_ID", col=targetColor[i], pch=24, ...)
				if ( ! is.null(ans)) {
					f1 <- c( f1, ans$Fold.Ref)
					f2 <- c( f2, ans$Fold.Target)
					nQusage <- length( ans$Fold.Ref)
					cat( "\r  Add QuSage:  ", thisGrp, nQusage)
					thisDetails <- rbind( thisDetails, data.frame( "Class"="QuSage", ans, stringsAsFactors=F))
				}
			}
		}

		# OK, we have all the data we need to make a correlation call
		cc <- cor.test( f1, f2)
		outCor[i] <- cc$estimate
		outPval[i] <- cc$p.value
		outNgene[i] <- nGene
		outNdensity[i] <- nDensity
		outNqusage[i] <- nQusage
		outN[i] <- nGene + nDensity + nQusage

		# prep the details
		ord <- order( thisDetails$Delta.Fold, thisDetails$Class, thisDetails$Name)
		thisDetails <- thisDetails[ ord, ]
		rownames(thisDetails) <- 1:nrow(thisDetails)
		outDetails[[ i]] <- thisDetails
	}

	# after every group is done, redraw more randomly
	comparePlotFinish( targetOrder=targetOrder, legend.cex=legend.cex, ...)

	# package up the final results
	out <- data.frame( "Group"=groupNames, "Pearson.R"=round(outCor,digits=3), 
				"Pvalue"=as.numeric( formatC( outPval, format="e", digits=3)), 
				"N_Points"=outN, "N_Genes"=outNgene, "N_Density"=outNdensity,
				"N_QuSage"=outNqusage, stringsAsFactors=F)
	ord <- diffExpressRankOrder( out$Pearson.R, out$Pvalue)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	# also make the details visible
	names(outDetails) <- groupNames
	DE_SIG_DETAILS <<- outDetails
	return( out)
}


`comparePlotSetup` <- function( refLabel, targetLabel, clip.fold=3.0, ...) {

	checkX11()
	mainText <- paste( "DE Signature Compare:    ", refLabel, "  .vs.  ", targetLabel)

	plot( 1, 1, type="n", main=mainText, xlim=c(-clip.fold,clip.fold), ylim=c(-clip.fold,clip.fold),
		xlab=paste( "Fold Change in", refLabel), ylab=paste( "Fold Change in Groups of", targetLabel))

	# set up a proxy for perfect agreement
	lines( clip.fold*c(-2,2), clip.fold*c(-2,2), col=1, lty=2, lwd=3)
	text( clip.fold, clip.fold, paste( "Perfect Correlation:\n", refLabel, " "), col=1, pos=2, cex=0.8)

	# set some storage for replotting after
	COR_SIG_X <<- COR_SIG_Y <<- COR_SIG_COL <<- COR_SIG_PCH <<- COR_SIG_GRP <<- vector()
	return()
}


compareOneSignature <- function( refDF, targetDF, groupName, mode=c("intersect","union"), idColumn="GENE_ID", 
				col=1, pch=1, pt.cex=1, ...) {

	# grab the identifiers from both
	refID <- refDF[[ idColumn]]
	targetID <- targetDF[[ idColumn]]

	# get the set of IDs in common, missing will be set to zero, zero
	# and grap the observed fold change data
	mode <- match.arg( mode)
	if (mode == "union") {
		allID <- sort( union( refID, targetID))
	} else {
		allID <- sort( intersect( refID, targetID))
	}
	NID <- length(allID)

	# find and fill the 2 vectors of fold change data
	refFold <- targetFold <- rep.int( 0, NID)
	whRef <- match( allID, refID, nomatch=0)
	refFold[ whRef > 0] <- refDF$LOG2FOLD[ whRef]
	whTarget <- match( allID, targetID, nomatch=0)
	targetFold[ whTarget > 0] <- targetDF$LOG2FOLD[ whTarget]

	points( refFold, targetFold, col=col, pch=pch, cex=pt.cex)
	dev.flush()

	COR_SIG_X <<- c( COR_SIG_X, refFold)
	COR_SIG_Y <<- c( COR_SIG_Y, targetFold)
	COR_SIG_COL <<- c( COR_SIG_COL, rep.int( col, NID))
	COR_SIG_PCH <<- c( COR_SIG_PCH, rep.int( pch, NID))
	COR_SIG_GRP <<- c( COR_SIG_GRP, rep.int( groupName, NID))

	out <- data.frame( "Name"=allID, "Fold.Ref"=refFold, "Fold.Target"=targetFold, 
			"Delta.Fold"=abs(refFold-targetFold), stringsAsFactors=F)
	return( out)
}


`comparePlotFinish` <- function( pt.cex=1, legend.cex=1, targetOrder=NULL, ...) {

	# redraw all the groups in random order
	ord <- sample( length(COR_SIG_X))
	points( COR_SIG_X[ord], COR_SIG_Y[ord], col=COR_SIG_COL[ord], pch=COR_SIG_PCH[ord], cex=pt.cex)

	# and redraw the final correlation lines
	myGroups <- sort( unique( COR_SIG_GRP))
	myColors <- COR_SIG_COL[ match( myGroups, COR_SIG_GRP)]
	nGroups <- length( myGroups)
	myCors <- vector()
	for ( i in 1:nGroups) {
		myGrp <- myGroups[i]
		myCol <- myColors[i]
		who <- which( COR_SIG_GRP == myGrp)
		myX <- COR_SIG_X[who]
		myY <- COR_SIG_Y[who]
		myCors[i] <- cor( myX, myY)
		abline( lsfit( myX, myY), col=myCol, lwd=3, lty=1)
	}

	if ( is.null( targetOrder)) targetOrder <- order( myCors, decreasing=T)
	ord <- targetOrder
	legendText <- paste( padStringWidth(myGroups[ord]), "  R=", formatC( myCors[ord], format="f", digits=3, flag="+"))
	legend( "topleft", legendText, lwd=3, lty=1, col=myColors[ord], bg='white', cex=legend.cex)
	dev.flush()
}

