# pipe.qusage -- wrapper to the 'QuSage' gene set enrichment package

`pipe.QuSage` <- function( sampleIDset, speciesID=getCurrentSpecies(), annotationFile="Annotation.txt",
				optionsFile="Options.txt", results.path=NULL, folderName=NULL,
				tool=c( "MetaResults", "DESeq", "EdgeR", "RankProduct", "RoundRobin", "SAM"),
				m=NULL, contrast=NULL, useLog=TRUE,
				groupColumn="Group", geneSets=defaultGeneSets(speciesID), descriptor="CombinedGeneSets", 
				minGenesPerSet=if (speciesID %in% MAMMAL_SPECIES) 4 else 2, 
				min.rpkm=1, min.variance=1, offset=2,
				addCellTypes=(speciesID %in% MAMMAL_SPECIES), addLifeCycle=(speciesID %in% PARASITE_SPECIES), 
				doQuSage= TRUE, verbose=TRUE, n.points=2^12, ...)
{

	require( qusage)

	if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	# we need a subfolder for where to stuff these results
	if ( is.null( folderName)) stop( "Required argument 'folderName' not specified")

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\nStarting pipe 'QuSage' on Sample Set: \n")
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

	# small chance we onle want to revisit existing results...
	if (doQuSage) {
		# get the matrix of gene expression
		if ( is.null(m)) {
			if ( speciesID != getCurrentSpecies()) setCurrentSpecies( speciesID)
			fset <- file.path( results.path, "transcript", paste( sids, prefix, "Transcript.txt", sep="."))
			if (verbose) cat( "\nGathering Gene Expression matrix..")
			m <- expressionFileSetToMatrix( fset, sids)
		} else {
			# verify the matrix matches the samples given
			if ( ! all( make.names(sids) == colnames(m))) stop( "Column names of 'm' are not 'sampleIDset'")
		}
		rownames(m) <- shortGeneName( rownames(m), keep=1)
		mOrig <- m
	
		# QuSage dies on NA
		isNA <- apply( m, MARGIN=1, function(x) any( is.na(x)))
		drops <- which( isNA)
		if ( length( drops)) {
			m <- m[ -drops, ]
			cat( "\nDropping rows with 'NA' values.   N=", length(drops))
		}

		# QuSage suggest removal of low expression, and wants log2 data
		bigV <- apply( m, MARGIN=1, max, na.rm=T)
		drops <- which( bigV < min.rpkm)
		if ( length( drops)) {
			m <- m[ -drops, ]
			cat( "\nDropping low expression rows below: ", min.rpkm, " RPKM:   N=", length(drops))
			bigV <- bigV[ -drops]
		}
		# QuSage also dies if the variance is too small
		smallV <- apply( m, MARGIN=1, min, na.rm=T)
		rangeV <- bigV - smallV
		drops <- which( rangeV < min.variance)
		if ( length( drops)) {
			m <- m[ -drops, ]
			cat( "\nDropping low variance rows below: ", min.variance, " RPKM:   N=", length(drops))
		}

		# QuSage chokes on duplicate IDs
		drops <- which( duplicated( rownames(m)))
		if ( length(drops)) m <- m[ -drops, ]

		# debugging test,  QuSage crashes on some datasets..  
		# Authors think perhaps flat 'all zeros / no variance' data might cause that.  
		# Try to remove that possibility
		#isZero <- ( m < 0.5)
		#if ( nZero <- sum(isZero)) {
			#jitterOffset <- runif( nZero, 0.1, 0.25)
			#m[isZero] <- m[isZero] + jitterOffset
			#cat( "\nJittering zero expression to prevent QuSage crash: ", nZero, "zero values")
		#}
	
		# do the log transform with offset
		if (useLog) {
			cat( "\nApplying log2 transform after adding linear offset of", offset, "RPKM")
			m <- log2( m + offset)
		}
	}

	# allow several way of giving gene sets...
	GSanswer <- gatherGeneSets( geneSets, descriptor, mode="combined")
	bigSetOfGeneSets <- GSanswer$geneSets
	myGeneSets <- bigSetOfGeneSets[[1]]
	if ( is.character( geneSets) && length(geneSets) == 1) descriptor=geneSets[1]

	# force the check/drop of gene sets with too few genes in them, for this matrix of genes
	gsDrops <- vector()
	for ( i in 1:length(myGeneSets)) {
		nNow <- length( intersect( myGeneSets[[i]], rownames(m)))
		if ( nNow < minGenesPerSet) gsDrops <- c( gsDrops, i)
	}
	if ( length( gsDrops)) myGeneSets <- myGeneSets[ -gsDrops]
	myGeneSetNames <- names(myGeneSets)
	myGeneSetShortNames <- names(myGeneSets) <- cleanGeneSetModuleNames( myGeneSetNames, wrapParentheses=FALSE)

	# the contrast string
	if ( ! is.null( contrast)) {
		EXPLICIT_CONTRAST <- TRUE
		terms <- strsplit( contrast, split="-", fixed=T)[[1]]
		if ( ! all( terms %in% allGroups)) {
			cat( "\nError:  Explicit Contrast string does not match given group names.")
			cat( "\nContrast String:  ", contrast)
			cat( "\nGroups defined:   ", allGroups)
			stop()
		}
		# notation is "case - control", so let's talk about UP/DOWN in 'Case'
		group1 <- terms[1]
		group2 <- terms[2]
	} else {
		EXPLICIT_CONTRAST <- FALSE
	}

	# ready to call QuSage
	cat( "\n")
	tool <- match.arg( tool)
	path <- file.path( results.path, tool, paste( prefix,folderName,sep="."), "QuSage")
	if ( ! file.exists( path)) dir.create( path, recursive=T)

	if (EXPLICIT_CONTRAST) {
		cat( "\nContrast Grouping descriptor:  ", contrast, "\n")
		ans <- do.QuSage( eset=m, labels=grps, contrast=contrast, geneSets=myGeneSets, descriptor=descriptor,
				group1=group1, group2=group2, path=path, prefix=prefix, makeDownHTML=TRUE,
				shortNames=myGeneSetShortNames, longNames=myGeneSetNames, addCellTypes=addCellTypes,
				addLifeCycle=addLifeCycle, doQuSage=doQuSage, n.points=n.points, ...)
		out <- list( ans, "Matrix"=mOrig)
		return(out)

	} else {
		# do all pairs of groups, allow it to be multi-threaded
		multicore.qusage <- function( g) {
				group1 <- g
				group2 <- paste( "Not", g, sep="")
				myGroups <- grps
				myGroups[ myGroups != group1] <- group2
				myContrast <- paste( group1, group2, sep="-")
				cat( "\nContrast Grouping descriptor:  ", myContrast, "\n")
				do.QuSage( eset=m, labels=myGroups, contrast=myContrast, geneSets=myGeneSets, descriptor=descriptor,
					group1=group1, group2=group2, path=path, prefix=prefix, makeDownHTML=(Ngroups>2),
					shortNames=myGeneSetShortNames, longNames=myGeneSetNames, addCellTypes=addCellTypes,
					addLifeCycle=addLifeCycle, doQuSage=doQuSage, n.points=n.points, ...)
				return(NULL)
			}
		multicore.lapply( allGroups, multicore.qusage)
	}

	if (verbose) {
		cat( "\n\n=============================")
		cat( "\n\nFinished pipe 'QuSage' on Sample Set:     ", unlist(sampleIDset), 
			"\nSpecies: ", speciesID,"\nGroups: ", allGroups, "\n")
	}
	return()
}


do.QuSage <- function( eset, labels, contrast, geneSets, descriptor="GeneSets", group1="Group1", group2="Group2", 
			path=".", prefix=getCurrentSpeciesFilePrefix(), makeDownHTML=TRUE, 
			shortNames=names(geneSets), longNames=names(geneSets), addCellTypes=FALSE, addLifeCycle=FALSE,
			doQuSage=TRUE, n.points=2^12, ...) {

	# having less than 2 items in any group will break qusage.  Check
	grpTable <- table( labels)
	if ( any( grpTable < 2)) {
		bad <- names(grpTable)[ grpTable < 2]
		cat( "\nWarning:  QuSage needs at least 2 samples per group.")
		cat( "\nGroups with less than 2 samples:  ", bad)
		cat( "\nQuSage not run for:  ", contrast)
		return( NULL)
	}

	if (doQuSage) {

		ans <- qusage( eset=eset, labels=labels, contrast=contrast, geneSets=geneSets, n.points=n.points, ...)

		# modify and clean up the results
		pval <- pdf.pVal(ans)
		fdr <- p.adjust( pval, method="fdr")
		vif <- round( as.numeric( ans$vif), digits=3)
		pathSize <- as.numeric( ans$path.size)
		pathName <- names( ans$path.size)
		gNames <- names( ans$mean)
		pathGenes <- sapply( ans$pathway, function(x) paste( sort( unique( gNames[x])), collapse=";"))
		pathFold <- round( as.numeric( ans$path.mean), digits=3)
		pival <- round( piValue( pathFold, pval), digits=3)
	
		# we gave QuSage the short names, put the longer full names back
		where <- match( pathName, shortNames)
		pathName <- longNames[ where]
		if ( addCellTypes) {
			cellType <- getGeneSetCellType( pathName)
			nKeep <- 8
			out <- data.frame( "PATHWAY"=pathName, "CellType"=cellType, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "FDR"=fdr, "VIF"=vif, "PIVALUE"=pival, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		} else if ( addCellTypes) {
			lifeCycle <- getGeneSetLifeCycle( pathName)
			nKeep <- 8
			out <- data.frame( "PATHWAY"=pathName, "LifeCycle"=lifeCycle, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "FDR"=fdr, "VIF"=vif, "PIVALUE"=pival, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		} else {
			nKeep <- 7
			out <- data.frame( "PATHWAY"=pathName, "LOG2FOLD"=pathFold, 
					"PVALUE"=pval, "FDR"=fdr, "VIF"=vif, "PIVALUE"=pival, "N_GENES"=pathSize, 
					"GENE_LIST"=pathGenes, stringsAsFactors=FALSE)
		}
		ord <- diffExpressRankOrder( out$LOG2FOLD, out$PVALUE, wt.pvalue=2)
		out <- out[ ord, ]
		rownames(out) <- 1:nrow(out)
		outfile <- file.path( path, paste( group1, prefix, "UP.QuSage", descriptor, "csv", sep="."))
		write.table( out, outfile, sep=",", quote=T, row.names=F)
	} else {
		cat( "\nUsing previously calculated QuSage results..")
		outfile <- file.path( path, paste( group1, prefix, "UP.QuSage", descriptor, "csv", sep="."))
		out <- read.csv( outfile, as.is=T)
		nKeep <- if (addCellTypes) 8 else 7
	}

	otherGrpString <- paste( setdiff( sort( unique( labels)), group1), collapse=" + ")
	otherGrpString <- paste( "{", otherGrpString, "}", sep="")
	html1 <- subset( out, LOG2FOLD >= 0.1 & PVALUE <= 0.05)[ , 1:nKeep]
	if ( nrow(html1)) {
		htmlfile <- sub( "csv$", "html", outfile)
		htmltitle <- paste( "QuSage:  '", descriptor, "'  Pathways UP in: &nbsp; ", group1, 
					" &nbsp; vs &nbsp; ", otherGrpString, sep="")
		table2html( html1, htmlfile, title=addSpeciesToHtmlTitle(htmltitle), linkColumnNames=NULL)
	}
	html2 <- subset( out, LOG2FOLD <= -0.1 & PVALUE <= 0.05)[ , 1:nKeep]
	if ( makeDownHTML && nrow(html2)) {
		html2 <- html2[ rev( 1:nrow(html2)), ]
		rownames(html2) <- 1:nrow(html2)
		htmlfile <- sub( "csv$", "html", outfile)
		htmlfile <- sub( "UP", "DOWN", htmlfile)
		htmltitle <- paste( "QuSage:  '", descriptor, "'  Pathways DOWN in: &nbsp; ", group1, 
					" &nbsp; vs &nbsp; ", otherGrpString, sep="")
		table2html( html2, htmlfile, title=addSpeciesToHtmlTitle(htmltitle), linkColumnNames=NULL)
	}

	# we can also make an enrichment table of cell types
	if (addCellTypes) {
		out1 <- subset( out, LOG2FOLD >= 0.1 & PVALUE <= 0.05)
		if ( nrow(out1)) {
			enrich <- cellTypeEnrichment( out1$CellType, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( group1, prefix, "UP.QuSage.CellTypeEnrichment.csv", sep=".")
			f <- file.path( path, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
		out2 <- subset( out, LOG2FOLD <= -0.1 & PVALUE <= 0.05)
		if ( makeDownHTML && nrow(out2)) {
			enrich <- cellTypeEnrichment( out2$CellType, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( group1, prefix, "DOWN.QuSage.CellTypeEnrichment.csv", sep=".")
			f <- file.path( path, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
	}
	if (addLifeCycle) {
		out1 <- subset( out, LOG2FOLD >= 0.1 & PVALUE <= 0.05)
		if ( nrow(out1)) {
			enrich <- lifeCycleEnrichment( out1$LifeCycle, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( group1, prefix, "UP.QuSage.LifeCycleEnrichment.csv", sep=".")
			f <- file.path( path, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
		out2 <- subset( out, LOG2FOLD <= -0.1 & PVALUE <= 0.05)
		if ( makeDownHTML && nrow(out2)) {
			enrich <- lifeCycleEnrichment( out2$LifeCycle, upOnly=F, minEnrich=1, 
							maxPvalue=1, correct=T, verbose=F)
			f <- paste( group1, prefix, "DOWN.QuSage.LifeCycleEnrichment.csv", sep=".")
			f <- file.path( path, f)
			write.table( enrich, f, sep=",", quote=T, row.names=F)
		}
	}
	
	return( list( "QSarray"=ans, "Result"=out))
}

