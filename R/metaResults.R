# metaResults.R -- combine several ranked files of gene results

 
metaResults <- function( targetGroup, results.path="results", speciesID=getCurrentSpecies(), 
			geneColumn="GENE_ID", subfolderName="Group", 
			tools=c("RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), altGeneMapLabel=NULL,
			rank.average.FUN=sqrtmean, value.average.FUN=mean, keepIntergenics=FALSE,
			topFolders=NULL, other.DE.files=NULL, missingGenes=c("drop", "fill", "na"), nFDRsimulations=0,
			direction=c("UP", "DOWN"), otherGroup=paste("Not",targetGroup,sep=".")) {

	missingGenes <- match.arg( missingGenes)
	direction <- match.arg( direction)

	setCurrentSpecies( speciesID)
	prefix <- getCurrentSpeciesFilePrefix()

	cat( "\nChecking file names...")
	allFiles <- allNames <- vector()
	nfiles <- 0

	if ( is.null(topFolders)) {
		nFolders <- 1
	} else {
		nFolders <- length( topFolders)
	}
	for (ifolder in 1:nFolders) {
	for (tool in tools) {
		suffix <- NULL
		if( tool == "RoundRobin") suffix <- "RR" 
		if( tool == "RankProduct") suffix <- "RP" 
		if( tool == "SAM") suffix <- "SAM" 
		if( tool == "EdgeR") suffix <- "EdgeR" 
		if( tool == "DESeq") suffix <- "DESeq" 
		if ( is.null(topFolders)) {
			mypath <- file.path( results.path, tool)
		} else {
			mypath <- file.path( topFolders[ifolder], results.path, tool)
		}
		mysubfolder <- paste( prefix, subfolderName, sep=".")
		myfile <- paste( targetGroup, prefix, suffix, "Ratio.txt", sep=".")
		if ( ! is.null( altGeneMapLabel)) {
			myfile <- paste( targetGroup, prefix, altGeneMapLabel, suffix, "Ratio.txt", sep=".")
		}

		filename <- file.path( mypath, mysubfolder, myfile)
		if ( ! file.exists( filename)) {
			cat( "\nFile not found:  ", filename)
			next
		}
		nfiles <- nfiles + 1
		allFiles[ nfiles] <- filename
		allNames[ nfiles] <- tool
		if ( nFolders > 1) allNames[ nfiles] <- paste( names(topFolders)[ifolder], tool, sep=" ")
	}}

	if ( ! is.null( other.DE.files)) {
		for ( j in 1:length( other.DE.files)) {
			filename <- other.DE.files[j]
			if ( ! file.exists( filename)) {
				cat( "\nFile not found:  ", filename)
				next
			}
			nfiles <- nfiles + 1
			allFiles[ nfiles] <- filename
			allNames[ nfiles] <- names(other.DE.files)[j]
		}
	}
	allNames[ is.null(allNames)] <- 'unnamed file'
	if ( length( allFiles) < 2) {
		cat( "\nLess than 2 files found.  Not enough files to do 'metaResults'...")
		return( data.frame())
	}

	ans <- metaRanks( allFiles, allNames, geneColumn=geneColumn, rank.average.FUN=rank.average.FUN,
			value.average.FUN=value.average.FUN, keepIntergenics=keepIntergenics, 
			missingGenes=missingGenes, nFDRsimulations=nFDRsimulations,
			diffExpressPvalueCorrection=TRUE, direction=direction)

	# append the RPKM data for the two conditions
	ans <- appendMetaResultGroupValues( ans, targetGroup, files=allFiles, tools=allNames, geneColumn=geneColumn,
					otherGroup=otherGroup)

	return( ans)
}


`appendMetaResultGroupValues` <- function( tbl, targetGroup, files, tools, geneColumn="GENE_ID", 
					otherGroup=paste("Not",targetGroup,sep=".")) {

	# several (but not all) DE tools give their call for gene expression.
	# grab those to supplement the results
	value1 <- value2 <- gene <- vector()
	nNow <- 0

	# accumulate for all tools with gene expression data
	for (j in 1:length(files)) {
		thisTool <- tools[j]
		if ( ! (thisTool %in% c( "RoundRobin", "RankProduct", "SAM"))) next

		smlDF <- read.delim( files[j], as.is=T)
		N <- nrow(smlDF)

		# make all the tools use same naming...
		#if ( thisTool == "RoundRobin") {
		#	name1 <- "RPKM_1"
		#	name2 <- "RPKM_2"
		#	if ( length( grep( "RPKM", colnames(smlDF))) < 1) {
		#		name1 <- "VALUE_1"
		#		name2 <- "VALUE_2"
		#	}
		#} else if ( thisTool %in% c("RankProduct", "SAM")) {
		name1 <- targetGroup
		name2 <- otherGroup
		# if the group name starts with a digit, the column will start with an 'X'
		if ( substr(name1,1,1) %in% as.character( 0:9)) name1 <- paste( 'X', name1, sep="")
		if ( substr(name2,1,1) %in% as.character( 0:9)) name2 <- paste( 'X', name2, sep="")
		# also, certain characters will not be preserved in an R table column name
		name1 <- make.names(name1)
		name2 <- make.names(name2)

		if ( ! all( c( geneColumn, name1, name2) %in% colnames(smlDF))) {
			cat( "\nSome needed VALUE columns not found.  Wanted: ", geneColumn, name1, name2)
			cat( "\nFound:  ", colnames(smlDF))
			cat( "\nSkipping..")
			next
		}
		this1 <- smlDF[[name1]]
		this2 <- smlDF[[name2]]
		thisGene <- smlDF[[geneColumn]]

		# simple trap for negative expression
		if ( any( c(this1,this2) < 0)) {
			who <- which( pmin( this1, this2) < 0)
			cat( "\nTrap negative expression: ", thisTool, "  N =", length(who), "\n")
			tmpDF <- data.frame( "GeneID"=thisGene[who], "ValueSelf"=this1[who], "ValueOther"=this2[who])
			print( head( tmpDF, 10))
			this1[ this1 < 0] <- 0
			this2[ this2 < 0] <- 0
		}

		now <- (nNow+1) : (nNow+N)
		value1[now] <- this1
		value2[now] <- this2
		gene[now] <- thisGene
		nNow <- max( now)
	}

	# use the average of all found tools
	geneFac <- factor( gene)
	geneNames <- levels(geneFac)
	value1 <- tapply( value1, geneFac, FUN=sqrtmean, na.rm=T)
	value2 <- tapply( value2, geneFac, FUN=sqrtmean, na.rm=T)

	# stuff these where they go
	where <- match( geneNames, tbl[[geneColumn]], nomatch=0)
	outVALUE <- matrix( NA, nrow=nrow(tbl), ncol=2)
	colnames(outVALUE) <- c( targetGroup, otherGroup)
	outVALUE[ where, 1] <- value1[where > 0]
	outVALUE[ where, 2] <- value2[where > 0]
	
	out <- cbind( tbl, round( outVALUE, digits=2), stringsAsFactors=FALSE)
	return( out)
}


`metaResultsToHTML` <- function(m, fileout, title="", maxRows=100, 
				linkColumnNames=NULL, ...) {

	if ("PRODUCT" %in% colnames(m)) m$PRODUCT <- gsub( "   ", " &nbsp; ", m$PRODUCT)
	m$LOG2FOLD <- formatC( m$LOG2FOLD, format="f", digits=3, flag="+")
	m$AVG_PVALUE <- formatC( m$AVG_PVALUE, format="e", digits=2)
	m$AVG_RANK <- formatC( m$AVG_RANK, format="f", digits=1)
	colnames(m)[ match( "LOG2FOLD",colnames(m))] <- "Log2 Fold"
	colnames(m)[ match( "AVG_PVALUE",colnames(m))] <- "Avg Pvalue"
	colnames(m)[ match( "AVG_RANK",colnames(m))] <- "Avg Rank"
	colnames(m) <- sub( "RoundRobin", "Round Robin", colnames(m))
	colnames(m) <- sub( "RankProduct", "Rank Prod", colnames(m))
	colnames(m) <- sub( "EdgeR", "Edge R", colnames(m))
	colnames(m) <- sub( "DESeq", "DE Seq", colnames(m))
	if ( "E_VALUE" %in% colnames(m)) {
		m$E_VALUE <- formatC( m$E_VALUE, format="f", digits=4)
		m$FP_RATE <- formatC( m$FP_RATE, format="f", digits=4)
		colnames(m) <- sub( "FP_RATE", "FDR", colnames(m))
		colnames(m) <- sub( "E_VALUE", "E Value", colnames(m))
	}

	NC <- ncol(m)
	if ( NC > 5) colnames(m)[6:NC] <- gsub( "_", " ", colnames(m)[6:NC], fixed=T)
	if ( NC > 5) colnames(m)[6:NC] <- gsub( ".", " ", colnames(m)[6:NC], fixed=T)

	table2html( m, fileout, title=title, maxRows=maxRows, linkColumnNames=linkColumnNames,
			...) 
	cat( "\nWrote file:  ", fileout)

	return()
}


`plotMetaResults` <- function( file, columnNames=c( "RoundRobin", "RankProduct", "SAM", "EdgeR", "DESeq"), 
				columnColors=rainbow( length(columnNames), end=0.72),
				label="") {

	if ( ! file.exists( file)) {
		cat("\nMeta Result file not found: ", file)
		return()
	}
	tbl <- read.delim( file, as.is=T)
	if ( nrow(tbl) < 2) {
		cat( "\nNot enough rows in Meta Result file")
		return()
	}

	# find those rank columns
	colPtr <- match( columnNames, colnames(tbl), nomatch=0)
	if ( any( colPtr == 0)) {
		cat( "\nSome Meta Result columns names not found:  ", columnNames[ colPtr == 0])
	}
	columnNames <- columnNames[ colPtr > 0]
	columnColors <- columnColors[ colPtr > 0]
	colPtr <- colPtr[ colPtr > 0]
	if ( length( colPtr) < 2) {
		cat( "\nNot enough Meta Result columns found to plot")
		return()
	}

	# gather the found data
	NC <- length( colPtr)
	NR <- nrow(tbl)
	m <- matrix( NA, nrow=NR, ncol=NC)
	colnames(m) <- columnNames
	for ( i in 1:NC) m[ , i] <- as.numeric( tbl[[ colPtr[i]]])

	# now turn every rank into its deviation from the consensus
	lapply( 1:NR, function(i) { m[ i, ] <<- i - m[ i, ] })

	# now set up to draw
	ylim <- range( m, na.rm=T)
	# be symmetric about the smaller limit
	ylim <- c(-1,1) * min( abs( ylim))
	par( mai=c(1,1,0.8,0.4))

	mainText <- paste( "Meta Result Deviation Plot:  ", label, "\n", basename(file))
	plot( 1,1, type="n", main=mainText, xlim=c(1,NR), ylim=ylim, 
			xlab="Consensus Gene Rank", ylab="(Worse)     Deviation from Consensus Rank     (Better)",
			font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1)

	# pre-make the way we will draw as line segments
	y1 <- col <- rep.int( 0, NC)

	lapply( 1:NR, function(i) {
			v <- m[ i, ]
			# we will draw from longest to shortest, so the colors overlay right
			ord <- order( abs(v), decreasing=T)
			# make it into line segments
			x0 <- x1 <- rep.int( i, NC)
			y0 <- v[ord]
			col <- columnColors[ord]
			segments( x0, y0, x1, y1, col=col)
		})

	legend( "topleft", columnNames, fill=columnColors, cex=1.2, bg="white")
	printPlot( "MetaResultDeviation.Lineplot")
	Sys.sleep(5)


	# also as a box plot to show the average deviation
	ylim <- c( 0, max(m, na.rm=T)*0.65)
	boxplot( abs(m), col=columnColors, main=mainText, 
			ylab="(Better)       Absolute Deviation from Consensus     (Worse)",
			ylim=ylim, pars=list( font.axis=2, font.lab=2, cex.axis=1.1, cex.lab=1.1))
	printPlot( "MetaResultDeviation.Boxplot")

	return()
}

