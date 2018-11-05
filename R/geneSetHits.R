# geneSetHits.R -- make a table of top gene set names from a set of genes


geneSetHits <- function( genes, geneSets=c("GO.BiologicalProcess", "GO.MolecularFunction",
				"GO.CellularComponent", "KEGG.Pathways", "MetabolicPathways", "GeneProduct",
				"PBMC.GeneModules", "Blood.GeneModules", "MPMP.Pathways")) {
	
	# the 'geneSets' can be a vector of strings of geneSet datasets..., 
	# or the usual list of vectors of genes
	cat( "\nGathering GeneSets..")
	GSanswer <- gatherGeneSets( geneSets, descriptor="geneSets", mode="combined")
	geneSets <- GSanswer$geneSets[[1]]

	# we don't want any HTML links for this step...
	names(geneSets) <- trimGeneSetNameLink( names( geneSets))


	getPathwaySetData <- function( pathwaySets) {

		# given the list of all gene sets
		setNames <- names( pathwaySets)
		cnts <- sapply( pathwaySets, length)
		allGenes <- unlist( pathwaySets)
		allNames <- rep( setNames, times=cnts)
		pathTable <- data.frame( "GeneID"=allGenes, "PathwayID"= allNames, stringsAsFactors=F)
		rownames(pathTable) <- 1:nrow(pathTable)

		nGenesInTable <- length( unique( pathTable$GeneID))
		pathID.table <- table( pathTable$PathwayID)
		pathGenePcts <- pathID.table * 100 / nGenesInTable
		pathGeneCnts <- pathID.table
	
		cat( "\nFound ", length(pathID.table), " Gene Sets to evaluate")
		out <- list( "nRowsTotal"=nrow( pathTable), "nUniqueGenes"=nGenesInTable, 
				"PathwayTable"=pathTable, "idTable"=pathID.table, 
				"pctByGene"=pathGenePcts, "countByGene"=pathGeneCnts)
		return( out)
	}


	getGeneSetData <- function( genes, pathTable) {
	
		mySet <- subset( pathTable, GeneID %in% genes)
		if ( nrow(mySet) < 1) {
			cat( "\nWarning:  None of the given genes are found in any GeneSet")
			cat( "\nVerify the current species matches your genes.\n")
			out <- list( "nRowsTotal"=0, "nUniqueGenes"=0, 
					"PathwayTable"=data.frame(), "idTable"=vector(), 
					"pctByGene"=vector(), "countByGene"=vector())
			return( out)
		}
		rownames(mySet) <- 1:nrow(mySet)
	
		nGenesInSet <- length( unique( mySet$GeneID))
		pathID.table <- table( mySet$PathwayID)
		pathGenePcts <- pathID.table * 100 / nGenesInSet
		pathGeneCnts <- pathID.table
	
		out <- list( "nRowsTotal"=nrow( mySet), "nUniqueGenes"=nGenesInSet, 
				"PathwayTable"=mySet, "idTable"=pathID.table, 
				"pctByGene"=pathGenePcts, "countByGene"=pathGeneCnts)
		return( out)
	
	}
	
	
	setData <- getPathwaySetData( geneSets)
	geneData <- getGeneSetData( genes, setData$PathwayTable)

	return( geneData)
}


geneSetHitsToPie <- function( hits, min.genes=1, max.show=20, label="", txt.padX = 0, txt.padY = 0 ) {

	# take the output of 'geneSetHits()' and render as a Pie chart
	countData <- hits$countByGene
	keep <- which( countData >= min.genes)
	countData <- countData[ keep]
	ord <- order( countData, decreasing=T)
	countData <- countData[ ord]
	if ( length(countData) > max.show) length(countData) <- max.show

	txt.padX <- rep( txt.padX, length.out=length(countData))
	txt.padY <- rep( txt.padY, length.out=length(countData))
	
	mycolors <- rainbow( length(countData), end=0.9)

	# append a percentage to the names...
	nams <- names(countData)
	pcts <- as.percent( countData, big.value=sum(countData), digits=0)
	for ( j in 1:length(nams)) {
		nc <- nchar( nams[j])
		pad <- paste( rep.int( " ", round(nc/2)-2), collapse="")
		nams[j] <- paste( nams[j], "\n", pad, pcts[j], pad, sep="")
	}
	names(countData) <- nams

	geneSetPie( countData, col=mycolors, radius=0.7, init.angle=90, edges=300, border='white', 
			main=label, txt.padX=txt.padX, txt.padY=txt.padY)

}


geneSetPie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
			init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
			col = NULL, border = NULL, lty = NULL, main = NULL, txt.padX = 0, txt.padY = 0) 
{
	if (!is.numeric(x) || any(is.na(x) | x < 0)) 
	    stop("'x' values must be positive.")
	if (is.null(labels)) 
	    labels <- as.character(seq_along(x))
	else labels <- as.graphicsAnnot(labels)

	x <- c(0, cumsum(x)/sum(x))
	dx <- diff(x)
	nx <- length(dx)
	plot.new()
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1)
	if (pin[1L] > pin[2L]) 
	    xlim <- (pin[1L]/pin[2L]) * xlim
	else ylim <- (pin[2L]/pin[1L]) * ylim
	dev.hold()
	on.exit(dev.flush())
	plot.window(xlim, ylim, "", asp = 1)
	if (is.null(col)) 
	    col <- if (is.null(density)) 
	        c("white", "lightblue", "mistyrose", "lightcyan", 
	            "lavender", "cornsilk")
	    else par("fg")
	col <- rep(col, length.out = nx)
	border <- rep(border, length.out = nx)
	lty <- rep(lty, length.out = nx)
	angle <- rep(angle, length.out = nx)
	density <- rep(density, length.out = nx)
	twopi <- if (clockwise) 
	    -2 * pi
	else 2 * pi
	t2xy <- function(t) {
	    t2p <- twopi * t + init.angle * pi/180
	    list(x = radius * cos(t2p), y = radius * sin(t2p))
	}
	for (i in 1L:nx) {
	    n <- max(2, floor(edges * dx[i]))
	    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
	    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
	        border = border[i], col = col[i], lty = lty[i])
	    P <- t2xy(mean(x[i + 0:1]))
	    lab <- as.character(labels[i])
	    if (!is.na(lab) && nzchar(lab)) {
	        lines( (c(1, 1.15+txt.padX[i]) * P$x), (c(1, 1.25+txt.padY[i]) * P$y))
	        text( (1.2 * P$x) + txt.padX[i], (1.35 * P$y) + txt.padY[i], labels[i], xpd = TRUE, 
                adj = if (P$x < -0.15) 1 else if (P$x > 0.15) 0 else 0.5, font=2)
	    }
	}
	title(main = main)
	invisible(NULL)
}

