# vsaTools.R -- small items to work with VSA gene names

`getVSAgeneMap` <- function() {

	vsaGeneMap <- NULL
	data( list="vsaGeneMap", package="DuffyTools",  envir=environment())
	if ( is.null( vsaGeneMap)) return( data.frame())
	return( vsaGeneMap)
}


`getVSAdomainMap` <- function() {

	vsaDomainMap <- NULL
	data( list="vsaDomainMap", package="DuffyTools",  envir=environment())
	if ( is.null( vsaDomainMap)) return( data.frame())
	return( vsaDomainMap)
}


`vsaTo3D7` <- function( genes) {

	vsaMap <- getVSAgeneMap()
	vsaMap$FULL_ID <- paste( vsaMap$STRAIN, vsaMap$GENE_ID, sep="::")

	out <- genes

	# the VSA genes are of the form STRAIN::ID
	isVSA <- union( which( genes %in% vsaMap$FULL_ID), which( genes %in% vsaMap$GENE_ID))

	for ( j in isVSA) {
		where <- match( genes[j], vsaMap$FULL_ID, nomatch=0)
		if ( where == 0) where <- match( genes[j], vsaMap$GENE_ID, nomatch=0)
		if ( where > 0) out[j] <- vsaMap$BEST_3D7_ID[ where]
	}

	return( out)
}


`vsaToOrig3D7` <- function( genes) {

	vsaMap <- getVSAgeneMap()
	vsaMap$FULL_ID <- paste( vsaMap$STRAIN, vsaMap$GENE_ID, sep="::")

	out <- genes

	# the VSA genes are of the form STRAIN::ID
	isVSA <- union( which( genes %in% vsaMap$FULL_ID), which( genes %in% vsaMap$GENE_ID))

	for ( j in isVSA) {
		where <- match( genes[j], vsaMap$FULL_ID, nomatch=0)
		if ( where == 0) where <- match( genes[j], vsaMap$GENE_ID, nomatch=0)
		if ( where > 0) out[j] <- vsaMap$BEST_ORIG_ID[ where]
	}

	return( out)
}


`vsaGeneProduct` <- function( genes) {

	prevSpecies <- getCurrentSpecies()
	if ( prevSpecies != "Pf3D7") {
		on.exit( setCurrentSpecies( prevSpecies))
		setCurrentSpecies( "Pf3D7")
	}

	vmap <- getVargeneDomainMap()
	vsamap <- getVSAgeneMap()

	vsaNames <- genes
	pf3d7Names <- vsaTo3D7( vsaNames)
	orig3d7Names <- vsaToOrig3D7( vsaNames)
	prods <- gene2Product( pf3d7Names)
	cassettes <- vsaCassettes( vsaNames)

	isVSA <- which( genes %in% vsamap$GENE_ID)
	is3D7var <- which( genes %in% c( vmap$GENE_NAME, vmap$ORIG_NAME))
	# let 3D7 take precedence...
	isVSA <- setdiff( isVSA, is3D7var)

	for ( j in isVSA) {
		group <- ""
		where <- match( pf3d7Names[j], vmap$GENE_NAME, nomatch=0)
		if ( where > 0) group <- paste( ";  Group:", vmap$GENE_GROUP[where])
		thisDC <- cassettes[j]
		if (thisDC != "") thisDC <- paste( "Cassettes: ", thisDC, ";", sep="")
		thisprod <- paste( "Variant PfEMP1 (", thisDC, " Similar to: ", pf3d7Names[j],
				" / ", orig3d7Names[j], group, ")", sep="")
		prods[j] <- thisprod
	}

	for ( j in is3D7var) {
		group <- ""
		where <- match( pf3d7Names[j], vmap$GENE_NAME, nomatch=0)
		if ( where > 0) group <- paste( "  Group:", vmap$GENE_GROUP[where])
		thisDC <- cassettes[j]
		if (thisDC != "") thisDC <- paste( "Cassettes: ", thisDC, ";", sep="")
		thisprod <- paste( prods[j], "  (", thisDC, group, ")", sep="")
		prods[j] <- thisprod
	}

	return( prods)
}


`domainCassetteDescriptor` <- function( dc) {

	prevSpecies <- getCurrentSpecies()
	if ( prevSpecies != "Pf3D7") {
		on.exit( setCurrentSpecies( prevSpecies))
		setCurrentSpecies( "Pf3D7")
	}
	vmap <- getVSAdomainMap()

	dc <- toupper( dc)
	N <- length(dc)
	out <- rep( "", times=N)
	for ( i in 1:N) {
		oneDC <- dc[i]
		smlmap <- subset( vmap, (CASSETTE != "" & CASSETTE == oneDC))
		if (nrow(smlmap) < 1) {
			cat( "\nUnknown VSA domain cassette name: ", oneDC)
			next
		}
		domains <- unique( smlmap$DOMAIN_CAT)
		grps <- sub( "^UPS_", "", unique( smlmap$UPS_CAT))
		grps <- sort( setdiff( grps, "ND"))
		domaintext <- paste( domains, collapse=",")
		grptext <- paste( grps, collapse=",")
		outtext <- paste( "Domains: {",domaintext, "} from Groups: {", grptext, "}", sep="")
		out[i] <- outtext
	}
	out
}


`vsaCassettes` <- function( vsaGenes) {

	# given a vector of VSA or 3D7 var genes
	vsaGenes <- gene2OrigID( vsaGenes)
	vdmap <- getVSAdomainMap()
	domHits <- sapply( vsaGenes, function(x) {
				who <- which( vdmap$GENE_ID == x)
				if ( ! length(who)) return( "")
				myDoms <- vdmap$CASSETTE[ who]
				drops <- which( myDoms %in% c( "", "ND", "none"))
				if ( length(drops)) {
					myDoms <- myDoms[ -drops]
					if ( ! length(myDoms)) return( "")
				}
				domSet <- sort( table( myDoms), decreasing=T)
				return( paste( names(domSet), collapse=","))
			})
	return( domHits)
}

	
`vsaDomainCassette` <- function( gene, AAloc=NULL, DNAloc=NULL) {

	# given a location inside a VSA gene, return the domain cassette name
	naDC <- ""
	vmap <- getVSAdomainMap()
	where <- match( gene, vmap$GENE_ID, nomatch=0)
	if (where == 0) {
		cat( "\nNot a known VSA gene: ", gene)
		return( naDC)
	}
	mode <- ""
	if ( ! is.null( AAloc)) {
		mode <- "AA"
		loc <- AAloc
	}
	if ( ! is.null( DNAloc)) {
		mode <- "DNA"
		loc <- DNAloc
	}
	if ( mode == "") {
		cat( "\nMust specify a location in the gene, as either 'AAloc' or 'DNAloc'")
		return( naDC)
	}

	# isolate that gene's domain edges, in coding order
	smlmap <- subset( vmap, GENE_ID == gene)
	ord <- order( smlmap$DNA_START)
	smlmap <- smlmap[ ord, ]
	if ( mode == "AA") {
		allStarts <- c( smlmap$AA_START, smlmap$AA_STOP[nrow(smlmap)])
	} else {
		allStarts <- c( smlmap$DNA_START, smlmap$DNA_STOP[nrow(smlmap)])
	}
	N <- length( allStarts)

	if ( loc < allStarts[1]) return( naDC)
	if ( loc > allStarts[N]) return( naDC)

	hit <- findInterval( loc, allStarts, all.inside=T)
	ans <- smlmap$CASSETTE[hit]
	if ( is.na( ans)) return( naDC)
	return(ans)
}


`domainCassetteEnrichment` <- function( domainIDs, upOnly=TRUE) {

	# given a vector of domains IDs (as from DE), give the enrichment of DCs
	vmap <- getVargeneDomainMap()
	where <- match( domainIDs, vmap$GENE_ID, nomatch=0)
	if ( any( where == 0)) cat( "\nSome domain names not in varGene domain map...")
	
	myDCs <- vmap$CASSETTE[ where]
	myTbl <- table( myDCs)
	if ( names(myTbl)[1] == "") names(myTbl)[1] <- "none"

	vTbl <- table( vmap$CASSETTE)
	if ( names(vTbl)[1] == "") names(vTbl)[1] <- "none"

	nDC <- length( vTbl)
	nDomains <- nrow(vmap)
	nGiven <- length(domainIDs)
	
	nhit <- phit <- enrich <- expect <- vector()
	for ( i in 1:nDC) {
		thisDC <- names(vTbl)[i]
		who <- match( thisDC, names(myTbl), nomatch=0)
		nhit[i] <- if (who > 0) myTbl[who] else 0
	
		ans <- enrichment( "nMatch"=nhit[i], "nYourSet"=nGiven, 
						"nTotal"=nDomains, "nTargetSubset"=vTbl[i])

		expect[i] <- ans$nExpect
		enrich[i] <- (nhit[i]/nGiven) / (vTbl[i]/nDomains)
		phit[i] <- min( ans$P_atLeast_N, ans$P_atMost_N)
	}

	out <- data.frame( "Cassette"=names(vTbl), "N_Total"=as.numeric(vTbl),
					"N_Given"=nhit, "Expected"=expect, "Enrichment"=enrich,
					"P_Value"=phit, stringsAsFactors=FALSE)
	ord <- order( out$Enrichment, -(out$P_Value), decreasing=T)
	out <- out[ ord, ]
	rownames(out) <- 1:nrow(out)

	if ( upOnly) {
		out <- subset( out, Enrichment > 1.0)
	}

	return( out)
}

