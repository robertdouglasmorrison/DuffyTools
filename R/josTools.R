# josTools.R -- small items to work with Jose's novel var gene names

`getJOSgeneMap` <- function() {

	mapSet <- NULL
	data( "JOS.MapSet",  package="DuffyTools",  envir=environment())
	if ( is.null( mapSet)) return( data.frame())
	return( mapSet$geneMap)
}


`getJOSdomainMap` <- function() {

	josDomainMap <- NULL
	data( "josDomainMap",  package="DuffyTools",  envir=environment())
	if ( is.null( josDomainMap)) return( data.frame())
	return( josDomainMap)
}


`josGeneProduct` <- function( genes) {

	gmap <- getJOSgeneMap()

	prods <- gene2Product( genes)

	where <- match( genes, gmap$GENE_ID, nomatch=0)
	prods[ where > 0] <- gmap$PRODUCT[ where]
	return( prods)
}


`josCassettes` <- function( josGenes) {

	# given a vector of JOS var genes
	vdmap <- getJOSdomainMap()
	domHits <- sapply( josGenes, function(x) {
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

	

`josDomainCassette` <- function( gene, AAloc=NULL, DNAloc=NULL) {

	# given a location inside a JOS gene, return the domain cassette name
	naDC <- ""
	vmap <- getJOSdomainMap()
	where <- match( gene[1], vmap$GENE_ID, nomatch=0)
	if (where == 0) {
		cat( "\nNot a known JOS gene: ", gene)
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

