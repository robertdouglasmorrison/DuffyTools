# AA3mers.R -- tools to work with amino acid 3-mer frequencies



`getAA3merFreqMap` <- function() {

	prefix <- getCurrentSpeciesFilePrefix()
	f <- paste( prefix, "AA3merFreqMap", sep=".")
	AA3merMap <- NULL
	data( list=f, envir=environment()) 
	if ( is.null( AA3merMap)) return(NULL)
	return( AA3merMap)
}


as.AA3mers <- function( protein) {

	PASTE <- base::paste
	MAPPLY <- base::mapply

	proteinV <- strsplit( toupper(protein), split="")[[1]]
	nc <- length( proteinV)
	if ( nc < 3) return("")

	beg <- 1:(nc-2)
	end <- beg + 2
	N3mers <- length(beg)
	triplets <- MAPPLY( beg, end, FUN=function(b,e) PASTE( proteinV[b:e], collapse=""),
			USE.NAMES=FALSE)
	triplets
}


AA3merFreq <- function( protein, AA3merFreqMap=NULL) {

	my3mers <- as.AA3mers( protein)
	if ( is.null(AA3merFreqMap)) AA3merFreqMap <- getAA3merFreqMap()

	myFreq <- rep.int( 0, length(my3mers))
	where <- match( my3mers, AA3merFreqMap$AA3MER, nomatch=0)
	myFreq[ where > 0] <- AA3merFreqMap$FREQUENCY[where]

	# make it be the same length as the protein
	out <- c( myFreq[1], myFreq, myFreq[length(myFreq)])
	out
}


`buildAA3merFreqMap` <- function( fastaFile, speciesID=getCurrentSpecies()) {

	MATCH <- base::match
	TABLE <- base::table

	fa <- loadFasta( fastaFile)
	data( list="PAM70MS", envir=environment())
	aaSet <- sort( unique( setdiff( colnames(PAM70MS), c("B","J","U","X","Z","*"))))
	nAA <- length(aaSet)
	aa1 <- rep( aaSet, each=(nAA^2))
	aa2 <- rep( rep(aaSet, each=nAA), times=nAA)
	aa3 <- rep( aaSet, times=(nAA^2))
	aa3mer <- paste( aa1, aa2, aa3, sep="")

	counts <- rep.int( 0, length(aa3mer))

	for ( i in 1:length( fa$desc)) {
		id <- fa$desc[i]
		seq <- fa$seq[i]
		my3mers <- as.AA3mers( seq)
		tbl <- TABLE( my3mers)
		where <- MATCH( names(tbl), aa3mer, nomatch=0)
		counts[ where] <- counts[ where] + tbl[where > 0]
		cat( "\r", i, id, sum(tbl), " ")
	}

	freqs <- counts * 100 / sum(counts)
	out <- data.frame( "AA3MER"=aa3mer, "COUNT"=counts, "FREQUENCY"=freqs,
			stringsAsFactors=F)
	out
}
