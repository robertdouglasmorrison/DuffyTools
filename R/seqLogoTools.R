# seqLogoTools.R -- functions to implement SeqLogo type visuals of MSA Alignments


`ALNtoInformationContent` <- function( aln) {

	# we may be given the top level ALN object or just the aligment matrix
	# or even just the filename
	if ( length(aln) == 1 && is.character(aln) && file.exists(aln)) {
		aln <- readALN( aln, verbose=F)
	}
	if ( "alignment" %in% names(aln)) {
		aln <- aln$alignment
	}

	# OK, we have one matrix of characters from a MSA
	aln <- toupper(aln)
	nch <- ncol(aln)
	nseq <- nrow(aln)

	# deduce DNA from AA
	charSet <- sort( unique( as.vector(aln)))
	isDNA <- all( charSet %in% c("A","C","G","T","U","N","-"))
	N <- if ( isDNA) 4 else 20

	# tally up the various metrics
	# 1) the proportions of each base at each location
	pctM <- apply( aln, 2, function(x) table( factor(x,levels=charSet)) / length(x))

	# 2) the IC in bits, using both Wikipedia and ggseqlogo package (Omar Wagih)
	Hi <- - apply( pctM, 2, function(x) sum( x * log2(x), na.rm=T))
	en <- ( 1/log(2)) * ((N-1)/(2*nseq))
	Ric <- pmax( log2(N) - (Hi + en), 0)

	# 3) lastly, the heights of each base
	htM <- pctM
	for ( j in 1:nch) htM[ , j] <- pctM[ ,j] * Ric[j]

	out <- list( "proportion"=pctM, "info.content"=Ric, "height"=htM)
	return( out)
}


