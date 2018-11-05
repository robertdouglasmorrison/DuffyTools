`cgFraction` <- function(dna) {

	# get the first element of the list returned from "strsplit"
	letters <- strsplit( toupper( dna), split="", fixed=TRUE)
	out <- sapply( letters, function(x) {
				at <- sum( x %in% c( "A","T"))
				cg <- sum( x %in% c( "C","G"))
				return( cg / (at+cg))
			})
	names(out) <- names(dna)
	out
}

