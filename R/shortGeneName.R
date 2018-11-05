# shortGeneName.R

shortGeneName <- function( genes, split=":", keep=2) {

	if ( keep == 1) {
		return( base::sub( base::paste( split,".+",sep=""), "", genes))
	}

	terms <- strsplit( genes, split=split, fixed=TRUE)
	out <- base::sapply( terms, function(x) { 
			n <- min( c(length(x), keep)); 
			base::paste( x[1:n], collapse=split)
		})
	return( out)
}
