`objectSize` <- function( objects=ls( envir=.GlobalEnv), N=length(objects), 
				sortOrder=c("size", "name", "none")) {

	sizes <- sapply( objects, FUN=function(x) object.size( getAnywhere(x)))

	sortOrder <- match.arg( sortOrder)
	if( sortOrder == "name") {
		ord <- order( names(sizes))
		sizes <- sizes[ord]
	}
	if( sortOrder == "size") {
		ord <- order( sizes, decreasing=TRUE)
		sizes <- sizes[ord]
	}
	
	if ( N < length(sizes)) sizes <- sizes[ 1:N]
	sizes
}
