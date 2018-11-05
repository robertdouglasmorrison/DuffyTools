# fastFindInterval.R

#  bare bones version of 'findInterval()' for speed
`fastFindInterval` <- function( x, vec) {

	majorVersion <- as.numeric( R.version$major)
	minorVersion <- as.numeric( R.version$minor)

	if ( majorVersion < 3) {
		nx <- length(x)
		index <- integer(nx)
		.C( "find_interv_vec", nt=as.double(vec), n=length(vec), x=as.double(x), nx=nx, FALSE, FALSE,
				index, DUP=FALSE, NAOK=TRUE, PACKAGE="base")
		return( index)
	}

	if ( majorVersion == 3) {
		if ( minorVersion >= 4.0) {
   			return( .Internal( findInterval( as.double(vec), as.double(x), rightmost.closed=FALSE, 
       		 		all.inside=FALSE, left.open=FALSE)))
		} else {
   			return( .Internal( findInterval( as.double(vec), as.double(x), rightmost.closed=FALSE, 
       		 		all.inside=FALSE)))
		}
	}

	warning( "No implementation of 'fastFindInterval' for R version", majorVersion, minorVersion)
	return( findInterval( x, vec))
}
