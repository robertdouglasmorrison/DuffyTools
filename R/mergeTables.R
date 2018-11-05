# mergeTables.R   - combine two 1-D tables that may or may not have common categories


`mergeTables` <- function( t1, t2) {

	nam1 <- names(t1);   nam2 <- names(t2);
	#allNames <- .Internal( unique( c( nam1, nam2), FALSE, FALSE))
	allNames <- unique.default( c( nam1, nam2))
	ord <- .Internal( order( TRUE, FALSE, allNames))
	allNames <- allNames[ ord]
	len <- length( allNames)
	allCounts <- rep( 0, times=len)
	where1 <- .Internal( match( nam1, allNames, NA, NULL))
	allCounts[ where1] <- allCounts[ where1] + t1
	where2 <- .Internal( match( nam2, allNames, NA, NULL))
	allCounts[ where2] <- allCounts[ where2] + t2
	names(allCounts) <- allNames
	if ( any( is.na( allNames))) length(allCounts) <- length(allCounts) - 1
	class(allCounts) <- "table"
	allCounts
}


`mergeIntegerTables` <- function( t1, t2) {

	allNames <- sort.int( as.integer( unique.default( c( names(t1), names(t2)))))
	len <- length( allNames)
	allCounts <- rep( 0, times=len)
	where1 <- base::match( as.integer(names(t1)), allNames, nomatch=0)
	allCounts[ where1] <- allCounts[ where1] + t1
	where2 <- base::match( as.integer(names(t2)), allNames, nomatch=0)
	allCounts[ where2] <- allCounts[ where2] + t2

	dim(allCounts) <- len
	dimnames(allCounts) <- list( allNames)
	class(allCounts) <- "table"

	return( allCounts)
}

