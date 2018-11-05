# oneDtable.R  - a quicker 1-D table function

oneDtable <- function( x) {

	#levels <- .Internal( unique( x, FALSE, FALSE))
	levels <- unique.default( x)
	levels <- levels[ sort.list( levels)]
	f <- .Internal( match( x, levels, NA_integer_, NULL))
	tbl <- tabulate( f, nbins=length(levels))
	names(tbl) <- levels
	if ( any( is.na( levels))) length(tbl) <- length(tbl) - 1
	class(tbl) <- "table"
	tbl
}
