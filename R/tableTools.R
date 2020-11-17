# tableTools.R -- hopefully faster verison of counting 1-D table entries

table.nosort <- function(x) {

	s <- x[ ! duplicated(x)]
	wh <- match( x, s)
	cnt <- base::tabulate(wh)
	names(cnt) <- s
	cnt
}
