# tableTools.R -- hopefully faster verison of counting 1-D table entries

table.nosort <- function(x) {

	s <- x[ ! base::duplicated.default(x)]
	wh <- base::match( x, s)
	cnt <- base::tabulate(wh)
	names(cnt) <- s
	cnt
}
