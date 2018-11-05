# as.rankPercentile.R

`as.rankPercentile` <- function( x) {

	N <- length(x)
	isNA <- is.na(x)
	NnotNA <- sum( ! isNA)
	ranks <- rank( x[ !isNA], na.last=T, ties.method="average")
	rankPct <- vector( mode="numeric", length=N)
	rankPct[ isNA] <- NA
	rankPct[ !isNA] <- ( ranks * 100 / NnotNA)
	return( rankPct)
}
