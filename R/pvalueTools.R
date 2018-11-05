# pvalueTools.R	-- assorted tools for combining P values


`p.valid` <- function( x, min.p.value=1e-300, max.p.value=1.0) {

	# zeros and negative break most p-values
	# trap at smallest positive seen
	useX <- x
	useX[ x < 0] <- NA
	useX[ x < min.p.value] <- min.p.value
	useX[ x > max.p.value] <- max.p.value
	return( useX)
}


p.combine <- function( x) {

	x <- p.valid(x, max.p.value=0.99999)
	ansFisher <- p.combine.fisher( x, check.valid=F)
	ansStouffer <- p.combine.stouffer( x, check.valid=F)
	ansLogmean <- logmean( x, na.rm=T)
	return( max( c( ansLogmean, ansFisher, ansStouffer), na.rm=T))
}


p.combine.fisher <- function( x, check.valid=TRUE) {

	if (check.valid) x <- p.valid(x)
	return( pchisq( (sum( log(x)) * -2), df=length(x)*2, lower.tail=F))
}


p.combine.stouffer <- function( x, check.valid=TRUE) {

	if (check.valid) x <- p.valid(x, max.p.value=0.99999)
	return( pnorm( sum( qnorm( x)) / sqrt( length( x))))
}
