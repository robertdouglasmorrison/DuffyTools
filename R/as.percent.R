# as.percent.R

`as.percent` <- function( x, digits=2, big.value=1.0, percentSign=TRUE, sep=" ") {

	scaleFactor <- 100 / big.value[1]
	x <- unclass(x)
	out <- base::paste( formatC( (scaleFactor * x), width=(4+digits), digits=digits, 
			format="f"))
	if( percentSign) out <- base::paste( out, "%", sep=sep)
	if ( ! is.null( names(x))) {
		names(out) <- names(x)
	}

	return( out)
}

`test.as.percent` <- function() {
	checkEquals( as.percent( 0.01), "  1.00 %")
	checkEquals( as.percent( 21, digits=0, big.value=42), "  50 %")
}
