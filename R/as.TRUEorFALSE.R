# as.TRUEorFALSE.R

`as.TRUEorFALSE` <- function( str) {

	if ( length(str) > 1) {
		warning( "the condition has length > 1 and only the first element will be used")
		str <- str[1]
	}

	out <- NA

	# allow environment variable substitution
	str <- env.sub( str)
	if ( str == TRUE) out <- TRUE
	if ( str == FALSE) out <- FALSE

	# standardize whats there
	str <- toupper( gsub( " ", "", str))
	# obvious yes
	if ( str %in%  c( "TRUE", "T", "Y", "YES", "ON")) out <- TRUE
	# obvious no
	if ( str %in%  c( "FALSE", "F", "N", "NO", "OFF")) out <- FALSE

	# empty string or zero is false
	if ( str %in%  c( "", "0")) out <- FALSE

	# numeric?
	if ( is.numeric(str) && str > 0) out <- TRUE

	# last choice:  some text string -- call it true
	if ( is.na( out)) out <- TRUE
	return(out)
}
