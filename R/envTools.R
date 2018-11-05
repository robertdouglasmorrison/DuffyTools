# envTools.R

#  little utilities to do environment variable substitutions in character strings

env.sub <- function( x) {

	out <- x
	for ( i in 1:length(x)) {

		thisX <- x[ i]
		dollarsigns <- gregexpr( "$", thisX, fixed=T)[[1]]
		nDollars <- length( dollarsigns)
		if ( nDollars == 0) next

		# at least one potential substitution
		newstr <- base::substr( thisX, 1, (dollarsigns[1] - 1))
		for ( j in 1:nDollars) {
			start <- dollarsigns[j] + 1
			if ( base::substr( thisX, start, start) == "{") {
				start <- start + 1
				endPattern <- "}"
				dropFront <- 2
				dropRear <- 1
			} else {
				endPattern <- c("/"," ", "$")
				dropFront <- 1
				dropRear <- 0
			}
			# allow a zero length insert
			end <- start - 1
			repeat {
				if ( end == base::nchar( thisX)) break
				if ( base::substr( thisX, end+1, end+1) %in% endPattern) break
				end <- end + 1
			}
			# when we break on a "}", consume it
			nextStart <- end + 1 + dropRear

			# ask System about this posible var
			varname <- base::substr( thisX, start, end)
			ans <- Sys.getenv( varname)
			if ( ans != "") {
				# got a hit
				varvalue <- ans
				newstr <- base::paste( newstr, varvalue, sep="")
			} else {
				# if it was {...} style, drop it all...  otherwise leave it all...
				if ( dropRear == 0)
				newstr <- base::paste( newstr, base::substr( thisX, start-1, nextStart-1), sep="")
			}
			if ( j < nDollars) {
				newstr <- base::paste( newstr, base::substr( thisX, nextStart, (dollarsigns[j+1] - 1)), sep="")
			}
		}

		# append the last bit...
		newstr <- base::paste( newstr, base::substr( thisX, nextStart, base::nchar( thisX)), sep="")
		out[i] <- newstr
	}
	
	return( out)
}
