`repeatFraction` <-
function(dna) {
	# get the first element of the list returned from "strsplit"
	letters <- (strsplit(dna[1], split="", fixed=TRUE))[[1]]
	n <- length(letters)
	nRepeat <- sum( letters[1:(n-1)] == letters[ 2:n])
	return( as.double(nRepeat) / as.double(n) )
}

