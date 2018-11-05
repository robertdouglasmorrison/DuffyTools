# levenshtein.R  -- edit distance tools

levenshtein <- function(string1, string2, case=TRUE, damerau=FALSE, map=NULL, motif=FALSE, fastMode=T) {
	
	########
	#
	# levenshtein algorithm in R
	#
	# Author  : Hans-Joerg Bibiko
	# Date    : 09/07/2006
	#
	# Contact : bibiko@eva.mpg.de
	#
	########
	#
	# string1, string2 := strings to compare
	#
	# case = TRUE := case sensitivity; case = FALSE := case insensitivity
	#
	# damerau = TRUE := Damerau enhancement: "al" to "la" costs only 1 instead of 2
	#
	# map := character vector of c(regexp1, replacement1, regexp2, replacement2, ...)
	#
	#   example:
	#      map <- c("[aeiou]","V","[^aeiou]","C") := replaces all vowels with V and all others with C
	#
	#      levenshtein("Bank","Bond", map=map)   =>  0
	#
	# motif:  if TRUE, then single replacement costs 2, remove and insert...
	#
	########
	
	
	if(!is.null(map)) {
		m <- matrix(map, ncol=2, byrow=TRUE)
		s <- c(ifelse(case, string1, tolower(string1)), ifelse(case, string2, tolower(string2)))
		for(i in 1:dim(m)[1]) s <- gsub(m[i,1], m[i,2], s)
		string1 <- s[1]
		string2 <- s[2]
	}
 	if  ( ! case) {
		string1 <- base::tolower( string1)
		string2 <- base::tolower( string2)
	}
	if( string1 == string2) return(0)

	# TRY THE FAST VERSION...
	if ( fastMode) return( fastLeven.C( string1, string2))

	chars <- strsplit( c( string1, string2), split="")
 	s1 <- chars[[1]]
 	s2 <- chars[[2]]
	rm( chars)

	# we can leave out tails that match exactly...
	mystart <- 1
	ncheck <- min( length(s1), length(s2))
	myREV <- base::rev
	while( mystart < ncheck && s1[mystart] == s2[mystart]) mystart <- mystart + 1
	ss1 <- myREV( s1[ mystart:length(s1)])
	ss2 <- myREV( s2[ mystart:length(s2)])
	mystart <- 1
	ncheck <- min( length(ss1), length(ss2))
	while( mystart < ncheck && ss1[mystart] == ss2[mystart]) mystart <- mystart + 1
	s1 <- c( " ", myREV( ss1[ mystart:length(ss1)]))
	s2 <- c( " ", myREV( ss2[ mystart:length(ss2)]))
	
	l1 <- length(s1)
	l2 <- length(s2)
	
	# the distance matrix
	d <- matrix(nrow = l1, ncol = l2)
	#for(i in 1:l1) d[i,1] <- i-1
	d[ ,1] <- 0:(l1-1)
	#for(i in 1:l2) d[1,i] <- i-1
	d[ 1, ] <- 0:(l2-1)

	# the equality of characters matrix...
	cost <- ifelse( motif, 2, 1)
	s1EQs2Cost <- matrix( cost, nrow = l1, ncol = l2)
	#for(i in 1:l1) {
	#for(j in 1:l2) {
		#if( s1[i] == s2[j]) s1EQs2Cost[i,j] <- 0
 	#}}
	for(i in 1:l1) {
		hits <- which( s1[i] == s2)
		if ( length(hits) > 0) s1EQs2Cost[ i, hits] <- 0
	}

	for(i in 2:l1) {
		im1 <- i - 1
		for(j in 2:l2) {
		jm1 <- j - 1
		# d[i,j] <- min((d[im1,j]+1) , (d[i,jm1]+1) , (d[im1,jm1]+ifelse(s1[i] == s2[j], 0, 1)))
	 	d[i,j] <- min((d[im1,j]+1) , (d[i,jm1]+1) , (d[im1,jm1]+ s1EQs2Cost[i,j]))
		#if(damerau && i>1 && j>1 && s1[i]==s2[j-1] && s1[i-1]==s2[j])
		#	d[i,j] <- min(d[i,j] , d[i-2,j-2]+ifelse(s1[i] == s2[j], 0, 1))
		#	d[i,j] <- min(d[i,j] , d[i-2,j-2]+ s1EQs2Cost[i,j])
	}}
	
	return( d[l1,l2])
}


levenshteinDistanceMatrix <- function( stringVector, case=TRUE, damerau=FALSE, map=NULL, motif=FALSE,
					diagonalValue=0, verbose=TRUE) {

	N <- length( stringVector)
	strs <- gsub( " ","", stringVector)
	m <- matrix( diagonalValue, nrow=N, ncol=N)

	if (verbose) cat( "\n")

	if ( multicore.currentCoreCount() == 1) {
		for ( i in 1:(N-1)) {
			if (verbose) cat( "\r", i, "of", N)
			for ( j in (i+1):N) {
				d <- levenshtein( strs[i], strs[j], case=case, damerau=damerau, map=map, 
						motif=motif, fastMode=TRUE)
				m[i,j] <- m[j,i] <- d
			}	
		}
	} else {
		for ( i in 1:(N-1)) {
			if (verbose) cat( "\r", i, "of", N)
			dVec <- multicore.lapply( strs[(i+1):N], FUN=levenshtein, strs[i], case=case, 
						damerau=damerau, map=map, motif=motif, fastMode=TRUE)
			for ( j in (i+1):N) {
				m[i,j] <- m[j,i] <- dVec[[ j-i]]
			}
		}
	}
	return( m)
}


`fastLeven.C` <- function( string1, string2){

	chInt1 <- utf8ToInt(string1)
	chInt2 <- utf8ToInt(string2)

	ans <- .C( "fastLevenshtein", 
			as.integer(chInt1), 
			as.integer(length(chInt1)), 
			as.integer(chInt2), 
			as.integer(length(chInt2)), 
			score=integer(1))
	return( ans$score[1])
}

