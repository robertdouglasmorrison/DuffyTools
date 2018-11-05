# phredScoreTools.R

# assess read qualities

`getMaxPhredScore` <- function() { return( 60)}


`solexaCharScoreToSolexaIntScore` <- function( str) {

	# turn the compact Solexa character string of scores into their older string of integer scores
	base64 <- 64
	return( base::paste( as.character( utf8ToInt( str) - base64), collapse=" "))
}


`solexaToPhred` <- function( str, scoreType="Phred33") {

	# turn the text string of integer base quality scores into the Phred-type compact ASCII string
	str <- base::sub( "^ +", "", str)
	#while(  base::substr( str, 1,1) == " ") str <- base::substr( str, 2, base::nchar(str))
	scoretmp <- strsplit( str, split=" ", fixed=TRUE)[[1]]

	# turn integers into Phred characters
	phredScores <- solexaScore2phredScore( as.numeric( scoretmp))

	# turn these into a character string where 0 == " "
	offSet <- 33
	if (scoreType == "Phred64") offSet <- 64
	outStr <- intToUtf8( as.integer( phredScores + offSet))
	return( outStr)
}


`oldPhredScoreStringToInt` <- function( txt, scoreType="Phred33") {

	offSet <- 33
	if (scoreType == "Phred64") offSet <- 64

	# turn each cryptic ASCII string of scores into a vector of integers...	
	fac <- factor( txt)
	uniqueStrs <- levels(fac)
	facPtrs <- tapply( 1:length(txt), fac, FUN=NULL)
	ncSet <- base::nchar( uniqueStrs)
	cat( " PhredInt Speedup=", formatC( (length(txt)/length(uniqueStrs)), digits=2, format="E"))

	tmpOut <- matrix( 0, nrow=length(uniqueStrs), ncol=max(ncSet))
	for( i in 1:length( uniqueStrs)) {
		tmpOut[ i, (1:ncSet[i])] <- utf8ToInt( uniqueStrs[i]) - offSet
	}

	out <- tmpOut[ facPtrs, ]

	return( out)
}


`phredScoreStringToInt` <- function( txt, scoreType="Phred33") {

	offSet <- 33
	if (scoreType == "Phred64") offSet <- 64

	# turn the vector of cryptic ASCII string of scores into a matrix (and vector) of integers...	

	# to gaurantee a matrix, these must be equal length strings...
	alllens <- base::nchar( txt)
	maxlen <- max( alllens)
	if ( any( alllens < maxlen)) {
		# slower...
		mChars <- t( sapply( txt, function( x) { v <- strsplit( x, split="",fixed=T)[[1]]; n <- base::nchar(x); 
							if (n < maxlen) v[ (n+1):maxlen] <- " "; return(v) },
						 USE.NAMES=FALSE))
	} else {
		mChars <- t( sapply( txt, function( x) { strsplit( x, split="",fixed=T)[[1]]}, USE.NAMES=FALSE))
	}

	vChars <- as.vector( mChars)
	fac <- factor( vChars)
	facPtrs <- tapply( 1:length(vChars), INDEX=fac, FUN=NULL)

	# get the integer value for each unique score character
	oneStringOfAllScoreChars <- base::paste( levels(fac), collapse="") 
	intScores <- utf8ToInt( oneStringOfAllScoreChars) - offSet
	#cat( " PhredInt Speedup=", formatC( (length(vChars)/length(intScores)), digits=2, format="f"))

	vOut <- intScores[ facPtrs]
	mOut <- matrix( vOut, nrow=nrow(mChars), ncol=ncol(mChars))

	return( mOut)
}



`solexaScore2phredScore` <- function( x) {

	# map from Solexa scores (which can be negative) to Phred scores (always > 0)
	x <- as.numeric( x)
	out <- 10 * log( 1 + 10^( x/10.0)) / log(10)
	return( out)
}


`solexaCharScoresToPhredScores` <- function( x) {

	myScores <- x
	fac <- factor( myScores)
	uniqueScores <- levels( fac)
	facPtrs <- tapply( 1:length(myScores), fac, FUN=NULL)
	cat( "  Solexa Conv Speedup=", formatC( (length(myScores)/length(uniqueScores)), digits=2, format="f"))
	uniqueSolexaStrs <- sapply( uniqueScores, FUN=solexaCharScoreToSolexaIntScore)
	myScores <- uniqueSolexaStrs[ facPtrs]

	fac <- factor( myScores)
	uniqueScores <- levels( fac)
	facPtrs <- tapply( 1:length(myScores), fac, FUN=NULL)
	cat( "  Phred Conv Speedup=", formatC( (length(myScores)/length(uniqueScores)), digits=2, format="f"))
	uniquePhredStrs <- sapply( uniqueScores, FUN=solexaToPhred)
	myScores <- uniquePhredStrs[ facPtrs]

	return( myScores)
}

