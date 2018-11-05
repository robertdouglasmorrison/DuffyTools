# stringTools.R  -- tools for strings

`clipLongString` <- function(str, max.length=256, pct.front=0.5) {

	nc <- nchar( str)
	needClip <- which( nc > max.length)
	if ( ! length(needClip)) return(str)

	ncWantFront <- floor( max.length * pct.front) - 3
	ncWantRear <- max.length - ncWantFront - 3
	fronts <- substr( str[needClip], 1, ncWantFront)
	if (ncWantRear > 0) {
		backs <- substr( str[needClip], (nc-ncWantRear+1)[needClip], nc[needClip])
	} else {
		backs <- rep.int( "", length(needClip))
	}
	newStr <- paste( fronts, backs, sep="...")
	out <- str
	out[ needClip] <- newStr
	out
}
