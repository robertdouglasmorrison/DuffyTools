# baseDepthTables.R   - tools to turn a set of begin:end base positions into a table of depths,
#			and other such operations...


`baseDepthVector` <- function( start, stop, weight=1) {

	# given a set of begin:end base locations, return a 1D table where the names are the
	# locations and and the depth at each base is the value
	if ( length( start) < 1) return( EMPTY_BASE_DEPTH_VECTOR)

	# trap any NAs...
	if ( any( drops <- which( is.na(start)))) {
		start <- start[ -drops]
		stop <- stop[ -drops]
		if ( length( start) < 1) return( EMPTY_BASE_DEPTH_VECTOR)
	}
	if ( any( drops <- which( is.na(stop)))) {
		start <- start[ -drops]
		stop <- stop[ -drops]
		if ( length( start) < 1) return( EMPTY_BASE_DEPTH_VECTOR)
	}
	# trap any negative size ranges.  Tey should never occur, but if present throw off the size estimate
	if ( any( drops <- which( start > stop))) {
		start <- start[ -drops]
		stop <- stop[ -drops]
		if ( length( start) < 1) return( EMPTY_BASE_DEPTH_VECTOR)
	}

	bigN <- sum( stop - start + 1, na.rm=T)
	if ( bigN < 1) return( EMPTY_BASE_DEPTH_VECTOR)

	allWts <- vector( mode="numeric", length=bigN)
	allLocs <- vector( mode="integer", length=bigN)
	nout <- 0

	base::mapply( FUN=function( start, stop, weight) {
				locs <- start : stop
				nnow <- length(locs)
				who <- (nout+1) : (nout+nnow)
				allLocs[ who] <<- locs
				allWts[ who] <<- rep.int( weight[1], nnow)
				nout <<- nout + nnow
				return( NULL)
			}, start=start, stop=stop, weight=weight, MoreArgs=NULL, SIMPLIFY=FALSE, USE.NAMES=FALSE)

	length(allWts) <- length(allLocs) <- nout
	fac <- factorBaseLocations( allLocs)
	totalWts <- base::sapply( base::split.default( allWts, fac), FUN=sum)
	return( totalWts)
}


`baseDepthTable` <- function( start, stop, weight=1) {

	# given a set of begin:end base locations, return a 2D table of {begin, end, depth}
	# like wanted for wiggle tracks

	bases1D <- baseDepthVector( start, stop, weight)

	bases2D <- baseDepthVectorToTable( bases1D)

	return( bases2D)
}


`emptyBaseDepthVector` <- function() {

	return( vector( mode="numeric", length=0))
}


`emptyBaseDepthTable` <- function() {

	nodata <- vector( mode="numeric", length=0)
	noloc <- vector( mode="numeric", length=0)
	return( data.frame( "START"=noloc, "STOP"=noloc, "DEPTH"=nodata))
}



`baseDepthTableSubset` <- function( baseTbl, start, stop) {

	# given a base depth table, extract the subset of rows that contain the start:stop region
	# like a subset of wiggle tracks
	if ( is.null(baseTbl)) return( EMPTY_BASE_DEPTH_TABLE)
	N <- nrow( baseTbl)
	if ( N < 1) return( EMPTY_BASE_DEPTH_TABLE)

	# try to be as fast as possible...

	# get the pointers into these rows that contain this start:stop region
	locstart <- fastFindInterval( x=start[1], vec=baseTbl$STOP) + 1
	locstop <- fastFindInterval( x=stop[1], vec=baseTbl$START)

	# OK, we have our range of rows... is it non-empty?
	if ( locstop - locstart < 0) return( EMPTY_BASE_DEPTH_TABLE)

	outX <- locstart : locstop 
	out <- data.frame( "START"=baseTbl$START[ outX], 
			"STOP"=baseTbl$STOP[ outX], 
			"DEPTH"=baseTbl$DEPTH[ outX])

	# tweak the end points to exactly match this interval
	if ( out$START[1] < start) out$START[1] <- start
	N <- nrow(out)
	if ( out$STOP[N] > stop) out$STOP[N] <- stop
	return( out)
}


`sum.baseDepthTable` <- function( baseTbl) {

	if ( is.null( baseTbl)) return( 0)
	segLengths <- baseTbl$STOP - baseTbl$START + 1
	counts <- segLengths * baseTbl$DEPTH
	return( sum( counts))
}


`sum.baseDepthTableSubset` <- function( baseTbl, start, stop) {

	# given a base depth table, extract the subset of rows that contain the start:stop region
	# like a subset of wiggle tracks
	if ( is.null( baseTbl)) return( 0)
	N <- nrow( baseTbl)
	if ( N < 1) return( 0)

	# try to be as fast as possible...

	# get the pointers into these rows that contain this start:stop region
	locstart <- fastFindInterval( x=start[1], vec=baseTbl$STOP) + 1
	locstop <- fastFindInterval( x=stop[1], vec=baseTbl$START)

	# OK, we have our range of rows... is it non-empty?
	if ( locstop - locstart < 0) return( 0)

	outX <- locstart : locstop 
	myStarts <- baseTbl$START[ outX]
	myStops <- baseTbl$STOP[ outX]
	myDepths <- baseTbl$DEPTH[ outX]

	# tweak the end points to exactly match this interval
	if ( myStarts[1] < start) myStarts[1] <- start
	N <- length( outX)
	if ( myStops[N] > stop) myStops[N] <- stop

	# count up those bases
	segLengths <- myStops - myStarts + 1
	counts <- segLengths * myDepths
	return( sum( counts))
}


`baseDepthVectorToTable` <- function( baseVec) {

	# given a 1D table, where the names are the positions, and the values are the depths
	# we want to make 1 row for each consectutive set of bases at a given depth...

	N <- length( baseVec)
	if ( N < 2) return( EMPTY_BASE_DEPTH_TABLE)

	baseLocs <- as.numeric( names( baseVec))
	baseDepths <- as.numeric( baseVec)
	isConsectutive <- ( diff( baseLocs) == 1)
	isSameDepth <- ( baseDepths[1:(N-1)] == baseDepths[2:N] )
	stairSteps <- which( ! ( isConsectutive & isSameDepth))

	# we always do a stair step at the very first and very last
	allFroms <- c( 1, stairSteps+1)
	allTos <- c( stairSteps, N)
	outStart <- baseLocs[ allFroms]
	outStops <- baseLocs[ allTos]
	outDepth <- baseDepths[ allFroms]

	out <- data.frame( "START"=as.numeric(outStart), "STOP"=as.numeric(outStops), "DEPTH"=outDepth)
	rownames(out) <- 1:nrow(out)

	return( out)
}


`baseDepthTableToVector` <- function( baseTbl) {

	# given a 2D table, where the rows are {start, stop, depth},
	# we want to make a 1D vector where the names are base locations and the values are the depths...

	if ( is.null( baseTbl)) return( EMPTY_BASE_DEPTH_VECTOR)
	if ( nrow( baseTbl) < 1) return( EMPTY_BASE_DEPTH_VECTOR)

	inStart <- baseTbl$START

	# make sure the result is in increasing order...
	ord <- base::order( inStart)
	if ( ! all( ord == 1:nrow(baseTbl))) {
		baseTbl <- baseTbl[ ord, ]
		inStart <- baseTbl$START
	}
	inStop <- baseTbl$STOP
	inDepth <- baseTbl$DEPTH

	bigN <- sum( inStop - inStart + 1)
	outDepth <- vector( mode="numeric", length=bigN)
	outName <- vector( mode="numeric", length=bigN)
	nout <- 0

	base::mapply( FUN=function( start, stop, depth) {
				x <- start[1] : stop[1]
				nnow <- length(x)
				who <- (nout+1) : (nout+nnow)
				outDepth[ who] <<- rep.int( depth[1], nnow)
				outName[ who] <<- x
				nout <<- nout + nnow
				return(NULL)
			}, start=inStart, stop=inStop, depth=inDepth, MoreArgs=NULL, SIMPLIFY=FALSE, USE.NAMES=FALSE)

	length(outDepth) <- length(outName) <- nout

	names( outDepth) <- outName
	return( outDepth)
}


`merge.baseDepthTables` <- function( baseTbl1, baseTbl2) {

	# given 2 2D tables, add their depths
	if ( nrow( baseTbl1) < 1) return( baseTbl2)
	if ( nrow( baseTbl2) < 1) return( baseTbl1)

	baseVec1 <- baseDepthTableToVector( baseTbl1)
	baseVec2 <- baseDepthTableToVector( baseTbl2)

	ans <- mergeIntegerTables( baseVec1, baseVec2)

	return( baseDepthVectorToTable( ans))
}


`factorBaseLocations` <- function (x) {
 
    y <- unique.default(x)
    y <- sort.int(y)
    f <- base::match(x, y)
    levels(f) <- y
    class(f) <- "factor"
    f
}


`baseDepthTableSqueeze` <- function( baseTbl, squeezeFactor, squeeze.FUN=max) {

	# given a base depth table, reduce the number of rows by 'Squeezing' multiple bases
	# into a single row -- typically to speed up rendering in gene plots

	#  squeezeFactor is how many real base should resolve down to one row when done
	# ignore if minimal improvement
	# FUN is the function to apply to those depths at each new point

	N <- nrow( baseTbl)
	if ( N < 1) return( EMPTY_BASE_DEPTH_TABLE)

	if ( squeezeFactor < 3) return( baseTbl)

	# try to be as fast as possible...
	midBase <- (baseTbl$START + baseTbl$STOP) / 2
	newBase  <- as.integer( round( midBase / squeezeFactor) * squeezeFactor)
	baseFac <- factor( newBase)
	newStart <- newStop <- newDepth <- vector()
	nnow <- 0

	# build a new table, using all rows that will map to one new row
	tapply( 1:N, baseFac, function(x) {
		nnow <<- nnow + 1
		newStart[nnow] <<- baseTbl$START[x[1]]
		newStop[nnow] <<- baseTbl$STOP[x[length(x)]]
		newDepth[nnow] <<- squeeze.FUN( baseTbl$DEPTH[x])
		})

	out <- data.frame( "START"=newStart, "STOP"=newStop, "DEPTH"=newDepth)

	return( out)
}

