# chromatogramPairTools.R -- pieces to try to cross-validate and clean a pair of Fwd / Rev chromatograms
#				with the aim of finding & fixing errors that make the 2 inconsistent


`prepChromatogramPair` <- function( fwdFile, revFile) {

	# read up 2 AB1 files from disk... into our standard format of trace, peaks, & sequence
	# and write them as R objects and .CSV tables for editting
	require( "sangerseqR")

	# read them in
	if ( ! file.exists( fwdFile)) {
		cat( "\nForward Chromatogram file not found:  ", fwdFile)
		return( NULL)
	}
	fwdObj <- loadChromatogram( fwdFile)
	if ( ! file.exists( revFile)) {
		cat( "\nReverse Chromatogram file not found:  ", revFile)
		return( NULL)
	}
	revObj <- loadChromatogram( revFile)

	# prep file name roots for other files we will create
	fwdFileRoot <- sub( "\\.ab1$", "", fwdFile)
	revFileRoot <- sub( "\\.ab1$", "", revFile)

	# for the forward one, it is ready as is
	fwdTbl <- chromatogramToTable( fwdObj)
	fwdTblFile <- paste( fwdFileRoot, "ChromatogramDetails.csv", sep=".")
	write.table( fwdTbl, fwdTblFile, sep=",", quote=T, row.names=F, na="")
	fwdObjFile <- paste( fwdFileRoot, "ChromatogramDetails.rda", sep=".")
	save( fwdObj, file=fwdObjFile)
	fwdDNA <- fwdObj$DNA_Calls[1]

	# for the Reverse chromatogram, we need to make the RevComp for most activities we will perform
	revObj <- revCompChromatogram( revObj)

	# now ready to save the reverse one
	revTbl <- chromatogramToTable( revObj)
	revTblFile <- paste( revFileRoot, "ChromatogramDetails.csv", sep=".")
	write.table( revTbl, revTblFile, sep=",", quote=T, row.names=F, na="")
	revObjFile <- paste( revFileRoot, "ChromatogramDetails.rda", sep=".")
	save( revObj, file=revObjFile)
	revDNA <- revObj$DNA_Calls[1]

	# ideally we have some good overlap of these two DNA chunks
	#cat( "\nFwd: ", nchar(fwdDNA), "|", fwdDNA)
	#cat( "\nRev: ", nchar(revDNA), "|", revDNA)

	pa <- pairwiseAlignment( fwdDNA, revDNA, type="local")
	fwdOut <- as.character( alignedPattern( pa))
	revOut <- as.character( alignedSubject( pa))
	fwdStart <- start( pattern( pa))
	revStart <- start( subject( pa))

	#PA <<- pa

	out <- data.frame( "Chromat"=c( "Fwd", "Rev"), "Length"=nchar( c(fwdDNA,revDNA)), 
			"Starts"=c(fwdStart,revStart), "DNA_Alignment"=c(fwdOut, revOut), stringsAsFactors=F)
	out
}


`inspectChromatogramPair` <- function( fwdFile, revFile, fwdStart, revStart, fwdShow=20, revShow=fwdShow,
					fwdAA=1, revAA=fwdAA, substring=NULL, showTraceRowNumbers=T) {

	# visually inspect these 2 chromatograms at one region of overlapping interest
	# read them in
	fwdFileRoot <- sub( "\\.ab1$", "", fwdFile)
	revFileRoot <- sub( "\\.ab1$", "", revFile)
	fwdTblFile <- paste( fwdFileRoot, "ChromatogramDetails.csv", sep=".")
	fwdObjFile <- paste( fwdFileRoot, "ChromatogramDetails.rda", sep=".")
	revTblFile <- paste( revFileRoot, "ChromatogramDetails.csv", sep=".")
	revObjFile <- paste( revFileRoot, "ChromatogramDetails.rda", sep=".")
	if ( ! file.exists( fwdObjFile)) {
		cat( "\nForward Chromatogram R object file not found:  ", fwdObjFile)
		return( NULL)
	}
	load( fwdObjFile)
	if ( ! file.exists( revObjFile)) {
		cat( "\nReverse Chromatogram R object file not found:  ", revFObjile)
		return( NULL)
	}
	load( revObjFile)

	# we will draw
	checkX11()
	par( mfcol=c(2,1))

	# get that subset from both
	if ( is.null( substring)) {
		fwdSubObj <- subsetChromatogram( fwdObj, range=fwdStart:(fwdStart+fwdShow))
		revSubObj <- subsetChromatogram( revObj, range=revStart:(revStart+revShow))
	} else {
		fwdSubObj <- subsetChromatogram( fwdObj, substring=substring)
		revSubObj <- subsetChromatogram( revObj, substring=substring)
	}

	# draw them
	plotChromatogram( fwdSubObj, label="Forward", showAA=fwdAA, showTraceRowNumbers=showTraceRowNumbers)
	plotChromatogram( revSubObj, label="Reverse", showAA=revAA, showTraceRowNumbers=showTraceRowNumbers)

}
