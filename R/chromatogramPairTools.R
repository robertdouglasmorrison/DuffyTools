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
	writeCuratedChromatogram(fwdObj)
	fwdDNA <- fwdObj$DNA_Calls[1]

	if ( ! file.exists( revFile)) {
		cat( "\nReverse Chromatogram file not found:  ", revFile)
		return( NULL)
	}
	revObj <- loadChromatogram( revFile)

	# for the Reverse chromatogram, we need to make the RevComp for most activities we will perform
	revObj <- revCompChromatogram( revObj)
	writeCuratedChromatogram(revObj)
	revDNA <- revObj$DNA_Calls[1]

	# ideally we have some good overlap of these two DNA chunks
	pa <- pairwiseAlignment( fwdDNA, revDNA, type="local")
	fwdOut <- as.character( alignedPattern( pa))
	revOut <- as.character( alignedSubject( pa))
	fwdStart <- start( pattern( pa))
	revStart <- start( subject( pa))

	out <- data.frame( "Chromat"=c( "Fwd", "Rev"), "Length"=nchar( c(fwdDNA,revDNA)), 
			"Starts"=c(fwdStart,revStart), "DNA_Alignment"=c(fwdOut, revOut), stringsAsFactors=F)
	out
}


`inspectChromatogramPair` <- function( fwdFile, revFile, fwdStart, revStart, fwdShow=20, revShow=fwdShow,
					fwdAA=1, revAA=fwdAA, substring=NULL, showTraceRowNumbers=T, showConfidence=T) {

	# visually inspect these 2 chromatograms at one region of overlapping interest
	# read them in
	fwdObj <- loadCuratedChromatogram( fwdFile)
	revObj <- loadCuratedChromatogram( revFile)

	# we will draw
	checkX11()
	par( mfcol=c(2,1))
	on.exit( par( mfcol=c(1,1)))

	# get that subset from both
	if ( is.null( substring)) {
		fwdSubObj <- subsetChromatogramByRange( fwdObj, range=fwdStart:(fwdStart+fwdShow))
		revSubObj <- subsetChromatogramByRange( revObj, range=revStart:(revStart+revShow))
	} else {
		fwdSubObj <- subsetChromatogramBySequence( fwdObj, seq=substring)
		revSubObj <- subsetChromatogramBySequence( revObj, seq=substring)
	}

	# draw them
	plotChromatogram( fwdSubObj, label="Forward", showAA=fwdAA, showTraceRowNumbers=showTraceRowNumbers,
					showConfidence=showConfidence)
	plotChromatogram( revSubObj, label="Reverse", showAA=revAA, showTraceRowNumbers=showTraceRowNumbers,
					showConfidence=showConfidence)
}
