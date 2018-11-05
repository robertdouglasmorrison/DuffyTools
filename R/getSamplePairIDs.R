# getSamplePairIDs -- convert a sampleID into possibly a pair of IDs


`getSamplePairIDs` <- function( sampleID, annotationFile="Annotation.txt") {

	trueID <- originalSamplePairID( sampleID, annotationFile)

	paired <- getAnnotationTrue( annotationFile, key=trueID, columnArg="PairedEnd", verbose=FALSE)
	if (paired) {
		return( paste( trueID, c("1","2"), sep="_"))
	} else {
		return( trueID)
	}
}


`originalSamplePairID` <- function( sampleID, annotationFile="Annotation.txt") {

	# look carefully to find the sampleID
	annT <- readAnnotationTable( annotationFile)

	# find this sampleID in the set of known samples
	who <- match( sampleID[1], annT$SampleID, nomatch=0)
	if ( who == 0) {
		childPair <- TRUE
		tryID <- sub( "_[12]$", "", sampleID[1])
		who <- match( tryID, annT$SampleID, nomatch=0)
	} else {
		childPair <- FALSE
	}
	if ( who == 0) {
		stop( paste( "\nSampleID not found in Annotation Table: ", sampleID[1], "\n"))
	}
	trueID <- annT$SampleID[who]
	paired <- getAnnotationTrue( annT, key=trueID, columnArg="PairedEnd", verbose=FALSE)

	if (paired && childPair) {
		return( trueID)
	} else {
		return( trueID)
	}
}
