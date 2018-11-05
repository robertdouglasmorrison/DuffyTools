# test.DuffyTools.R

`test.DuffyTools` <- function() {

	# primitives, etc
	test.as.percent()
	cat( "\nPassed utility tests.")

	# simple, 'low-level', routines
	test.readFastqFile()
	test.writeFastqFile()
	cat( "\nPassed FASTQ file tests.")

	# higher level processing tests
	test.mapsetTools()
	test.speciesTools()
	cat( "\nPassed species mapSet tests.")
	
	test.codonTools()
	cat( "\nPassed DNA, AA codon tests.")

	# ... none yet ...

	cat( "\nPassed all DuffyTools self tests.\n")
}


# little bits to make any test data we may need

`remove.testFile` <- function( filename) { file.remove( filename) }


`build.testFastqFile` <- function() {

	filename <- "DuffyTools.test.fastq"
	got <- data( "testFastq.txt", package="DuffyTools")
	if ( got != "testFastq.txt") warning( "DuffyTools:  loading test data failed.")
	writeLines( testFastq.txt, con=filename)
	rm( testFastq.txt, envir=.GlobalEnv)
	return( filename)
}

