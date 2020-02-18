# fileSemaphore.R


# adding some 'semaphore' type construct to avoid race conditions...


`file.lock` <- function( filename, ID="", sleeptime=3, retries=60) {

	return(TRUE)
	lockFilename <- base::paste( filename, "LOCK", sep=".")

	sleeptime <- round( sleeptime);
	retries <- round( retries);
	timeout <- round( sleeptime * retries);

	lockArgs <- base::paste( "-", sleeptime, " -r ", retries, " -l ", timeout, " ", sep="")
	cmd <- base::paste( "lockfile ", lockArgs, lockFilename)

	ans <- catch.system( command=cmd, intern=FALSE)

	return( ans == 0)
}


`file.unlock` <- function( filename, ID) {

	return(TRUE)
	lockFilename <- base::paste( filename, "LOCK", sep=".")

	catch.system( command=base::paste( "rm -f ", lockFilename))

	return( TRUE)
}
