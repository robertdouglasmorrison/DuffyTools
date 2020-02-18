# systemTools.R  -- tools for wrapping system calls

`catch.system` <- function( command, intern=FALSE, ignore.stdout=FALSE,
			ignore.stderr=FALSE, wait=TRUE, input=NULL,
			timeout=0, ...) {

	# just call the normal 'system' command, but watch the status of the 
	# completed command, to correctly abort the workflow
	cmdStatus <- system( command=command, intern=intern, ignore.stdout=ignore.stdout,
			ignore.stderr=ignore.stderr, wait=wait, input=input, 
			timeout=timeout, ...)
	
	# how to know the command status depends on 'intern'
	if ( intern == TRUE) {
		ans <- as.character( cmdStatus)
		status <- if ( "status" %in% names(attributes(cmdStatus))) attributes(cmdStatus)$status else 0
	} else {
		status <- cmdStatus
		ans <- status
	}

	if ( status != 0) stop( paste( "Aborting.  'system()' returned non-zero status: ", status))

	if ( intern == TRUE) {
		return( ans)
	} else {
	}
		return( invisible( ans))
}
