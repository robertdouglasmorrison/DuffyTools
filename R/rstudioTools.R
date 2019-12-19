#  rstudioTools.R -- special routines to interact seemlessly with Rstudio


`is.Rstudio` <- function() {

	# determine if we are running inside an Rstudio session
	# there should be environment variables as a clue
	envSet <- Sys.getenv()
	nRstudio <- length( grep( "^RSTUDIO", names(envSet)))
	return( nRstudio > 0)
}

