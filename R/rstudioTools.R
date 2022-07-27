#  rstudioTools.R -- special routines to interact seemlessly with Rstudio


`is.Rstudio` <- function() {

	# determine if we are running inside an Rstudio session
	# there should be environment variables as a clue
	return( Sys.getenv("RSTUDIO") == "1")
}

