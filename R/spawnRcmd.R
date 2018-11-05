`spawnRcmd` <- function( rCommands, Rpackage="DuffyTools", logFile="spawnRcmd.log.txt",
			verbose=TRUE) {

	# submit an R command as a background process
	Rprogram <- Sys.getenv( "R_PROGRAM")
	how <- "Sys.getenv('R_PROGRAM')"
	if ( nchar( Rprogram) < 1) {
		Rprogram <- Sys.which( "R")
		how <- "Sys.which('R')"
	}
	if ( nchar( Rprogram) < 1) {
		Rprogram <- "R"
		how <- "default"
	}

	# turn the whole set into one compound expression
	firstCmd <- paste( "library(", Rpackage, ")")
	lastCmd <- quote( print( paste( "Spawned_R_session_finished:  ", date())))
	rCommands <- c( firstCmd, rCommands, lastCmd)
	RcommandLines <- paste( rCommands, collapse=";  ")
	RcommandLines <- paste( "{", RcommandLines, "}")

	# squeeze out spaces
	RcommandLines <- gsub( " ", "", RcommandLines)

	# build the UNIX command line that starts up R on this set of commands
	RcommandLine <- shQuote( RcommandLines)
	commandLine <- paste( Rprogram, " --quiet  --no-restore  -e ", RcommandLine,
				"  >> ", logFile, " 2>&1")

	# say 'what will happen' to the current process std out.
	cat( "\n\nSpawning R Command session: ")
	cat( "\n  Found R by: ", how, "  -> ", Rprogram)
	cat( "\n  Log File:   ", logFile)
	cat( "\n  Command:    ", commandLine, "\n")

	# create the file that will get appended to...
	file.delete( logFile)
	if (verbose) {
		sink( logFile)
		# say it again to the new process's std out...
		cat( "\nSpawned R Command Session \n  Time Start:      ", date())
		cat( "\n  Command:   ", commandLine, "\n\n\n")
		sink()
	}

	# run it!
	system( command=commandLine, wait=FALSE, ignore.stderr=FALSE)

	return()
}
