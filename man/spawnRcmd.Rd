\name{spawnRcmd}
\alias{spawnRcmd}
\title{
Launch R Command as a Child Process
}
\description{
launch one (or more) R commands as a separate child process.
}
\usage{
spawnRcmd(rCommands, Rpackage = "DuffyTools", logFile = "spawnRcmd.log.txt")
}
\arguments{
  \item{rCommands}{
one character string that is a series of valid R comands to execute.  Care
must be taken to have embedded character strings awways come through intact.
}
  \item{Rpackage}{
one R package to preload, via \code{\link{library}}
}
  \item{logFile}{
created text file that captures stdout and stderr of the child process.
}
}
\details{
This is a rather crude method to spawn off separate subtasks to speed up
overall processing.
}
\note{
Interrupting the main R session via ctrl-C will send a kill signal to 
all spawned child processes.
}
