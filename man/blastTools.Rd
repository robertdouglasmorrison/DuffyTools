\name{readBlastOutput}
\alias{readBlastOutput}
\alias{callBlast}
\alias{callBlastp}
\alias{callBlastn}
\title{
Call BLAST and/or read results from a BLAST query
}
\description{
Wrapper functions to implement a call to BLAST, and turn the results of a BLAST query into an R table.
}
\usage{
callBlast( fastafile, outfile = "blastOut.txt", program = "blastn", db = "nt",
	path = Sys.getenv("BLASTINDEX_PATH"), task = "", wordsize = 8, evalue = 1, 
	threads = 4, outfmt = 6, filter = "no", maxhits = 5, blastArgs = "")
callBlastn( fastafile, outfile = "blastOut.txt", program = "blastn", db = "NT/nt",
	path = Sys.getenv("BLASTINDEX_PATH"), task = "", wordsize = 8, evalue = 1, 
	threads = 4, outfmt = 6, filter = "no", maxhits = 5, blastArgs = "")
callBlastp( fastafile, outfile = "blastOut.txt", program = "blastp", db = "NR/nr",
	path = Sys.getenv("BLASTINDEX_PATH"), task = "", wordsize = 3, evalue = 1, 
	threads = 4, outfmt = 6, filter = "no", maxhits = 5, blastArgs = "")

readBlastOutput( infile, outfmt = 6, verbose = TRUE, nKeep = NULL)
}
\arguments{
  \item{fastafile}{
Character string, full pathname to an existing FASTA file of query sequences to submit to BLAST.
}
  \item{outfile}{
Character string, full pathname for the resulting output from BLAST.
}
  \item{program}{
Character string, full pathname for the executable version of BLAST.  Passed to \code{Sys.which} for 
lookup via PATH environment variable.
}
  \item{db}{
Character string, name of the BLAST database to search against.
}
  \item{path}{
Character string, full pathname to the folder containing the BLAST database files.
}
  \item{outfmt}{
Blast output mode.  Currently, only \code{outfmt = 6} is supported.
}
  \item{infile}{
Character string, full pathname to an existing file created by BLAST using output mode '-outfmt=6'.
}
  \item{nKeep}{
When NULL, keep all BLAST hits for each query sequence, else an integer specifying the maximum number of 
BLAST matches to return for each query.  Only the highest scoring \code{nKeep} hits will be returned for 
each query string.
}
}
\details{
The BLAST command line details and result format assume a BLAST version of at least 2.2.25
}
\value{
For the calls to BLAST, no return value, just a file of results written to disk.

For \code{reaedBlastOutput}, a data frame with columns:  
	PROBE_ID, SEQ_ID, PCT_MATCH, LEN_MATCH, MIS_MATCH, GAP, 
	P_FIRST, P_LAST, S_BEG, S_END, E_VALUE, SCORE, STRAND
that contains the BLAST results from the file.
}
