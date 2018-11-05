\name{enrichment}
\alias{enrichment}
\alias{enrichment.Nway}
\title{
Enrichment Compared To Chance
}
\description{
Calculate the enrichment and P-value using the hypergeometric distribution.
}
\usage{
enrichment(nMatch, nYourSet, nTotal, nTargetSubset)

enrichment.Nway(nSets = 2, nMatch, nDrawn, nTotal, 
		nSimulations = 1000000)

simulate.enrichment.Nway(nSets = 2, nMatch, nDrawn, nTotal, 
		nSimulations = 1000000)
}
\arguments{
  \item{nMatch}{
  the number of genes in common from 'yourSet' and the 'targetSet'
}
  \item{nYourSet}{
  the number of genes in your selected subset
}
  \item{nTotal}{
  the total number of gene
}
  \item{nTargetSubset}{
  the number of genes in the target subset
}
  \item{nSets}{
  for extending to a higher number of samples
}
  \item{nSimulations}{
  number of simulations to run
}
}
\details{
Based on the hypergeometric function.  See \code{\link{dhyper}}.  Given a set 
of 'nTotal' genes, that contains a particular subset of 'nTargetSubset' genes 
of interest.  This represents the theoretical 'urn of balls' with this subset 
representing the white balls.  By some selection criteria, you draw out a 
subset of 'nYourSet' genes from the urn, and observe that 'nMatch' of your set are from the 
target subset (i.e. are white).  This function returns the expected number of 
white balls drawn, and the P-value likelihood of drawing 'nMatch' by chance alone.

For higher numbers of sets, there is a simulation-based method (that is quite
slow).  Imagine 4 independent differential expression results, and you wish to
know how likely it is to see 17 genes in common from the top 200 DE genes in a
genome with 5000 genes.  Use: 
\code{enrichment.Nway( nSets=4, nMatch=17, nDrawn=200, nTotal=5000)}

\code{simulate.enrichment.Nway} is a convenience wrapper function that
runs the simulation and graphically fits the
results to a normal distribution, to help estimate the probability of having a
matching overlap of that size.
}
\value{
A list that restates the inputs, with additional terms:
  \item{nExpected}{The expected number of matches, if pure chance was the only influence}
  \item{P_atLeast_N}{Probability of pulling at least 'nMatch' of the particular subset by chance}
  \item{P_atMost_N}{Probability of pulling no more than 'nMatch' of the particular subset by chance}

Additionally for \code{enrichment.Nway}:
  \item{distribution }{The results of the simulation, showing how often each possible 
  outcome (number of matches) was observed}

Additionally for \code{simulate.enrichment.Nway}:
  a plot of the distribution, with a best fit normal curve that highlights the number of
  matches.
}
\seealso{
See \code{\link{dhyper}}
}
