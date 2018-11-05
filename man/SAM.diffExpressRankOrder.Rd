\name{SAM.diffExpressRankOrder}
\alias{SAM.diffExpressRankOrder}
\title{
Order Differentially Expressed Genes
}
\description{
Orders the results of a differential expression calculation, as a 
combination of the reported fold change, P-value, and SAM
distance for each gene.
}
\usage{
SAM.diffExpressRankOrder(folds, pvalues, dists, 
		wt.folds = 1, wt.pvalues = 1, wt.dists = 1, 
		two.sided.P.values = TRUE)
}
\arguments{
  \item{folds}{
numeric vector of fold change values.  Must be expresses as log2 
ratio values (no expression change equals 0)
}
  \item{pvalues}{
numeric vector of same length as \code{folds} giving the P-value 
associated with each fold change value
}
  \item{dists}{
numeric vector of distances 'D' as measured by SAM
}
  \item{wt.folds}{
scalar weight for the contribution of fold change rankings to the result
}
  \item{wt.pvalues}{
scalar weight for the contribution of P-value rankings to the result
}
  \item{wt.dists}{
scalar weight for the contribution of SAM distance to the result
}
  \item{two.sided.P.values}{
logical specifying whether the P-values came from a 2-sided or 1-sided 
test
}
}
\details{
Fold change is a way of describing the difference between two expression 
values as a ratio \code{ x1/x2 }.  For very low magnitudes, the fold 
change can be artifically high.  For very high magnitudes, the P-value can
be artifically small (especially in NextGen sequencing datasets).  Further
algorithms, such as SAM, have their own metrics of the statistical 
significance of differential expression.  This 
function calcualates the rank order by each criteria separately, and 
then uses a weighted average of those to produce the resulting gene rank order.
}
\value{
A vector of integer indices the same length as \code{folds}, giving the 
permutation order of the input from highest differentially expressed 
(most up-regulated) to lowest (most down-regulated).
}
\note{
In a usual 2-sided test, both the strongly up-regulated and strongly
down-regulated genes will have good P-values.  The function assumes this
about the given P-values.  If a 1-sided test that gives poor P-values for
the down-regulated genes was used, set \code{two.sided.P.value=FALSE}
}

\seealso{
  \code{\link{diffExpressRankOrder}}
}
