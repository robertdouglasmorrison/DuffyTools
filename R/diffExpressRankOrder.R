# diffExpressRankOrder.R

# find the best ordering for a set of genes, based on a combination of
# fold change and P-value

`diffExpressRankOrder` <- function( folds, pvalues, wt.folds=1, wt.pvalues=1, notDE.value=0, two.sided.P.values=TRUE) {

	N <- length( folds)
	if ( length( pvalues) != N) stop( "FoldChange and Pvalue vectors must be same length")

	# we sort them in 3 groups,  the ups, the downs, the rest...
	ranks <- 1:N
	whoUp <- which( folds > notDE.value)
	whoZero <- which( folds == notDE.value)
	whoDown <- which( folds < notDE.value)

	# now sort just the UPs based on Pvalue and fold..
	Nup <- length( whoUp)
	rankUp_P <- rankUp_F <- vector( length=Nup)
	ordUp_P <- base::order( pvalues[ whoUp], (-folds[whoUp]), decreasing=FALSE)
	ordUp_F <- base::order( folds[ whoUp], (-pvalues[whoUp]), decreasing=TRUE)
	rankUp_P[ ordUp_P] <- 1:Nup
	rankUp_F[ ordUp_F] <- 1:Nup

	#best rank is the average rank from both criteria
	ordUp <- base::order( ((rankUp_P*wt.pvalues) + (rankUp_F*wt.folds))/(wt.pvalues+wt.folds))
	outUp <- whoUp[ ordUp]

	# the zero folds may not all have P=1.0, so sort on those
	ordZero <- base::order( pvalues[ whoZero])
	outZero <- whoZero[ ordZero]

	# now sort just the DOWNs based on Pvalue and fold..
	Ndown <- length( whoDown)
	rankDown_P <- rankDown_F <- vector( length=Ndown)

	# if its not a two-sided Pvalue test, the very negative folds have Pvalues near 1
	if ( ! two.sided.P.values) pvalues <- 1.0 - pvalues
	ordDown_P <- base::order( pvalues[ whoDown], folds[whoDown], decreasing=TRUE)
	ordDown_F <- base::order( folds[ whoDown], pvalues[whoDown], decreasing=TRUE)
	rankDown_P[ ordDown_P] <- 1:Ndown
	rankDown_F[ ordDown_F] <- 1:Ndown
	#best rank is the average rank from both criteria
	ordDown <- base::order( ((rankDown_P*wt.pvalues) + (rankDown_F*wt.folds))/(wt.pvalues+wt.folds))
	outDown <- whoDown[ ordDown]

	out <- c( outUp, outZero, outDown)

	# make sure that we still have every index exactly once
	if ( length( setdiff( out, ranks)) != 0) warning( "diffExpessRankOrder:  ordering by Fold and Rank problem...")
	return( out)
}



`diffExpressDistanceRankOrder` <- function( folds, pvalues, dists, wt.folds=1, wt.pvalues=1, 
					wt.dists=1, notDE.value=0, two.sided.P.values=TRUE) {

	N <- length( folds)
	if ( length( pvalues) != N) stop( "FoldChange and Pvalue vectors must be same length")
	if ( length( dists) != N) stop( "FoldChange and Distance vectors must be same length")

	# we sort them in 3 groups,  the ups, the downs, the rest...
	ranks <- 1:N
	whoUp <- which( folds > notDE.value)
	whoZero <- which( folds == notDE.value)
	whoDown <- which( folds < notDE.value)

	# now sort just the UPs based on Pvalue and fold and distance..
	Nup <- length( whoUp)
	rankUp_P <- rankUp_F <- rankUp_D <- vector( length=Nup)
	ordUp_P <- base::order( pvalues[ whoUp], (-folds[whoUp]), decreasing=FALSE)
	ordUp_F <- base::order( folds[ whoUp], (-pvalues[whoUp]), decreasing=TRUE)
	ordUp_D <- base::order( dists[ whoUp], (-pvalues[whoUp]), decreasing=TRUE)
	rankUp_P[ ordUp_P] <- 1:Nup
	rankUp_F[ ordUp_F] <- 1:Nup
	rankUp_D[ ordUp_D] <- 1:Nup

	#best rank is the average rank from all criteria
	ordUp <- base::order( ((rankUp_P*wt.pvalues) + (rankUp_F*wt.folds) + (rankUp_D*wt.dists))
				/ (wt.pvalues+wt.folds+wt.dists))
	outUp <- whoUp[ ordUp]

	# the zero folds may not all have P=1.0, so sort on those
	ordZero <- base::order( pvalues[ whoZero])
	outZero <- whoZero[ ordZero]

	# now sort just the DOWNs based on Pvalue and fold and distance..
	Ndown <- length( whoDown)
	rankDown_P <- rankDown_F <- rankDown_D <- vector( length=Ndown)

	# if its not a two-sided Pvalue test, the very negative folds have Pvalues near 1
	if ( ! two.sided.P.values) pvalues <- 1.0 - pvalues
	ordDown_P <- base::order( pvalues[ whoDown], folds[whoDown], decreasing=TRUE)
	ordDown_F <- base::order( folds[ whoDown], pvalues[whoDown], decreasing=TRUE)
	ordDown_D <- base::order( dists[ whoDown], pvalues[whoDown], decreasing=TRUE)
	rankDown_P[ ordDown_P] <- 1:Ndown
	rankDown_F[ ordDown_F] <- 1:Ndown
	rankDown_D[ ordDown_D] <- 1:Ndown
	#best rank is the average rank from all criteria
	ordDown <- base::order( ((rankDown_P*wt.pvalues) + (rankDown_F*wt.folds) + (rankDown_D+wt.dists))
				/ (wt.pvalues+wt.folds+wt.dists))
	outDown <- whoDown[ ordDown]

	out <- c( outUp, outZero, outDown)

	# make sure that we still have every index exactly once
	if ( length( setdiff( out, ranks)) != 0) warning( "diffExpessDistanceRankOrder:  ordering by Fold,Pvalue,Distance problem...")
	return( out)
}


`diffExpressRankRankOrder` <- function( folds, pvalues, ranks, wt.folds=1, wt.pvalues=1, 
					wt.ranks=1, notDE.value=0, two.sided.P.values=TRUE) {

	N <- length( folds)
	if ( length( pvalues) != N) stop( "FoldChange and Pvalue vectors must be same length")
	if ( length( ranks) != N) stop( "FoldChange and Ranks vectors must be same length")

	# we sort them in 3 groups,  the ups, the downs, the rest...
	myRanks <- 1:N
	whoUp <- which( folds > notDE.value)
	whoZero <- which( folds == notDE.value)
	whoDown <- which( folds < notDE.value)

	# now sort just the UPs based on Pvalue and fold and Rank..
	Nup <- length( whoUp)
	rankUp_P <- rankUp_F <- rankUp_R <- vector( length=Nup)
	ordUp_P <- base::order( pvalues[ whoUp], (-folds[whoUp]), decreasing=FALSE)
	ordUp_F <- base::order( folds[ whoUp], (-pvalues[whoUp]), decreasing=TRUE)
	ordUp_R <- base::order( ranks[ whoUp], (-folds[whoUp]), decreasing=FALSE)
	rankUp_P[ ordUp_P] <- 1:Nup
	rankUp_F[ ordUp_F] <- 1:Nup
	rankUp_R[ ordUp_R] <- 1:Nup

	#best rank is the average rank from all 3 criteria
	ordUp <- base::order( ((rankUp_P*wt.pvalues) + (rankUp_F*wt.folds) + 
				(rankUp_R*wt.ranks))/(wt.pvalues+wt.folds+wt.ranks))
	outUp <- whoUp[ ordUp]

	# the zero folds may not all have P=1.0, so sort on those
	ordZero <- base::order( pvalues[ whoZero], ranks[ whoZero])
	outZero <- whoZero[ ordZero]

	# now sort just the DOWNs based on Pvalue and fold and Rank..
	Ndown <- length( whoDown)
	rankDown_P <- rankDown_F <- rankDown_R <- vector( length=Ndown)

	# if its not a two-sided Pvalue test, the very negative folds have Pvalues near 1
	if ( ! two.sided.P.values) pvalues <- 1.0 - pvalues
	ordDown_P <- base::order( pvalues[ whoDown], folds[whoDown], decreasing=TRUE)
	ordDown_F <- base::order( folds[ whoDown], pvalues[whoDown], decreasing=TRUE)
	ordDown_R <- base::order( ranks[ whoDown], folds[whoDown], decreasing=FALSE)
	rankDown_P[ ordDown_P] <- 1:Ndown
	rankDown_F[ ordDown_F] <- 1:Ndown
	rankDown_R[ ordDown_R] <- 1:Ndown
	#best rank is the average rank from both criteria
	ordDown <- base::order( ((rankDown_P*wt.pvalues) + (rankDown_F*wt.folds) + (rankDown_R*wt.ranks))/(wt.pvalues+wt.folds+wt.ranks))
	outDown <- whoDown[ ordDown]

	out <- c( outUp, outZero, outDown)

	# make sure that we still have every index exactly once
	if ( length( setdiff( out, myRanks)) != 0) warning( "diffExpessRankOrder:  ordering by Fold and Rank problem...")
	return( out)
}


`diffExpressMetaResultOrder` <- function( folds, pvalues, ranks, wt.folds=1, wt.pvalues=1, 
					wt.ranks=1, notDE.value=0) {

	N <- length( folds)
	if ( length( pvalues) != N) stop( "FoldChange and Pvalue vectors must be same length")
	if ( length( ranks) != N) stop( "FoldChange and Ranks vectors must be same length")

	# we sort them in 3 groups,  the ups, the downs, the rest...
	myRanks <- 1:N
	whoUp <- which( folds > notDE.value)
	whoZero <- which( folds == notDE.value)
	whoDown <- which( folds < notDE.value)

	# now sort just the UPs based on Pvalue and fold and Rank..
	Nup <- length( whoUp)
	rankUp_P <- rankUp_F <- rankUp_R <- vector( length=Nup)
	ordUp_P <- base::order( pvalues[ whoUp], (-folds[whoUp]), decreasing=FALSE)
	ordUp_F <- base::order( folds[ whoUp], (-pvalues[whoUp]), decreasing=TRUE)
	ordUp_R <- base::order( ranks[ whoUp], (-folds[whoUp]), decreasing=FALSE)
	rankUp_P[ ordUp_P] <- 1:Nup
	rankUp_F[ ordUp_F] <- 1:Nup
	rankUp_R[ ordUp_R] <- 1:Nup

	#best rank is the average rank from all 3 criteria
	ordUp <- base::order( ((rankUp_P*wt.pvalues) + (rankUp_F*wt.folds) + 
				(rankUp_R*wt.ranks))/(wt.pvalues+wt.folds+wt.ranks))
	outUp <- whoUp[ ordUp]

	# the zero folds may not all have P=1.0, so sort on those
	ordZero <- base::order( pvalues[ whoZero], ranks[ whoZero])
	outZero <- whoZero[ ordZero]

	# now sort just the DOWNs based on Pvalue and fold and Rank..
	Ndown <- length( whoDown)
	rankDown_P <- rankDown_F <- rankDown_R <- vector( length=Ndown)

	# Meta Results will always have one sided P-values and close to one when fold is negative
	ordDown_P <- base::order( pvalues[ whoDown], -folds[whoDown], decreasing=FALSE)
	ordDown_F <- base::order( folds[ whoDown], -pvalues[whoDown], decreasing=TRUE)
	ordDown_R <- base::order( ranks[ whoDown], -folds[whoDown], decreasing=FALSE)
	rankDown_P[ ordDown_P] <- 1:Ndown
	rankDown_F[ ordDown_F] <- 1:Ndown
	rankDown_R[ ordDown_R] <- 1:Ndown
	#best rank is the average rank from both criteria
	ordDown <- base::order( ((rankDown_P*wt.pvalues) + (rankDown_F*wt.folds) + (rankDown_R*wt.ranks))/(wt.pvalues+wt.folds+wt.ranks))
	outDown <- whoDown[ ordDown]

	out <- c( outUp, outZero, outDown)

	# make sure that we still have every index exactly once
	if ( length( setdiff( out, myRanks)) != 0) warning( "diffExpessRankOrder:  ordering by Fold and Rank problem...")
	return( out)
}

