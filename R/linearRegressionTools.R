# linearRegressionTools.R


`lmSlopeDifference` <- function( lmAns1, lmAns2, coef1=2, coef2=2) {

# Two different citations for the identical formula:
#       formula for the test statistic:  Z = (b1-b2) / sqrt( SEb1^2 + SEb2^2)
#       for 2 slopes 'b1, b2' with standard errors 'SEb1, SEb2'
#
#  1)   http://stats.stackexchange.com/questions/55501/test-a-significant-difference-between-two-slope-values
#	cites: 	Paternoster, R., Brame, R., Mazerolle, P., & Piquero, A. R. (1998). 
#		Using the Correct Statistical Test for the Equality of Regression Coefficients. Criminology, 36(4), 859–866.
#
#  2)   http://www.real-statistics.com/regression/hypothesis-testing-significance-regression-line-slope/
#			comparing-slopes-two-independent-samples/
#	cites:	David Howells textbook entitled Statistical Methods for Psychology, Wadsworth CENGAGE Learning, 2010.
#

	# extract the slope and standard error from the summary objects of two separate liner regression model objects
	lmSummary1 <- summary( lmAns1)
	modelCoefs1 <- coef( lmSummary1)
	lmSummary2 <- summary( lmAns2)
	modelCoefs2 <- coef( lmSummary2)

	# now extract the wanted coefficients (slope and standard error) from each model
	# (in most any LM model, the origin is row #1, and our coefficient of interest is in row #2)
	# (the slope is in column 1, and the standard error is in column 2)
	slope1 <- modelCoefs1[ coef1, 1]
	se1 <- modelCoefs1[ coef1, 2]
	slope2 <- modelCoefs2[ coef2, 1]
	se2 <- modelCoefs2[ coef2, 2]

	# difference in slopes is 'model2 - model1'
	dslope <- slope2 - slope1

	# the test statistic 'Z'
	t.value <- dslope / sqrt( se1*se1 + se2*se2)

	# the combined degrees of freedom
	degFree <- lmSummary1$df[2] + lmSummary2$df[2]

	# calculated the P-value for seeing that large a T statistic ( 2-sided -- 'are the slopes equal')
	p.value <- pt( abs(t.value), df=degFree, lower.tail=FALSE) * 2

	return( list( "difference"=dslope, "t.value"=t.value, "p.value"=p.value))
}

