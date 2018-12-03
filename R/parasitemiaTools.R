# parasitemiaTools.R -- conversion functions for bloodsmear data


`BScountToPctParasitemia` <- function( BScount, WBCdenominator=300, WBCperML=8000, RBCperML=5000000) {

	bsPct <- BScount / WBCdenominator
	if (is.na( WBCperML)) WBCperML <- 8000
	if (is.na( RBCperML)) RBCperML <- 5000000
	rbcFactor <- 1 / RBCperML
	pctParasitemia <- bsPct * WBCperML * rbcFactor * 100
	return( round( pctParasitemia, digits=4))
}


`pctParasitemiaToBScount` <- function( pctParasitemia, WBCdenominator=300, WBCperML=8000, RBCperML=5000000) {

	if (is.na( WBCperML)) WBCperML <- 8000
	if (is.na( RBCperML)) RBCperML <- 5000000
	wbcFactor <- 1 / WBCperML
	bsCnt <- (pctParasitemia / 100) * RBCperML * wbcFactor * WBCdenominator
	useDigits <- ifelse( bsCnt >= 2, 0, 1)
	return( round( bsCnt, digits=useDigits))
}


`BScountPerML` <- function( BScount, WBCdenominator=300, WBCperML=8000) {

	bsPct <- BScount / WBCdenominator
	if (is.na( WBCperML)) WBCperML <- 8000
	perML <- bsPct * WBCperML
	return( round( perML, digits=2))
}

