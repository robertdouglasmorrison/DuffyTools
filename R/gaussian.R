# gaussian.R - a set of probability density functions for peak fitting


gaussian <- function( x, center=0, width=1, height=NULL, floor=0) {

	# adapted from Earl F. Glynn;  Stowers Institute for Medical Research, 2007
	twoVar <- 2 * width * width
	sqrt2piVar <- sqrt( pi * twoVar)
	y <- exp( -( x - center)^2 / twoVar) / sqrt2piVar

	# by default, the height is such that the curve has unit volume
	if ( ! is.null (height)) {
		scalefactor <- sqrt2piVar
		y <- y * scalefactor * height
	}
	y + floor
}


fit.gaussian <- function( x, y, start.center=NULL, start.width=NULL, start.height=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best gaussian to fit the given data

	# make some rough estimates from the values of Y
	who.max <- which.max(y)
	if ( is.null( start.center)) start.center <- x[ who.max]
	if ( is.null( start.height)) start.height <- y[ who.max]
	if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2

	# call the Nonlinear Least Squares, either fitting the floor too or not
	controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	if ( ! fit.floor) {
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
		nlsAns <- try( nls( y ~ gaussian( x, center, width, height), start=starts, control=controlList))
	} else {
		if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height,
				"floor"=start.floor)
		nlsAns <- try( nls( y ~ gaussian( x, center, width, height, floor), start=starts, control=controlList))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		centerAns <- start.center
		widthAns <- start.width
		heightAns <- start.height
		floorAns <- if ( fit.floor) start.floor else 0
		yAns <- gaussian( x, centerAns, widthAns, heightAns, floorAns)
		residualAns <- y - yAns
	} else {
		coefs <-coef(nlsAns)
		centerAns <- coefs[1]
		widthAns <- coefs[2]
		heightAns <- coefs[3]
		floorAns <- if ( fit.floor) coefs[4] else 0
		yAns <- fitted( nlsAns)
		residualAns <- residuals( nlsAns)
	}

	# always report the SD as a possitive value
	widthAns <- abs( widthAns)

	out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}

	return( out)
}


lorentzian <- function( x, center=0, width=1, height=NULL, floor=0) {

	widSq <- width * width
	y <-  width / ( pi * (( x - center)^2 + widSq))

	# by default, the height is such that the curve has unit volume
	if ( ! is.null (height)) {
		scalefactor <- pi * width
		y <- y * scalefactor * height
	}
	y + floor
}


fit.lorentzian <- function( x, y, start.center=NULL, start.width=NULL, start.height=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best lorentzian to fit the given data

	# make some rough estimates from the values of Y
	who.max <- which.max(y)
	if ( is.null( start.center)) start.center <- x[ who.max]
	if ( is.null( start.height)) start.height <- y[ who.max]
	if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2

	# call the Nonlinear Least Squares, either fitting the floor too or not
	controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	if ( ! fit.floor) {
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
		nlsAns <- try( nls( y ~ lorentzian( x, center, width, height), start=starts, control=controlList))
	} else {
		if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height,
				"floor"=start.floor)
		nlsAns <- try( nls( y ~ lorentzian( x, center, width, height, floor), start=starts, control=controlList))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		centerAns <- start.center
		widthAns <- start.width
		heightAns <- start.height
		floorAns <- if ( fit.floor) start.floor else 0
		yAns <- lorentzian( x, centerAns, widthAns, heightAns, floorAns)
		residualAns <- y - yAns
	} else {
		coefs <-coef(nlsAns)
		centerAns <- coefs[1]
		widthAns <- coefs[2]
		heightAns <- coefs[3]
		floorAns <- if ( fit.floor) coefs[4] else 0
		yAns <- fitted( nlsAns)
		residualAns <- residuals( nlsAns)
	}

	# always report the SD as a possitive value
	widthAns <- abs( widthAns)

	out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}

	return( out)
}


gumbel <- function( x, center=0, width=1, height=NULL, floor=0) {

	terms <- ( x - center) / width
	expTerms <- exp( terms)
	y <-  exp( terms - expTerms) / width

	# by default, the height is such that the curve has unit volume
	if ( ! is.null (height)) {
		scalefactor <- exp(1) * width
		y <- y * scalefactor * height
	}
	y + floor
}


fit.gumbel <- function( x, y, start.center=NULL, start.width=NULL, start.height=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best gumbel to fit the given data

	# make some rough estimates from the values of Y
	who.max <- which.max(y)
	if ( is.null( start.center)) start.center <- x[ who.max]
	if ( is.null( start.height)) start.height <- y[ who.max]
	if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2

	# call the Nonlinear Least Squares, either fitting the floor too or not
	controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	if ( ! fit.floor) {
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
		nlsAns <- try( nls( y ~ gumbel( x, center, width, height), start=starts, control=controlList))
	} else {
		if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height,
				"floor"=start.floor)
		nlsAns <- try( nls( y ~ gumbel( x, center, width, height, floor), start=starts, control=controlList))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		centerAns <- start.center
		widthAns <- start.width
		heightAns <- start.height
		floorAns <- if ( fit.floor) start.floor else 0
		yAns <- gumbel( x, centerAns, widthAns, heightAns, floorAns)
		residualAns <- y - yAns
	} else {
		coefs <-coef(nlsAns)
		centerAns <- coefs[1]
		widthAns <- coefs[2]
		heightAns <- coefs[3]
		floorAns <- if ( fit.floor) coefs[4] else 0
		yAns <- fitted( nlsAns)
		residualAns <- residuals( nlsAns)
	}

	# the width for a Gumbel keeps its sign!

	out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}

	return( out)
}
