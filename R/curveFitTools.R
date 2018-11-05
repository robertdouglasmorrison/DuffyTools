# gaussian.R - a set of probability density functions for peak fitting

normal <- function( x, mean=0, sd=1, height=NULL, floor=0) {

	y <- dnorm( x, mean=mean, sd=sd)
	# by default, the height is such that the curve has unit volume
	if ( ! is.null (height)) {
		scalefactor <- height / max(y)
		y <- y * scalefactor
	}
	y + floor
}


fit.normal <- function( x, y, start.mean=NULL, start.sd=NULL, start.height=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best normal curve to fit the given data
	# make some rough estimates from the values of Y
	who.max <- which.max(y)
	if ( is.null( start.mean)) start.mean <- x[ who.max]
	if ( is.null( start.height)) start.height <- y[ who.max]
	if ( is.null( start.sd)) start.sd <- sum( y > (start.height/2)) / 2

	# call the Nonlinear Least Squares, either fitting the floor too or not
	controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	if ( ! fit.floor) {
		starts <- list( "mean"=start.mean, "sd"=start.sd, "height"=start.height)
		nlsAns <- try( nls( y ~ normal( x, mean, sd, height), start=starts, control=controlList))
	} else {
		if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
		starts <- list( "mean"=start.mean, "sd"=start.sd, "height"=start.height,
				"floor"=start.floor)
		nlsAns <- try( nls( y ~ normal( x, mean, sd, height, floor), start=starts, control=controlList))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		meanAns <- start.mean
		sdAns <- start.sd
		heightAns <- start.height
		floorAns <- if ( fit.floor) start.floor else 0
		yAns <- normal( x, meanAns, sdAns, heightAns, floorAns)
		residualAns <- y - yAns
	} else {
		coefs <-coef(nlsAns)
		meanAns <- coefs[1]
		sdAns <- coefs[2]
		heightAns <- coefs[3]
		floorAns <- if ( fit.floor) coefs[4] else 0
		yAns <- fitted( nlsAns)
		residualAns <- residuals( nlsAns)
	}

	# always report the SD as a positive value
	sdAns <- abs( sdAns)
	out <- list( "mean"=meanAns, "sd"=sdAns, "height"=heightAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}
	return( out)
}


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

	# always report the SD as a positive value
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

	# always report the SD as a positive value
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
	on.exit( options( warn=0)); options( warn= -1);
	on.exit( options( show.error.messages=TRUE)); options( show.error.messages=FALSE);
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


sine <- function( x, period=max(x), phase=0, amplitude=1, floor=0) {

	# scale factor to convert user units to/from radians
	radianFactor <- (2*pi) / period
	# convert to radians
	xr <- x * radianFactor
	phaser <- phase * radianFactor
	# make that sine wave
	y <- sin( xr + phaser)
	# convert to output scaling
	yy <- (y * amplitude) + floor
	yy
}


fit.sine <- function( x, y, start.period=NULL, start.phase=NULL, start.amplitude=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best sine curve to fit the given data
	# make some rough estimates from the values of X and Y
	who.max <- which.max(y)
	who.min <- which.min(y)
	if ( is.null( start.period)) start.period <- abs( x[who.max] - x[who.min]) * 2
	if ( is.null( start.phase)) start.phase <- x[who.max] - start.period/4
	if ( is.null( start.amplitude)) start.amplitude <- y[ who.max]

	# call the Nonlinear Least Squares, either fitting the floor too or not
	controlList <- nls.control( maxiter=200, minFactor=1/1024, warnOnly=TRUE)
	# put lower bounds on things
	lowerBounds <- rep.int( 0, 3)
	if ( ! fit.floor) {
		starts <- list( "period"=start.period, "phase"=start.phase, "amplitude"=start.amplitude)
		nlsAns <- try( nls( y ~ sine( x, period, phase, amplitude), start=starts, control=controlList,
				lower=lowerBounds, algorithm="port"))
	} else {
		if (is.null( start.floor)) start.floor <- median( y, na.rm=T)
		starts <- list( "period"=start.period, "phase"=start.phase, "amplitude"=start.amplitude,
				"floor"=start.floor)
		lowerBounds <- rep.int( 0, 4)
		nlsAns <- try( nls( y ~ sine( x, period, phase, amplitude, floor), start=starts, control=controlList,
				lower=lowerBounds, algorithm="port"))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		periodAns <- start.period
		phaseAns <- start.phase
		amplitudeAns <- start.amplitude
		floorAns <- if ( fit.floor) start.floor else 0
		yAns <- sine( x, periodAns, phaseAns, amplitudeAns, floorAns)
		residualAns <- y - yAns
	} else {
		saveNLS <<- nlsAns
		saveCOEFs <<- coefs <-coef(nlsAns)
		periodAns <- coefs[1]
		phaseAns <- coefs[2]
		amplitudeAns <- coefs[3]
		floorAns <- if ( fit.floor) coefs[4] else 0
		yAns <- fitted( nlsAns)
		residualAns <- residuals( nlsAns)
	}

	out <- list( "period"=periodAns, "phase"=phaseAns, "amplitude"=amplitudeAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}
	return( out)
}


triangleWave <- function( x, period=max(x), phase=0, amplitude=1, floor=0) {

	# scale factor to convert user units to/from radians
	radianFactor <- (2*pi) / period
	# convert to radians
	xr <- x * radianFactor
	phaser <- phase * radianFactor
	# make that triangle wave
	xradians <- xr + phaser
	y <- (2/pi) * asin( sin( xradians))
	# convert to output scaling
	yy <- (y * amplitude) + floor
	yy
}


fit.triangleWave <- function( x, y, start.period=NULL, start.phase=NULL, start.amplitude=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best triangle curve to fit the given data
	# make some rough estimates from the values of X and Y
	who.max <- which.max(y)
	who.min <- which.min(y)
	if ( is.null( start.period)) start.period <- abs( x[who.max] - x[who.min]) * 2
	if ( is.null( start.phase)) start.phase <- x[who.max] - start.period/4
	if ( is.null( start.amplitude)) start.amplitude <- y[ who.max]

	# call the Nonlinear Least Squares, either fitting the floor too or not
	controlList <- nls.control( maxiter=100, minFactor=1/1048, warnOnly=TRUE)
	# put lower bounds on things
	lowerBounds <- rep.int( 0, 3)
	if ( ! fit.floor) {
		starts <- list( "period"=start.period, "phase"=start.phase, "amplitude"=start.amplitude)
		nlsAns <- try( nls( y ~ triangleWave( x, period, phase, amplitude), start=starts, control=controlList,
				lower=lowerBounds, algorithm="port"))
	} else {
		if (is.null( start.floor)) start.floor <- median( y, na.rm=T)
		starts <- list( "period"=start.period, "phase"=start.phase, "amplitude"=start.amplitude,
				"floor"=start.floor)
		lowerBounds <- rep.int( 0, 4)
		nlsAns <- try( nls( y ~ triangleWave( x, period, phase, amplitude, floor), start=starts, 
				control=controlList, lower=lowerBounds, algorithm="port"))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		periodAns <- start.period
		phaseAns <- start.phase
		amplitudeAns <- start.amplitude
		floorAns <- if ( fit.floor) start.floor else 0
		yAns <- triangleWave( x, periodAns, phaseAns, amplitudeAns, floorAns)
		residualAns <- y - yAns
	} else {
		saveNLS <<- nlsAns
		saveCOEFs <<- coefs <-coef(nlsAns)
		periodAns <- coefs[1]
		phaseAns <- coefs[2]
		amplitudeAns <- coefs[3]
		floorAns <- if ( fit.floor) coefs[4] else 0
		yAns <- fitted( nlsAns)
		residualAns <- residuals( nlsAns)
	}

	out <- list( "period"=periodAns, "phase"=phaseAns, "amplitude"=amplitudeAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}
	return( out)
}


pulse <- function( x, center=0, width=1, height=NULL, floor=0) {

	# simple step function -- straight up - over - down
	N <- length(x)
	# assume no peak at all
	yValue <- 1/N
	y <- rep.int( yValue, N)
	# the set of points clearly inside the pulse
	leftEdge <- center - width
	rightEdge <- center + width
	whoLeft <- whoRight <- NA
	who <- which( x >= leftEdge & x <= rightEdge)
	if ( NN <- length(who)) {
		# by default, the height is such that the curve has unit volume
		# so everyone outside is zero
		y <- rep.int( 0, N)
		yValue <- 1/NN
		y[who] <- yValue
		whoLeft <- who[1]
		whoRight <- who[NN]
	}

	# let's do a better job of being cuntinuously differentiable, by giving no-zero values to the "shoulder" points
	if ( ! is.na(whoLeft) && whoLeft > 1) {
		whoEdge <- whoLeft - 1
		edgePct <- x[whoLeft] - leftEdge
		y[whoEdge] <- yValue * edgePct
	}
	if ( ! is.na(whoRight) && whoRight < N) {
		whoEdge <- whoRight + 1
		edgePct <- rightEdge - x[whoRight]
		y[whoEdge] <- yValue * edgePct
	}

	# rescale to force exactly unit volume
	y <- y / sum(y)
	if ( ! is.null( height)) {
		scalefactor <- height / max(y)
		y <- y * scalefactor
	}
	y + floor
}


fit.pulse <- function( x, y, start.center=NULL, start.width=NULL, start.height=NULL,
			start.floor=NULL, fit.floor=FALSE) {

	# try to find the best pulse curve to fit the given data
	# make some rough estimates from the values of Y
	who.max <- which.max(y)
	if ( is.null( start.center)) start.center <- x[ who.max]
	if ( is.null( start.height)) start.height <- y[ who.max]
	if ( is.null( start.width)) start.width <- sum( y > (start.height/2)) / 2

	# call the Nonlinear Least Squares, either fitting the floor too or not
	on.exit( options( warn=0)); options( warn= -1);
	on.exit( options( show.error.messages=TRUE)); options( show.error.messages=FALSE);
	controlList <- nls.control( maxiter=100, minFactor=1/512, warnOnly=TRUE)
	if ( ! fit.floor) {
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height)
		nlsAns <- try( nls( y ~ pulse( x, center, width, height), start=starts, control=controlList))
	} else {
		if (is.null( start.floor)) start.floor <- quantile( y, seq(0,1,0.1))[2]
		starts <- list( "center"=start.center, "width"=start.width, "height"=start.height,
				"floor"=start.floor)
		nlsAns <- try( nls( y ~ pulse( x, center, width, height, floor), start=starts, control=controlList))
	}

	# package up the results to pass back
	if ( class( nlsAns) == "try-error") {
		centerAns <- start.center
		widthAns <- start.width
		heightAns <- start.height
		floorAns <- if ( fit.floor) start.floor else 0
		#cat( "\nDebug:  fit error: ", centerAns, widthAns, heightAns, floorAns, "\n")
		yAns <- pulse( x, centerAns, widthAns, heightAns, floorAns)
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

	# always report the SD as a positive value
	widthAns <- abs( widthAns)
	out <- list( "center"=centerAns, "width"=widthAns, "height"=heightAns, "y"=yAns,
			"residual"=residualAns)
	if ( fit.floor) {
		out <- c( out, "floor"=floorAns)
	}
	return( out)
}

