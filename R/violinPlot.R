# violinPlot.R -- various attempts at violin plots


# try 1:  wrapper around the ggplot geom_violin plot function
`violinPlot` <- function( df, aes, ..., horiz=FALSE, boxwid=NULL, facet=NULL,
			main="", xlab="X", ylab="Y", log="") {

	require( ggplot2)

	p <- ggplot( df, aes) + geom_violin( ...)
	if (horiz) p <- p + coord_flip()
	if ( ! is.null( facet)) {
		p <- p + facet_wrap( paste( "~", as.character(facet)), scales="free_x")
	}

	# add boxplot?
	if ( ! is.null(boxwid)) {
		p <- p + geom_boxplot( width=boxwid)
	}

	# any X,Y log scaling?
	if ( regexpr( 'Y', toupper(log)) > 0) p <- p + scale_y_log10()
	if ( regexpr( 'X', toupper(log)) > 0) p <- p + scale_x_log10()

	p <- p + labs( title=main, x=xlab, y=ylab)
	p
}


# try 2::  a more generic version, starting from the code at UsingR, by John Verzani <verzani@math.csi.cuny.edu>
# 	   original courtesy of boxplot from base

`violinplot` <- function(x, ...) UseMethod("violinplot")

`violinplot.default` <- function(x, orientation=c("vertical","horizontal"), at=NULL, viowid=0.8,
		  		 boxwid=0.2, names=NULL, pars=NULL, col=NULL, border=par('fg'),
		  		 n.density=128, cut.density=0, bw.density="nrd0",
		  		 xlim=NULL, ylim=NULL, log="", main=NULL, sub=NULL, 
		  		 xlab=NULL, ylab=NULL, ann=par("ann"), axes=TRUE, frame.plot=axes, 
				 pch=NULL, lwd=1, lty=1, ...) {
	 
	 # deduce which args get passed where	 		 
	args <- list(x, ...)
	namedargs <- if ( !is.null( attributes(args)$names)) (attributes(args)$names != "") 
			else rep(FALSE, length = length(args))
	pars <- c( args[namedargs], pars)

	# turn what we were given into the groups to plot
	groups <- if ( is.list(x)) x else if ( is.matrix(x)) as.list(as.data.frame(x)) 
			else if ( is.data.frame(x)) as.list(x) else args[!namedargs]
	if ( ! (nGroups <- length(groups))) stop("invalid first argument")
	if ( length(class(groups))) groups <- unclass(groups)
	if ( !missing(names)) attr(groups, "names") <- names else {
		if ( is.null(attr(groups, "names"))) attr(groups, "names") <- 1:nGroups
		names <- attr(groups, "names")
	}

	## set up storage for the density for each group
	xvals <- matrix( NA, nrow=n.density, ncol=nGroups)
	yvals <- matrix( NA, nrow=n.density, ncol=nGroups)
	if ( is.null(at)) center <- 1:nGroups else {
		center <- as.numeric( at)[1:nGroups]
		if ( any( is.na(center))) stop( "length of 'at' must match number of groups")
	}
	
	orientation <- match.arg( orientation)
	needLog <- (orientation == "vertical" && log == "y") || (orientation == "horizontal" && log == "x")
	for(i in 1:nGroups) {
		myX <- groups[[i]]
		myX <- myX[ !is.na(myX)]
		if ( ! length(myX)) next
		if ( needLog) myX <- log10(myX)
		tmp.dens <- density( myX, bw=bw.density, cut=cut.density, n=n.density)
		xvals[,i] <- tmp.dens$x
		if (needLog) xvals[,i] <- 10 ^ tmp.dens$x
		yvals.needtoscale <- tmp.dens$y
		## scale the density curve to not touch neighbors
		yvals.scaled <- (viowid/2) * yvals.needtoscale / max(yvals.needtoscale)
		yvals[,i] <- yvals.scaled
	}

	## now ready to plot
	## need to first make the plot range, depends on horizontal or vertical
	if ( orientation == "vertical") {
		xrange <- range(center) + c(-0.5,0.5)
		yrange <- range(xvals)
	} else {
		xrange <- range(xvals)
		yrange <- range(center) + c(-0.5,0.5) 
	}
	if ( is.null(xlim)) xlim <- xrange
	if ( is.null(ylim)) ylim <- yrange
	
	# make that new plot window, borrowed from 'plot.default()'
	localAxis <- function(..., col, bg, pch, cex, lty, lwd) Axis(...)
	localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
	localWindow <- function(xlim, ylim, log, ..., col, bg, pch, cex, lty, lwd) plot.window( xlim, ylim, log, ...)
	localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
	dev.hold()
	on.exit(dev.flush())
	plot.new()
	localWindow( xlim=xlim, ylim=ylim, log=log, ...)
	if (frame.plot) localBox(...)
	if (ann) localTitle(main = main, sub = sub, xlab = xlab, ylab = ylab, ...)

	# ready to make actual violins
	border <- rep( border, length.out=nGroups)
	if ( is.null(col)) col <- par('bg')
	col <- rep( col, length.out=nGroups)
	if ( ! is.null(pch)) pch <- rep( pch, length.out=nGroups)
	lty <- rep( lty, length.out=nGroups)

	# store the location data to return
	vioPts <- matrix( NA, nrow=5, ncol=nGroups)
	for (i in 1:nGroups) {
		if ( all( is.na( xvals[,i]))) next
		ans1 <- vlnplt( xvals[,i], yvals[,i], center[i], border=border[i], col=col[i],
				 orientation=orientation, lwd=lwd, lty=lty[i], ...)
		vioPts[1,i] <- ans1[1]
		vioPts[5,i] <- ans1[2]
		if (boxwid > 0 || !is.null(pch)) {
			ans2 <- vlnbxp( groups[[i]], center[i], border=border[i], width=boxwid,
				 	orientation=orientation, pch=if (is.null(pch)) NULL else pch[i], 
					lwd=lwd, lty=lty[i], ...)
			vioPts[2:4,i] <- ans2
		}
	}
	
	if ( orientation == "vertical") {	
		axis( 1, at=center, labels=names, ...)
		axis( 2, ...)
	} else {
		axis( 1, ...)
		axis( 2, at=center, labels=names, ...)
	}
	return( invisible( vioPts))
}		  


##' formula interface for violinplot
`violinplot.formula` <- function(formula, data = NULL, ..., subset) {

	if(missing(formula) || (length(formula) != 3)) stop("formula missing or incorrect")
	m <- match.call(expand.dots = FALSE)
	if(is.matrix(eval(m$data, parent.frame()))) m$data <- as.data.frame(data)
	m$... <- NULL
	m[[1]] <- as.name("model.frame")
	mf <- eval(m, parent.frame())
	response <- attr(attr(mf, "terms"), "response")
	violinplot( split( mf[[response]], mf[-response]), ...)
}


## Make a simple violin plot call from violinplot. values are x,y to plot
`vlnplt` <- function( x, y, center, orientation=c("vertical","horizontal"),
			col=NA , border='black', ...) {

	orientation <- match.arg( orientation)

	## double up first
	x <- c(x,rev(x))
	y <- c(y,-rev(y))
	y <- y + center 
	if (orientation == "vertical") {
		# swtich x and y
		tmp=x; x=y; y=tmp;
	}
	polygon( x, y, border=border, col=col, ...)
	return( range(x))
}


## Make a simple boxplot on the violin
`vlnbxp` <- function( x, center, orientation=c("vertical","horizontal"),
 			border='black', width=0.2, lwd=1, lty=1, pch=NULL, ...) {

	orientation <- match.arg( orientation)
	hw <- width / 2
	## get the box stats, and crop them to the window
	z <- boxplot.stats(x)$stats
	usr <- par('usr')
	if ( par('xlog')) usr[1:2] <- 10 ^ usr[1:2]
	if ( par('ylog')) usr[3:4] <- 10 ^ usr[3:4]
	if (orientation == "vertical") {
		z[1] <- max( z[1], usr[3])
		z[2] <- max( z[2], usr[3])
		z[3] <- max( z[3], usr[3])
		z[3] <- min( z[3], usr[4])
		z[4] <- min( z[4], usr[4])
		z[5] <- min( z[5], usr[4])
	} else {
		z[1] <- max( z[1], usr[1])
		z[2] <- max( z[2], usr[1])
		z[3] <- max( z[3], usr[1])
		z[3] <- min( z[3], usr[2])
		z[4] <- min( z[4], usr[2])
		z[5] <- min( z[5], usr[2])
	}

	if ( hw > 0) {
		if (orientation == "vertical") {
			rect( center-hw, z[2], center+hw, z[4], border=border, lwd=lwd, lty=lty, ...)
			lines( c( center, center, NA, center, center), c( z[1], z[2], NA, z[4], z[5]), 
					col=border, lwd=lwd, lty=lty, ...)
			lines( c(center-hw, center+hw), c( z[3], z[3]), lwd=lwd*3, lty=lty, col=border, ...)
		} else {
			rect( z[2], center-hw, z[4], center+hw, border=border, lwd=lwd, lty=lty, ...)
			lines( c( z[1], z[2], NA, z[4], z[5]), c( center, center, NA, center, center), 
					col=border, lwd=lwd, lty=lty, ...)
			lines( c( z[3], z[3]), c(center-hw, center+hw), lwd=lwd*3, lty=lty, col=border, ...)
		}
	}
	if ( ! is.null( pch)) {
		ctr <- rep.int( center, length(x))
		if (orientation == "vertical") {
			points( jitter(ctr,amount=hw*0.5), x, pch=pch, col=border, ...)
		} else {
			points( x, jitter(ctr,amount=hw*0.5), pch=pch, col=border, ...)
		}
	}
	return( z[2:4])
}  
 
