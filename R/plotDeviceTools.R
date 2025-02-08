# plotDeviceTools.R

# a set of routines to operate on plot images


`getPlotDeviceType` <- function( optT="Options.txt") {

	return( getOptionValue( optT, arg="plot.device", notfound="pdf", verbose=F))
}


# write the current plot image to a named file, using hints from the Options.txt file

`printPlot` <- function( filename, device=NULL, width=NULL, height=NULL, optT="Options.txt", ...) {

	if ( is.null(device)) device <- getOptionValue( optT, arg="plot.device", notfound="pdf", verbose=F)
	if ( is.null(width)) width <- getOptionValue( optT, arg="plot.width", notfound="12", verbose=F)
	if ( is.null(height)) height <- getOptionValue( optT, arg="plot.height", notfound="8", verbose=F)

	# allow using the cairo_pdf device when available
	if ( is.character(device) && (device == "pdf") && capabilities("cairo")) device <- "cairo_pdf"

	# turn device name into a function refereence
	dev.FUN <- if ( is.character(device)) get( device) else device
	dev.name <- if ( is.function(device)) as.character( substitute( device)) else device
	if ( ! is.function( dev.FUN)) stop( paste( "Given 'plot.device' is not an R function: ", device))
	width <- as.numeric( width)
	if ( ! is.numeric( width)) stop( paste( "Given 'plot.width' is not numeric: ", width))
	height <- as.numeric( height)
	if ( ! is.numeric( height)) stop( paste( "Given 'plot.height' is not numeric: ", height))

	# Force the filename suffix to match the device type
	dev.name <- sub( "cairo_", "", dev.name)
	suffixPattern <- paste( "\\.", dev.name, "$", sep="")
	if ( ! grepl( suffixPattern, filename)) {
		if ( dev.name == "png") filename <- sub( "\\.pdf$", ".png", filename)
		if ( dev.name == "pdf") filename <- sub( "\\.png$", ".pdf", filename)
	}
	# still not good, just apppend
	if ( ! grepl( suffixPattern, filename)) filename <- paste( filename, dev.name, sep=".")

	# OK, print that 
	dev.print( dev.FUN, filename, width=width, height=height, ...)
	return( filename)
}


# open a new plot device for future plotting, using hints from the Options.txt file

`openPlot` <- function( filename, device=NULL, width=NULL, height=NULL, optT="Options.txt", ...) {

	if ( is.null(device)) device <- getOptionValue( optT, arg="plot.device", notfound="pdf", verbose=F)
	if ( is.null(width)) width <- getOptionValue( optT, arg="plot.width", notfound="12", verbose=F)
	if ( is.null(height)) height <- getOptionValue( optT, arg="plot.height", notfound="8", verbose=F)

	# turn device name into a function refereence
	dev.FUN <- if ( is.character(device)) get( device) else device
	dev.name <- if ( is.function(device)) as.character( substitute( device)) else device
	if ( ! is.function( dev.FUN)) stop( paste( "Given 'plot.device' is not an R function: ", device))
	width <- as.numeric( width)
	if ( ! is.numeric( width)) stop( paste( "Given 'plot.width' is not numeric: ", width))
	height <- as.numeric( height)
	if ( ! is.numeric( height)) stop( paste( "Given 'plot.height' is not numeric: ", height))

	# Force the filename suffix to match the device type
	suffixPattern <- paste( "\\.", device, "$", sep="")
	if ( ! grepl( suffixPattern, filename)) {
		if ( dev.name == "png") filename <- sub( "\\.pdf$", ".png", filename)
		if ( dev.name == "pdf") filename <- sub( "\\.png$", ".pdf", filename)
	}
	# still not good, just apppend
	if ( ! grepl( suffixPattern, filename)) filename <- paste( filename, dev.name, sep=".")

	# OK, open a new plot device.  Caller will close later with dev.off()
	dev.FUN( filename, width=width, height=height, ...)

	return( filename)
}

