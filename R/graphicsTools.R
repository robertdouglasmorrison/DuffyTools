# graphicsTools.R -- various drawing routines


# abline that gets clipped to a line segment

abline.segment <- function( a=NULL, b=NULL, reg=NULL, x1, x2, origin=0, log="", ...) {

	intercept <- slope <- NA

	# grab slope / intercept
	if ( ! is.null(a)) intercept <- a
	if ( ! is.null(b)) slope <- b
	if ( ! is.null(reg)) {
		intercept <- coef( reg)[1]
		slope <- coef( reg)[2]
	}
	if (any( is.na( c( intercept, slope)))) stop( "unable to determine slope/intercept from a,b,reg arguments")

	y1 <- slope * (x1-origin) + intercept
	y2 <- slope * (x2-origin) + intercept
	if ( log != "") {
		if ( grepl("x10", log)) {
			y1 <- slope * (log10(x1)-origin) + intercept
			y2 <- slope * (log10(x2)-origin) + intercept
		}
		if ( grepl("x2",log)) {
			y1 <- slope * (log2(x1)-origin) + intercept
			y2 <- slope * (log2(x2)-origin) + intercept
		}
		if ( grepl("y2",log)) {
			y1 <- 2 ^ y1
			y2 <- 2 ^ y2
		}
		if ( grepl("y10",log)) {
			y1 <- 10 ^ y1
			y2 <- 10 ^ y2
		}
	}

	lines( c( x1,x2), c(y1,y2), ...)
}


# pad strings with blank space to make them similar width on plot

padStringWidth <- function( strs, pad=c("right","left")) {

	pad <- match.arg( pad)

	# get how wide all the strings are
	lens <- strwidth( strs)
	lens[ is.na(lens)] <- 0
	bigLen <- max( lens, na.rm=T)
	spaceLen <- strwidth( " ")

	# how many blanks do we need to add?
	nSpaceAdd <- round( (bigLen - lens) / spaceLen)

	# stretch those that need it
	out <- strs
	for ( i in 1:length(strs)) {
		if ( nSpaceAdd[i]) {
			extra <- paste( rep.int( " ", nSpaceAdd[i]), collapse="")
			if ( pad ==  "right") {
				out[i] <- paste( strs[i], extra, sep="")
			} else {
				out[i] <- paste( extra, strs[i], sep="")
			}
		}
	}
	return(out)
}


# my_size_n_color.R -- adapted from package PLOTRIX

my_size_n_color <- function (x = NULL, y, size, sizefun = "sqrt", col, main = "", 
    xlim = NA, xlab = "", xat = NULL, xaxlab = NULL, xcex = 1, 
    xlas = 0, xgrid = FALSE, ylim = NA, ylab = "", yat = NULL, 
    yaxlab = NULL, ycex = 1, ylas = 1, ygrid = TRUE, mar = c(5, 
        4, 4, 2), boxit = TRUE, add = FALSE, border = NA, lwd=1, ...) 
{
    if (!is.na(sizefun)) 
        size <- do.call(sizefun, list(size))
    if (is.matrix(size)) {
        dimsize <- dim(size)
        if (is.null(x)) {
            x <- matrix(rep((1:dimsize[2]) * max(size * 2), each = dimsize[1]), 
                ncol = dimsize[2])
            y <- matrix(rep((dimsize[1]:1) * max(size * 2), dimsize[2]), 
                ncol = dimsize[2])
            if (is.null(xat)) 
                xat <- x[1, ]
            if (is.null(yat)) 
                yat <- y[, 1]
        }
        else {
            if (is.null(xat)) 
                xat <- 1:dimsize[2]
            if (is.null(yat)) 
                yat <- 1:dimsize[1]
        }
    }
    else {
        dimsize = c(length(size), 1)
        if (is.null(x)) 
            x <- 1:length(size)
    }
    xylim <- par("usr")
    aspect_ratio <- (xylim[2] - xylim[1])/(xylim[4] - xylim[3])
    maxsize <- max(size)
    if (is.na(xlim[1])) 
        xlim <- c(min(x) - maxsize * aspect_ratio, max(x) + maxsize * 
            aspect_ratio)
    if (is.na(ylim[1])) 
        ylim <- c(min(y) - maxsize/aspect_ratio, max(y) + maxsize/aspect_ratio)
    if (!add) {
        oldmar <- par(mar = mar)
        plot(x, y, main = main, xlab = xlab, ylab = ylab, xlim = xlim, 
            ylim = ylim, axes = FALSE, type = "n", ...)
    	xylim <- par("usr")
        if (xgrid) 
            segments(xat, xylim[3], xat, xylim[4], col = "lightgray", 
                lty = 2)
        if (ygrid) 
            segments(xylim[1], yat, xylim[2], yat, col = "lightgray", 
                lty = 2)
        axis(1, at = xat, labels = xaxlab, las = xlas, cex.axis = xcex)
        axis(2, at = yat, labels = yaxlab, las = ylas, cex.axis = ycex)
        if (boxit) 
            box()
    }
    if (is.matrix(size)) {
        if (is.null(dim(col))) 
            col = matrix(col, nrow = dimsize[1], ncol = dimsize[2])
        if (is.null(dim(border))) 
            border = matrix(border, nrow = dimsize[1], ncol = dimsize[2])
        for (row in 1:dimsize[1]) {
            for (column in 1:dimsize[2]) draw.circle(x[row, column], 
                y[row, column], size[row, column], border = border[row, column], 
                col = col[row, column], lwd=lwd)
        }
    }
    else {
        if (length(col) < length(size)) 
            col = rep(col, size)
        if (length(border) < length(size)) 
            border = rep(border, size)
        for (index in 1:length(size)) draw.circle(x[index], y[index], 
            size[index], border = border[index], col = col[index])
    }
    if (!add) 
        par(mar = oldmar)
}
