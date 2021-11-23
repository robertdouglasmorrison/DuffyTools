# radial.plot.R -- slight tweaks to the 'PLOTRIX' routine

radial.plot <- function (lengths, radial.pos = NULL, labels = NA, label.pos = NULL, 
    radlab = FALSE, start = 0, clockwise = FALSE, rp.type = "r", 
    label.prop = 1.15, main = "", xlab = "", ylab = "", line.col = par("fg"), cex = par( "cex.axis"),
    lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2), 
    show.grid = TRUE, show.grid.labels = 4, show.radial.grid = TRUE, 
    grid.col = "gray", grid.bg = "transparent", grid.left = FALSE, 
    grid.unit = NULL, point.symbols = NULL, point.col = NULL, 
    show.centroid = FALSE, radial.lim = NULL, radial.labels = NULL, 
    boxed.radial = TRUE, poly.col = NULL, add = FALSE, ...) 
{

    # set up args based on defaults
    if (is.null(radial.lim)) radial.lim <- range(lengths)
    length.dim <- dim(lengths)
    if (is.null(length.dim)) {
        npoints <- length(lengths)
        nsets <- 1
        lengths <- matrix(lengths, nrow = 1)
    } else {
        npoints <- length.dim[2]
        nsets <- length.dim[1]
        lengths <- as.matrix(lengths)
    }
    lengths <- lengths - radial.lim[1]
    lengths[lengths < 0] <- NA
    if (is.null(radial.pos[1])) radial.pos <- seq(0, pi * (2 - 2/npoints), length.out = npoints)
    radial.pos.dim <- dim(radial.pos)
    if (is.null(radial.pos.dim)) {
    	radial.pos <- matrix(rep(radial.pos, nsets), nrow = nsets, byrow = TRUE)
    } else {
    	radial.pos <- as.matrix(radial.pos)
    }
    if (clockwise) radial.pos <- -radial.pos
    if (start) radial.pos <- radial.pos + start
    if (show.grid) {
        if (length(radial.lim) < 3) 
            grid.pos <- pretty(radial.lim)
        else grid.pos <- radial.lim
        if (grid.pos[1] < radial.lim[1]) grid.pos <- grid.pos[-1]
        maxlength <- max(grid.pos - radial.lim[1])
        angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
    } else {
        grid.pos <- NA
        maxlength <- diff(radial.lim)
    }

    # we will modify several drawing parameters
    oldpar <- par("xpd", "mar", "pty")
    on.exit( par( oldpar))

    # OK, start to draw
    if (!add) {

    	# first is the radial grid
        par(mar = mar, pty = "s")
        plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
            type = "n", axes = FALSE, main = main, xlab = xlab, 
            ylab = ylab)
        if (show.grid) {
            for (i in seq(length(grid.pos), 1, by = -1)) {
                xpos <- cos(angles) * (grid.pos[i] - radial.lim[1])
                ypos <- sin(angles) * (grid.pos[i] - radial.lim[1])
                polygon(xpos, ypos, border = grid.col, col = grid.bg)
		if (grid.pos[i] == 0) {
			border.col <- "gray10"
                	polygon(xpos, ypos, border=border.col, col=NA, lwd=2, lty=2)
		}
            }
        }
    }

    # now the colored lines
    par(xpd = TRUE)
    if (length(line.col) < nsets) 
        line.col <- 1:nsets
    if (length(rp.type) < nsets) 
        rp.type <- rep(rp.type, length.out = nsets)
    if ( !is.null(point.symbols) && length(point.symbols) < nsets) 
        point.symbols <- rep(point.symbols, length.out = nsets)
    if ( !is.null(point.col) && length(point.col) < nsets) 
        point.col <- rep(point.col, length.out = nsets)
    if ( !is.null(poly.col) && length(poly.col) < nsets) 
        poly.col <- rep(poly.col, length.out = nsets)
    if (length(lty) < nsets) 
        lty <- rep(lty, length.out = nsets)
    if (length(lwd) < nsets) 
        lwd <- rep(lwd, length.out = nsets)
    for (i in 1:nsets) {
        if (nsets > 1) {
            linecol <- line.col[i]
            polycol <- poly.col[i]
            pointcol <- point.col[i]
            pointsymbols <- point.symbols[i]
            ltype <- lty[i]
            lwidth <- lwd[i]
        }
        else {
            linecol <- line.col
            polycol <- poly.col
            pointcol <- point.col
            pointsymbols <- point.symbols
            ltype <- lty
            lwidth <- lwd
        }
        rptype <- unlist(strsplit(rp.type[i], ""))
        if (match("s", rptype, 0)) {
            if (is.null(pointsymbols)) 
                pointsymbols <- i
            if (is.null(pointcol)) 
                pointcol <- i
        }
        xpos <- cos(radial.pos[i, ]) * lengths[i, ]
        ypos <- sin(radial.pos[i, ]) * lengths[i, ]
        if (match("r", rptype, 0)) 
            segments(0, 0, xpos, ypos, col = linecol, lty = ltype, 
                lwd = lwidth, ...)
        if (match("p", rptype, 0)) 
            polygon(xpos, ypos, border = linecol, col = polycol, 
                lty = ltype, lwd = lwidth, ...)
        if (match("s", rptype, 0)) 
            points(xpos, ypos, pch = pointsymbols, col = pointcol, 
                ...)
        if (show.centroid) 
            if (match("p", rptype, 0)) {
                nvertices <- length(xpos)
                polygonarea <- xpos[nvertices] * ypos[1] - xpos[1] * 
                  ypos[nvertices]
                for (vertex in 1:(nvertices - 1)) polygonarea <- polygonarea + 
                  xpos[vertex] * ypos[vertex + 1] - xpos[vertex + 
                  1] * ypos[vertex]
                polygonarea <- polygonarea/2
                centroidx <- (xpos[nvertices] + xpos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                centroidy <- (ypos[nvertices] + ypos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                for (vertex in 1:(nvertices - 1)) {
                  centroidx <- centroidx + (xpos[vertex] + xpos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                  centroidy <- centroidy + (ypos[vertex] + ypos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                }
                points(centroidx/(6 * polygonarea), centroidy/(6 * 
                  polygonarea), col = point.col[i], pch = point.symbols[i], 
                  cex = 2, ...)
            }
            else points(mean(xpos), mean(ypos), col = pointcol, 
                pch = pointsymbols, cex = 2, ...)
    }

    # now the labels for each dimension
    if (!add) {
        if (is.na(labels[1])) {
            label.pos <- seq(0, 1.8 * pi, length = 9)
            labels <- as.character(round(label.pos, 2))
        }
        lablen <- length(labels)
        if (is.null(label.pos[1])) {
            label.pos <- seq(0, pi * (2 - 2/lablen), length.out = lablen)
        }
        if (clockwise) 
            label.pos <- -label.pos
        if (start) 
            label.pos <- label.pos + start
        xpos <- cos(label.pos) * maxlength
        ypos <- sin(label.pos) * maxlength
        if (show.radial.grid) segments(0, 0, xpos, ypos, col = grid.col)

	# draw the grid axis labels before the dimension labels
        if (show.grid.labels) {
            if (show.grid.labels%%2) {
                ypos <- grid.pos - radial.lim[1]
                xpos <- rep(0, length(grid.pos))
                if (show.grid.labels == 1) 
                  ypos <- -ypos
            } else {
                xpos <- grid.pos - radial.lim[1]
                ypos <- rep(0, length(grid.pos))
                if (show.grid.labels == 2) 
                  xpos <- -xpos
            }
            if (is.null(radial.labels)) radial.labels = as.character(grid.pos)
            if (!is.null(grid.unit)) radial.labels[length(grid.pos)] <- paste(radial.labels[length(grid.pos)], grid.unit)
            if (boxed.radial) {
                boxed.labels(xpos, ypos, radial.labels, col=grid.col, border=FALSE, cex=par("cex.lab")*0.9, 
				bg="transparent")
            } else {
	        text(xpos, ypos, radial.labels, col=grid.col, cex=par("cex.lab")*0.9)
	    }
        }

	# let's trim the really long ones
	labels <- clipLongString( labels, max.length=75, pct.front=1.0)

        if (radlab) {
	    # when radial labels, not special shifting
            xpos <- cos(label.pos) * maxlength * label.prop
            ypos <- sin(label.pos) * maxlength * label.prop
            for (label in 1:length(labels)) {
                labelsrt <- (180 * label.pos[label]/pi) + 180 * 
                  (label.pos[label] > pi/2 && label.pos[label] < 
                    3 * pi/2)
                text(xpos[label], ypos[label], labels[label], 
                  cex = cex, srt = labelsrt)
            }
        } else {
	    # when horizontal labels, tweak the positions near pi/2 and 3/2pi...
	    label.prop <- rep.int( label.prop[1], lablen)
	    shiftUp <- which( abs(label.pos) >= (0.35*pi) & abs(label.pos) <= (0.65*pi))
	    if ( length(shiftUp)) label.prop[shiftUp] <- label.prop[shiftUp] + 0.02
	    shiftDown <- which( abs(label.pos) >= (1.35*pi) & abs(label.pos) <= (1.65*pi))
	    if ( length(shiftDown)) label.prop[shiftDown] <- label.prop[shiftDown] + 0.02

	    shiftUp <- which( abs(label.pos) >= (0.4*pi) & abs(label.pos) <= (0.6*pi))
	    if ( length(shiftUp)) label.prop[shiftUp] <- label.prop[shiftUp] + 0.04
	    shiftDown <- which( abs(label.pos) >= (1.4*pi) & abs(label.pos) <= (1.6*pi))
	    if ( length(shiftDown)) label.prop[shiftDown] <- label.prop[shiftDown] + 0.04

	    shiftUp <- which( abs(label.pos) >= (0.46*pi) & abs(label.pos) <= (0.54*pi))
	    if ( length(shiftUp)) label.prop[shiftUp] <- label.prop[shiftUp] + 0.06
	    shiftDown <- which( abs(label.pos) >= (1.46*pi) & abs(label.pos) <= (1.54*pi))
	    if ( length(shiftDown)) label.prop[shiftDown] <- label.prop[shiftDown] + 0.06

            xpos <- cos(label.pos) * maxlength * label.prop
            ypos <- sin(label.pos) * maxlength * label.prop
	    boxed.labels(xpos, ypos, labels, xpad=1.01, ypad = 0.7, border = FALSE, cex = cex,
	    		bg="transparent")
	}
    }
    invisible(oldpar)
}
