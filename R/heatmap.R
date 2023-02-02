# heatmap -- slightly modified version of 'heatmap()'


`heatmap` <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, distfun = dist, 
			hclustfun = hclust, reorderfun = function(d, w) reorder(d, w), add.expr, 
			symm = FALSE, revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), 
			na.rm = TRUE, margins = c(5, 5), heatColors = heat.colors(12), ColSideColors, RowSideColors, 
			cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
			colRowLab = par('fg'), colColLab = par('fg'), main = NULL, xlab = NULL, ylab = NULL, 
			RdendroWidth = 1, CdendroHeight = 1, heatWidth = 4, heatHeight = 4,
			show.values=FALSE, show.cex=0.75, show.digits=2,
			keep.dendro = FALSE, verbose = getOption("verbose"), ...) {

    scale <- if (symm && missing(scale)) "none" else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) stop("'x' must be a numeric matrix")

    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1) stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L) stop("'margins' must be a numeric vector of length 2")

    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv")) doCdend <- FALSE
    if (is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                x
            else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc

    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
        if (is.null(rownames(x))) 
            (1L:nr)[rowInd]
        else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
        if (is.null(colnames(x))) 
            (1L:nc)[colInd]
        else colnames(x)
    else labCol[colInd]
    if ( length(colRowLab) == nrow(x)) colRowLab <- colRowLab[rowInd]
    if ( length(colColLab) == ncol(x)) colColLab <- colColLab[colInd]

    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }

    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) RdendroWidth else 0.05, heatWidth)
    lhei <- c((if (doCdend) CdendroHeight else 0.05) + if (!is.null(main)) 0.2 else 0, heatHeight)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1L], 0.2, lhei[2L])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
            1), lmat[, 2] + 1)
        lwid <- c(lwid[1L], 0.2, lwid[2L])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
            "; lmat=\n")
        print(lmat)
    }

    dev.hold()
    on.exit(dev.flush())
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    layout(lmat, widths = lwid, heights = lhei, respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1L], 0, 0, 0.5))
        image(rbind(if (revC) 
            nr:1L
        else 1L:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2L]))
        image(cbind(1L:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1L], 0, 0, margins[2L]))
    if (!symm || scale != "none") 
        x <- t(x)
    if (revC) {
        iy <- nr:1
        if (doRdend) 
            ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1L:nr

    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), 
    	axes = FALSE, xlab = "", ylab = "", col=heatColors, ...)

	if (show.values) {
		for ( ix in 1:nr) text( rep.int(ix,nc), 1:nc, as.character( round(x[ix,],dig=show.digits)), cex=show.cex)
	}
	
    	nColors <- length( myColors <- unique(colColLab))
    	if ( nColors == 1) {
    		axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, cex.axis = cexCol,
    			col.axis=colColLab[1])
	} else {
		mtext( side=1, at=1:nc, text=labCol, las=2, cex=cexCol, col=colColLab)
#		for ( myCol in myColors) {
#			who <- which( colColLab == myCol)
#    			axis(1, who, labels = labCol[who], las = 2, line = -0.5, tick = 0, cex.axis = cexCol,
#    				col.axis=myCol)
#		}
	}
    	if (!is.null(xlab)) 
        	mtext(xlab, side = 1, line = margins[1L] - 1.25)

    	nColors <- length( myColors <- unique(colRowLab))
    	if ( nColors == 1) {
    		axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, cex.axis = cexRow,
    			col.axis=colRowLab[1])
	} else {
		mtext( side=4, at=iy, text=labRow, las=2, cex=cexRow, col=colRowLab)
#		for ( myCol in myColors) {
#			who <- which( colRowLab == myCol)
#    			axis(4, who, labels = labRow[who], las = 2, line = -0.5, tick = 0, cex.axis = cexRow,
#    				col.axis=myCol)
#		}
	}
    	if (!is.null(ylab)) 
        	mtext(ylab, side = 4, line = margins[2L] - 1.25)

    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    par(mar = c(margins[1L], 0, 0, 0))
    if (doRdend) 
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2L]))
    if (doCdend) 
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
        frame()
    if (!is.null(main)) {
        par(xpd = NA)
        title(main, cex.main = 1.5 * op[["cex.main"]], line= -1)
    }
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
        doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}


`heatMapLegend` <- function( colorRamp, at, values, label=NULL, border.lwd=3, ...) {

	par( mai=c(0.1,0.1,0.1,0.1))
	N <- length( colorRamp)
	# make a filled rectangle
	plot( 1,1, type="n", xaxt="n", yaxt="n", xlim=c(0,N),
		ylim=c(0,3), xlab=NA, ylab=NA, axes=F)

	for( i in 1:N) rect( i-1, 1, i, 2, col=colorRamp[i], border=colorRamp[i])
	rect( 0, 1, N, 2, border=1, lwd=border.lwd)
	text( at, 1, values, pos=1, ...)
	if ( ! is.null(label)) text( N/2, 2, label, pos=3)
}


`correlationHeatmap` <- function( m, Rowv=NA, Colv=NA, show.values=TRUE, show.digits=2, show.cex=0.75, ...) {

	# given a square matrix of correlation coefficients (-1...+1)
	require( gplots)
	require( heatmap.plus)

	# hard wire to a read green color ramp
	colorRamp <- rainbow( 201, start=0.05, end=0.4, cosineShiftMagnitude=0, rev=T)
	colorRampValues <- seq( -1, 1, by=0.01)

	# force CC matrix to -1/+1 just in case
	m[ m > 1] <- 1
	m[ m < -1] <- -1
	
	# when not using row dendrograms, invert to make it traditional with cell 1,1 in upper left
	if ( is.na( Rowv)) {
		mtmp <- m
		for (i in 1:ncol(m)) mtmp[,i] <- rev( m[,i])
		rownames(mtmp) <- rev( rownames(m))
		m <- mtmp
	}
	
	# plot that heat
	ans <- heatmap( m, heatColors=colorRamp, Rowv=Rowv, Colv=Colv, scale="none", 
			show.value=show.values, show.digits=show.digits, show.cex=show.cex, ...)
	return( invisible(ans))
}

