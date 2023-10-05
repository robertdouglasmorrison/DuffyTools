# survivalTools.R -- implement Kaplan-Meier curves

`kaplanMeier` <- function( groups, times, outcomes, eventOutcome="Y", col=2:(length(unique(groups))+1),
			makePlot=TRUE, xlab="Time  (weeks)", ylab="Survival", lwd=3, lty=1, legend.cex=1.1,
			main="Kaplan-Meier: ", nFDR=0, legend.bty="o", 
			xscale=1, yscale=1, mark.time=FALSE, pch=3, cumhaz=FALSE, ylim=c(0,1.03),
			xstagger=0, ystagger=0, show.pvalue=TRUE, legend.loc="topright", ...) {

	require( survival)

	# turn the outcome calls and followup time into censored survival data
	surv <- Surv( time=times, event=(outcomes == eventOutcome))
	survAns <- survfit( surv ~ groups)
	survDiff <- survdiff( surv ~ groups)

	# calc the P-value
	pval <- pchisq( survDiff$chisq, df=1, lower=F)

	# gather details for plotting
	grpNames <- if (is.factor(groups)) levels(groups) else sort( unique( groups))

	if (makePlot) {
		plotSurvFit( survAns, mark.time=mark.time, pch=pch, xscale=xscale, yscale=yscale, cumhaz=cumhaz,
			col=col, lwd=lwd, xlab=xlab, ylab=ylab, main=main, ylim=ylim, lty=lty, 
			xstagger=xstagger, ystagger=ystagger, ...)
		if (nchar(legend.loc)) legend( legend.loc, grpNames, col=col, lwd=lwd, lty=lty, bg='white', cex=legend.cex, bty=legend.bty)
		if (show.pvalue) legend( 'bottomleft', paste( "P-value =", round( pval, digits=4)), bg='white', cex=max(1,legend.cex),
			bty=legend.bty)
	}

	# perhaps doa FDR permutation test
	if (nFDR > 0) {
		pFDR <- vector( length=nFDR)
		N <- length(groups)
		for ( i in 1:nFDR) {
			fdr.groups <- groups[ sample(N)]
			survDiff2 <- survdiff( surv ~ fdr.groups)
			pFDR[i] <- pchisq( survDiff2$chisq, df=1, lower=F)
		}
		myFDR <- sum( pFDR <= pval) / nFDR
		legend( 'bottomright', paste( "FDR (Trials=",nFDR,") = ", round( myFDR, digits=4), sep=""), 
			bg='white', cex=legend.cex)
	}

	out <- list( "surv.diff"=survDiff, "p.value"=pval)
	return( invisible( out))
}



`plotSurvFit` <- function (x, conf.int, mark.time = FALSE, pch = 3, col = 1, lty = 1, 
    lwd = 1, cex = 1, log = FALSE, xscale = 1, yscale = 1, xlim, 
    ylim, xmax, fun, xlab = "", ylab = "", xaxs = "r", conf.times, 
    conf.cap = 0.005, conf.offset = 0.012, conf.type = c("log", 
        "log-log", "plain", "logit", "arcsin"), mark, noplot = "(s0)", 
    cumhaz = FALSE, firstx, ymin, xstagger=0, ystagger=0, ...) 
{
    dotnames <- names(list(...))
    if (any(dotnames == "type")) 
        stop("The graphical argument 'type' is not allowed")
    x <- survfit0(x, x$start.time)
    if (is.logical(log)) {
        ylog <- log
        xlog <- FALSE
        if (ylog) 
            logax <- "y"
        else logax <- ""
    }
    else {
        ylog <- (log == "y" || log == "xy")
        xlog <- (log == "x" || log == "xy")
        logax <- log
    }
    if (!missing(fun)) {
        if (is.character(fun)) {
            if (fun == "log" || fun == "logpct") 
                ylog <- TRUE
            if (fun == "cloglog") {
                xlog <- TRUE
                if (ylog) 
                  logax <- "xy"
                else logax <- "x"
            }
        }
    }
    if (missing(conf.int) && missing(conf.times)) 
        conf.int <- (!is.null(x$std.err) && prod(dim(x) == 1))
    if (missing(conf.times)) 
        conf.times <- NULL
    else {
        if (!is.numeric(conf.times)) 
            stop("conf.times must be numeric")
        if (missing(conf.int)) 
            conf.int <- TRUE
    }
    if (!missing(conf.int)) {
        if (is.numeric(conf.int)) {
            conf.level <- conf.int
            if (conf.level < 0 || conf.level > 1) 
                stop("invalid value for conf.int")
            if (conf.level == 0) 
                conf.int <- FALSE
            else if (conf.level != x$conf.int) {
                x$upper <- x$lower <- NULL
            }
            conf.int <- TRUE
        }
        else conf.level = 0.95
    }
    stime <- x$time
    std <- NULL
    yzero <- FALSE
    smat <- function(x) {
        dd <- dim(x)
        if (is.null(dd)) 
            as.matrix(x)
        else if (length(dd) == 2) 
            x
        else matrix(x, nrow = dd[1])
    }
    if (cumhaz) {
        if (is.null(x$cumhaz)) 
            stop("survfit object does not contain a cumulative hazard")
        if (is.numeric(cumhaz)) {
            dd <- dim(x$cumhaz)
            if (is.null(dd)) 
                nhazard <- 1
            else nhazard <- prod(dd[-1])
            if (cumhaz != floor(cumhaz)) 
                stop("cumhaz argument is not integer")
            if (any(cumhaz < 1 | cumhaz > nhazard)) 
                stop("subscript out of range")
            ssurv <- smat(x$cumhaz)[, cumhaz, drop = FALSE]
            if (!is.null(x$std.chaz)) 
                std <- smat(x$std.chaz)[, cumhaz, drop = FALSE]
        }
        else if (is.logical(cumhaz)) {
            ssurv <- smat(x$cumhaz)
            if (!is.null(x$std.chaz)) 
                std <- smat(x$std.chaz)
        }
        else stop("invalid cumhaz argument")
    }
    else if (inherits(x, "survfitms")) {
        i <- (x$states != noplot)
        if (all(i) || !any(i)) {
            ssurv <- smat(x$pstate)
            if (!is.null(x$std.err)) 
                std <- smat(x$std.err)
        }
        else {
            i <- which(i)
            if (length(dim(x$pstate)) == 3) {
                ssurv <- smat(x$pstate[, , i, drop = FALSE])
                if (!is.null(x$std.err)) 
                  std <- smat(x$std.err[, , i, drop = FALSE])
            }
            else {
                ssurv <- x$pstate[, i, drop = FALSE]
                if (!is.null(x$std.err)) 
                  std <- x$std.err[, i, drop = FALSE]
            }
        }
    }
    else {
        yzero <- TRUE
        ssurv <- as.matrix(x$surv)
        if (!is.null(x$std.err)) 
            std <- as.matrix(x$std.err)
        if (!missing(fun) && is.character(fun) && fun == "cumhaz") {
            cumhaz <- TRUE
            if (!is.null(x$cumhaz)) {
                ssurv <- as.matrix(x$cumhaz)
                if (!is.null(x$std.chaz)) 
                  std <- as.matrix(x$std.chaz)
            }
            else {
                ssurv <- as.matrix(-log(x$surv))
                if (!is.null(x$std.err)) {
                  if (x$logse) 
                    std <- as.matrix(x$std.err)
                  else std <- as.matrix(x$std.err/x$surv)
                }
            }
        }
    }
    if (is.null(x$strata)) {
        nstrat <- 1
        stemp <- rep(1, length(x$time))
    }
    else {
        nstrat <- length(x$strata)
        stemp <- rep(1:nstrat, x$strata)
    }
    ncurve <- nstrat * ncol(ssurv)
    conf.type <- match.arg(conf.type)
    if (conf.type == "none") 
        conf.int <- FALSE
    if (conf.int == "none") 
        conf.int <- FALSE
    if (conf.int == "only") {
        plot.surv <- FALSE
        conf.int <- TRUE
    }
    else plot.surv <- TRUE
    if (conf.int) {
        if (is.null(std)) 
            stop("object does not have standard errors, CI not possible")
        if (cumhaz) {
            if (missing(conf.type)) 
                conf.type = "plain"
            temp <- survfit_confint(ssurv, std, logse = FALSE, 
                conf.type, conf.level, ulimit = FALSE)
            supper <- as.matrix(temp$upper)
            slower <- as.matrix(temp$lower)
        }
        else if (is.null(x$upper)) {
            if (missing(conf.type) && !is.null(x$conf.type)) 
                conf.type <- x$conf.type
            temp <- survfit_confint(ssurv, x$std.err, logse = x$logse, 
                conf.type, conf.level, ulimit = FALSE)
            supper <- as.matrix(temp$upper)
            slower <- as.matrix(temp$lower)
        }
        else {
            supper <- as.matrix(x$upper)
            slower <- as.matrix(x$lower)
        }
    }
    else supper <- slower <- NULL
    if (!inherits(x, "survfitms") && !cumhaz & !missing(fun)) {
        yzero <- FALSE
        if (is.character(fun)) {
            tfun <- switch(tolower(fun), log = function(x) x, 
                event = function(x) 1 - x, cumhaz = function(x) -log(x), 
                cloglog = function(x) log(-log(x)), pct = function(x) x * 
                  100, logpct = function(x) 100 * x, identity = function(x) x, 
                f = function(x) 1 - x, s = function(x) x, surv = function(x) x, 
                stop("Unrecognized function argument"))
            if (tolower(fun) %in% c("identity", "s") && !inherits(x, 
                "survfitms") && !cumhaz) 
                yzero <- TRUE
        }
        else if (is.function(fun)) 
            tfun <- fun
        else stop("Invalid 'fun' argument")
        ssurv <- tfun(ssurv)
        if (!is.null(supper)) {
            supper <- tfun(supper)
            slower <- tfun(slower)
        }
    }
    if (missing(mark.time) & !missing(mark)) 
        mark.time <- TRUE
    if (missing(pch) && !missing(mark)) 
        pch <- mark
    if (length(pch) == 1 && is.character(pch)) 
        pch <- strsplit(pch, "")[[1]]
    pch <- rep(pch, length.out = ncurve)
    mcol <- rep(col, length.out = ncurve)
    if (is.numeric(mark.time)) 
        mark.time <- sort(mark.time)
    if (conf.int & is.null(conf.times)) {
        if (length(lty) == 1 && is.numeric(lty)) 
            lty <- rep(c(lty, lty + 1, lty + 1), ncurve)
        else if (length(lty) <= ncurve) 
            lty <- rep(rep(lty, each = 3), length.out = (ncurve * 
                3))
        else lty <- rep(lty, length.out = ncurve * 3)
        if (length(col) <= ncurve) 
            col <- rep(rep(col, each = 3), length.out = 3 * ncurve)
        else col <- rep(col, length.out = 3 * ncurve)
        if (length(lwd) <= ncurve) 
            lwd <- rep(rep(lwd, each = 3), length.out = 3 * ncurve)
        else lwd <- rep(lwd, length.out = 3 * ncurve)
    }
    else {
        col <- rep(col, length.out = ncurve)
        lty <- rep(lty, length.out = ncurve)
        lwd <- rep(lwd, length.out = ncurve)
    }
    if (!missing(xlim)) {
        if (!missing(xmax)) 
            warning("cannot have both xlim and xmax arguments, xmax ignored")
        if (!missing(firstx)) 
            stop("cannot have both xlim and firstx arguments")
    }
    if (!missing(ylim)) {
        if (!missing(ymin)) 
            stop("cannot have both ylim and ymin arguments")
    }
    if (!missing(xlim) && !is.null(xlim)) {
        tempx <- xlim
        if (xaxs == "S") 
            tempx[2] <- tempx[2] + diff(tempx) * 1.04
    }
    else {
        temp <- stime[is.finite(stime)]
        if (!missing(xmax) && missing(xlim)) 
            temp <- temp[temp <= xmax]
        if (xaxs == "S") {
            if (xlog) 
                tempx <- c(min(temp[temp > 0]), max(temp) + diff(temp) * 
                  1.04)
            else tempx <- c(min(temp), max(temp)) * 1.04
        }
        else if (xlog) 
            tempx <- range(temp[temp > 0])
        else tempx <- range(temp)
    }
    if (!missing(ylim) && !is.null(ylim)) 
        tempy <- ylim
    else {
        skeep <- is.finite(stime) & stime >= tempx[1] & stime <= 
            tempx[2]
        if (ylog) {
            if (!is.null(supper)) 
                tempy <- range(c(slower[is.finite(slower) & slower > 
                  0 & skeep], supper[is.finite(supper) & skeep]))
            else tempy <- range(ssurv[is.finite(ssurv) & ssurv > 
                0 & skeep])
            if (tempy[2] == 1) 
                tempy[2] <- 0.99
            if (any(c(ssurv, slower)[skeep] == 0)) {
                tempy[1] <- tempy[1] * 0.8
                ssurv[ssurv == 0] <- tempy[1]
                if (!is.null(slower)) 
                  slower[slower == 0] <- tempy[1]
            }
        }
        else {
            if (!is.null(supper)) 
                tempy <- range(c(supper[skeep], slower[skeep]), 
                  finite = TRUE, na.rm = TRUE)
            else tempy <- range(ssurv[skeep], finite = TRUE, 
                na.rm = TRUE)
            if (yzero) 
                tempy <- range(c(0, tempy))
        }
    }
    if (!missing(ymin)) 
        tempy[1] <- ymin
    temp <- if (xaxs == "S") 
        "i"
    else xaxs
    plot(range(tempx, finite = TRUE, na.rm = TRUE)/xscale, range(tempy, 
        finite = TRUE, na.rm = TRUE) * yscale, type = "n", log = logax, 
        xlab = xlab, ylab = ylab, xaxs = temp, ...)
    if (yscale != 1) {
        if (ylog) 
            par(usr = par("usr") - c(0, 0, log10(yscale), log10(yscale)))
        else par(usr = par("usr")/c(1, 1, yscale, yscale))
    }
    if (xscale != 1) {
        if (xlog) 
            par(usr = par("usr") - c(log10(xscale), log10(xscale), 
                0, 0))
        else par(usr = par("usr") * c(xscale, xscale, 1, 1))
    }
    if (xaxs == "i") 
        resetclip <- FALSE
    else resetclip <- !(missing(xlim) & missing(ylim) & missing(xmax) & 
        missing(firstx) & missing(ymin))
    if (resetclip) {
        if (xaxs == "S") 
            tempx <- c(tempx[1], temp[1])
        clip(tempx[1], tempx[2], tempy[1], tempy[2])
        options(plot.survfit = list(plotclip = c(tempx, tempy), 
            plotreset = par("usr")))
    }
    else options(plot.survfit = NULL)
    dostep <- function(x, y) {
        keep <- is.finite(x) & is.finite(y)
        if (!any(keep)) 
            return()
        if (!all(keep)) {
            x <- x[keep]
            y <- y[keep]
        }
        n <- length(x)
        if (n == 1) 
            list(x = x, y = y)
        else if (n == 2) 
            list(x = x[c(1, 2, 2)], y = y[c(1, 1, 2)])
        else {
            temp <- rle(y)$lengths
            drops <- 1 + cumsum(temp[-length(temp)])
            if (n %in% drops) {
                xrep <- c(x[1], rep(x[drops], each = 2))
                yrep <- rep(y[c(1, drops)], c(rep(2, length(drops)), 
                  1))
            }
            else {
                xrep <- c(x[1], rep(x[drops], each = 2), x[n])
                yrep <- c(rep(y[c(1, drops)], each = 2))
            }
            list(x = xrep, y = yrep)
        }
    }
    drawmark <- function(x, y, mark.time, censor, cex, ...) {
        if (!is.numeric(mark.time)) {
            xx <- x[censor > 0]
            yy <- y[censor > 0]
            if (any(censor > 1)) {
                j <- pmax(1, which(censor > 1) - 1)
                i <- censor[censor > 0]
                yy[i > 1] <- (yy[i > 1] + y[j])/2
            }
        }
        else {
            xx <- mark.time
            yy <- approx(x, y, xx, method = "constant", f = 0)$y
        }
        points(xx, yy, cex = cex, ...)
    }
    type <- "s"
    c1 <- 1
    c2 <- 1
    xend <- yend <- double(ncurve)
    if (length(conf.offset) == 1) 
        temp.offset <- (1:ncurve - (ncurve + 1)/2) * conf.offset * 
            diff(par("usr")[1:2])
    else temp.offset <- rep(conf.offset, length = ncurve) * diff(par("usr")[1:2])
    temp.cap <- conf.cap * diff(par("usr")[1:2])
    for (j in 1:ncol(ssurv)) {
        for (i in unique(stemp)) {
            who <- which(stemp == i)
            if (is.null(x$n.censor)) 
                censor <- ifelse(x$n.event[who] == 0, 1, 0)
            else censor <- ifelse(x$n.censor[who] == 0, 0, 1 + 
                (x$n.event[who] > 0))
            xx <- stime[who]
            yy <- ssurv[who, j]
	    if ( any( c(xstagger,ystagger) != 0)) {
	    	mysign <- if (cumhaz) 1 else -1
		xjit <- (i-1) * xstagger * mysign
		yjit <- (i-1) * ystagger * mysign
		xx <- xx + xjit
		yy <- yy + yjit
	    }
            if (plot.surv) {
                if (type == "s") 
                  lines(dostep(xx, yy), lty = lty[c2], col = col[c2], 
                    lwd = lwd[c2])
                else lines(xx, yy, type = type, lty = lty[c2], 
                  col = col[c2], lwd = lwd[c2])
                if (is.numeric(mark.time) || mark.time) 
                  drawmark(xx, yy, mark.time, censor, pch = pch[c1], 
                    col = mcol[c1], cex = cex)
            }
            xend[c1] <- max(xx)
            yend[c1] <- yy[length(yy)]
            if (conf.int && !is.null(conf.times)) {
                x2 <- conf.times + temp.offset[c1]
                templow <- approx(xx, slower[who, j], x2, method = "constant", 
                  f = 1)$y
                temphigh <- approx(xx, supper[who, j], x2, method = "constant", 
                  f = 1)$y
                segments(x2, templow, x2, temphigh, lty = lty[c2], 
                  col = col[c2], lwd = lwd[c2])
                if (conf.cap > 0) {
                  segments(x2 - temp.cap, templow, x2 + temp.cap, 
                    templow, lty = lty[c2], col = col[c2], lwd = lwd[c2])
                  segments(x2 - temp.cap, temphigh, x2 + temp.cap, 
                    temphigh, lty = lty[c2], col = col[c2], lwd = lwd[c2])
                }
            }
            c1 <- c1 + 1
            c2 <- c2 + 1
            if (conf.int && is.null(conf.times)) {
                if (type == "s") {
                  lines(dostep(xx, slower[who, j]), lty = lty[c2], 
                    col = col[c2], lwd = lwd[c2])
                  c2 <- c2 + 1
                  lines(dostep(xx, supper[who, j]), lty = lty[c2], 
                    col = col[c2], lwd = lwd[c2])
                  c2 <- c2 + 1
                }
                else {
                  lines(xx, slower[who, j], lty = lty[c2], col = col[c2], 
                    lwd = lwd[c2], type = type)
                  c2 <- c2 + 1
                  lines(xx, supper[who, j], lty = lty[c2], col = col[c2], 
                    lwd = lwd[c2], type = type)
                  c2 <- c2 + 1
                }
            }
        }
    }
    lastx <- list(x = xend, y = yend)
    if (resetclip) {
        xx <- par("usr")
        clip(xx[1], xx[2], xx[3], xx[4])
    }
    invisible(lastx)
}

