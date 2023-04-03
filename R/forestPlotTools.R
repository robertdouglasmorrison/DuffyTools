# forestPlotTools.R -- methods to generate 'forest plot' that shows deviation from zero / chance


forestPlot <- function( )  {

	leftLabelLoc <- text1Loc <- text2Loc <- text3Loc <- pvalLoc <- 0
	xDataLim <- yMax <- text.cex <- 0
	textTable <- matrix( "", 1,5)
	annotateText <- TRUE


	setup <- function( xRange, yBig, main="Forest Plot", noDifference=0, symmetricX=TRUE, text.cex=0.75,
				annotateText=TRUE, meanDiffMode=TRUE, sub=NULL, dividerLines=FALSE) {
	
		# given the numeric range of values to be drawn, set up our bounds}
		myXlim <- range( xRange)
		if (symmetricX) {
			bigX <- max( abs(xRange))
			myXlim <- c( -bigX, bigX)
		}
		# because R stretches the plot limits a bit, do that here too for clipping later
		xDataLim <<- myXlim * 1.025

		# place all the items we will draw relative to this
		if (meanDiffMode) {
			xWidth <- diff( myXlim)
			leftLabelLoc <<- myXlim[1] - xWidth
			text1Loc <<- myXlim[1] - xWidth*0.68
			text2Loc <<- myXlim[1] - xWidth*0.28
			text3Loc <<- myXlim[2] + xWidth*0.28
			pvalLoc <<- myXlim[2] + xWidth*0.58
			xlim <- c( leftLabelLoc-xWidth*0.1, pvalLoc+xWidth*0.05)
			tableNames <- c( "Label", "Group1", "Group2", "Difference", "P_Value")
		} else {
			xWidth <- diff( myXlim) * 0.4
			leftLabelLoc <<- myXlim[1] - xWidth
			text1Loc <<- myXlim[1] - xWidth*0.35
			text2Loc <<- myXlim[1] - xWidth*0.15
			text3Loc <<- myXlim[2] + xWidth*0.33
			pvalLoc <<- myXlim[2] + xWidth*0.8
			xlim <- c( leftLabelLoc-xWidth*0.1, pvalLoc+xWidth*0.05)
			tableNames <- c( "Label", "N_Mean_SD", "Unused", "CI", "P_Value")
		}
		yMax <<- yBig
		ylim <- c( 0.25, yBig+1.5)
		text.cex <<- text.cex
		textTable <<- matrix( "", yMax,5)
		colnames(textTable) <<- tableNames 
		annotateText <<- annotateText

		# make the plot window
		plot( 1,1, type="n", main=main, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n",
			ylab=NA, xlab=NA, frame.plot=annotateText)
		axis( side=1, at=pretty( myXlim), cex.axis=text.cex*1.1)
		if ( ! is.null(sub)) title( sub=sub, line=2.25, cex=text.cex*1.1)

		lines( rep.int(noDifference,2), c(-1, yBig+0.5), col=1, lwd=1)
		if (annotateText) lines( xlim*2, rep.int(yBig+0.5,2), col=1, lwd=1)
		#lines( xlim*2, rep.int(0.5,2), col=1, lwd=1)

		# thin lines between the groups?
		if (dividerLines) {
			for (i in 1:yBig) {
				lines( xlim*3, rep.int(i-0.5,2), col='grey80', lwd=0.25)
			}	
		}
		# done with setup
		return()
	}


	meanDiff.header <- function( label1="Group 1", label2="Group 1", cex=1, prefix="Higher in",
					offset=1.2) {

		if (annotateText) text( text1Loc, yMax+offset, paste(label1,"N, Mean (SD)", sep="\n"), cex=text.cex*cex)
		textTable[ yMax, 2] <<- paste(label1,"N, Mean (SD)", sep="\n")
		if (annotateText) text( text2Loc, yMax+offset, paste(label2,"N, Mean (SD)", sep="\n"), cex=text.cex*cex)
		textTable[ yMax, 3] <<- paste(label2,"N, Mean (SD)", sep="\n")
		if (annotateText) text( text3Loc, yMax+offset, "Mean Difference\n(95% CI)", cex=text.cex*cex)
		textTable[ yMax, 4] <<- "Mean Difference\n(95% CI)"
		if (annotateText) text( pvalLoc, yMax+offset, " \nP-value", cex=text.cex*cex)
		textTable[ yMax, 5] <<- "P-value"

		# put UP in Grp 1 on the left side
		if (annotateText) text( xDataLim[1]*0.6, 0.15, paste( prefix, label1), cex=text.cex*cex)
		if (annotateText) text( xDataLim[2]*0.6, 0.15, paste( prefix, label2), cex=text.cex*cex)
		return()
	}


	mean.header <- function( label1="Group 1", label2="Group 2", cex=1, prefix="Higher in",
					offset=1.2, groupName="Cell Type") {

		text( leftLabelLoc, yMax+offset, groupName, cex=text.cex*cex)
		if (annotateText) text( text1Loc, yMax+offset, "N, Mean (SD)", cex=text.cex*cex)
		textTable[ yMax, 2] <<- "N, Mean (SD)"
		if (annotateText) text( text3Loc, yMax+offset, "(95% CI)", cex=text.cex*cex)
		textTable[ yMax, 4] <<- "(95% CI)"
		if (annotateText) text( pvalLoc, yMax+offset, "P-value", cex=text.cex)
		textTable[ yMax, 5] <<- "P-value"

		# put UP in Grp 1 on the right side
		if (annotateText) text( xDataLim[2]*0.6, 0, paste( prefix, " '", label1, "'", sep=""), cex=text.cex*cex)
		if (annotateText) text( xDataLim[1]*0.6, 0, paste( prefix, " '", label2, "'", sep=""), cex=text.cex*cex)
		return()
	}


	mean.line <- function( x1, yRow, label, digits=2, col=1, lwd=3, pch=23, cex=1, pt.cex=1.5) {
	
		if ( yRow < 1 || yRow > yMax) {
			stop( "Error:  Invalid 'yRow' -- outside bounds of forest plot.")
		}

		# if x1 is NULL, then just write the label as a heading
		if ( is.null(x1)) {
			x.offset <- if ( nchar(label) > 16) leftLabelLoc*((nchar(label)-16)/200) else 0
			if (annotateText) text( leftLabelLoc-x.offset, yRow, label, cex=text.cex*cex)
			textTable[ yRow, 1] <<- label
			return()
		}

		# given the data for one row / branch, calc all the terms & text
		m1 <- round( mean( x1, na.rm=T), digits=digits)
		sd1 <- round( sd( x1, na.rm=T), digits=digits)
		
		# do the mean calc
		ans <- t.test( x1)
		est <- ans$estimate
		pval <- ans$p.value
		ci <- round( ans$conf.int, digits=digits)
		meanAns <- round( est[1], digits=digits)

		# make the labels and write them
		if (annotateText) text( leftLabelLoc, yRow, label, cex=text.cex*cex)
		textTable[ yRow, 1] <<- label
		myText1 <- paste( length(x1), ", ", m1, " (", sd1, ")", sep="")
		if (annotateText) text( text1Loc, yRow, myText1, cex=text.cex*cex)
		textTable[ yRow, 2] <<- myText1
		myText3 <- paste( "(", ci[1], ", ", ci[2], ")", sep="")
		x.offset <- if ( nchar(myText3) > 16) leftLabelLoc*((nchar(myText3)-16)/200) else 0
		if (annotateText) text( text3Loc+x.offset, yRow, myText3, cex=text.cex*cex)
		textTable[ yRow, 4] <<- myText3
		pvalText <- if ( pval >= 0.0001) formatC( pval, format="f", digits=4) else "< 0.0001"
		if (annotateText) text( pvalLoc, yRow, pvalText, cex=text.cex*cex)
		textTable[ yRow, 5] <<- pvalText

		# now draw the line and marks
		smallX <- diff( range( xDataLim)) * 0.02
		xx1 <- ci[1]
		xx2 <- ci[2]
		clipLeft <- clipRight <- FALSE
		if ( xx1 < xDataLim[1]) {
			clipLeft <- TRUE
			xx1 <- xDataLim[1]
			xx2 <- max( xx1+smallX, xx2)
		}
		if ( xx2 > xDataLim[2]) {
			clipRight <- TRUE
			xx2 <- xDataLim[2]
			xx1 <- min( xx1, xx2-smallX)
		}
		lines( c(xx1,xx2), c(yRow,yRow), lwd=lwd, col=col)
		# show arror head to denote more data
		if ( any( c( clipLeft,clipRight))) {
			smallY <- 0.1
			if ( clipLeft) {
				lines( c(xx1,xx1+smallX), c(yRow,yRow+smallY), lwd=lwd, col=col)
				lines( c(xx1,xx1+smallX), c(yRow,yRow-smallY), lwd=lwd, col=col)
			}
			if ( clipRight) {
				lines( c(xx2,xx2-smallX), c(yRow,yRow+smallY), lwd=lwd, col=col)
				lines( c(xx2,xx2-smallX), c(yRow,yRow-smallY), lwd=lwd, col=col)
			}
		}
		# always show the mean mark.  Clip to limits if needed
		if ( meanAns < xDataLim[1]) meanAns <- xx1
		if ( meanAns > xDataLim[2]) meanAns <- xx2
		points( meanAns, yRow, pch=pch, col=col, bg=col, cex=pt.cex)

		# done
		return()
	}


	meanDiff.line <- function( x1, x2, yRow, label, digits=2, col=1, lwd=3, pch=23, cex=1, pt.cex=1.8) {
	
		if ( yRow < 1 || yRow > yMax) {
			stop( "Error:  Invalid 'yRow' -- outside bounds of forest plot.")
		}

		# if x1 is NULL, then just write the label as a heading
		if ( is.null(x1)) {
			x.offset <- if ( nchar(label) > 16) leftLabelLoc*((nchar(label)-16)/200) else 0
			if (annotateText) text( leftLabelLoc-x.offset, yRow, label, cex=text.cex*cex)
			textTable[ yRow, 1] <<- label
			return()
		}

		# given the data for one row / branch, calc all the terms & text
		m1 <- round( mean( x1, na.rm=T), digits=digits)
		sd1 <- round( sd( x1, na.rm=T), digits=digits)
		m2 <- round( mean( x2, na.rm=T), digits=digits)
		sd2 <- round( sd( x2, na.rm=T), digits=digits)
		
		# do the mean diff calc
		# changed to put Up in Grp1 on the left
		ans <- t.test( x2, x1)
		est <- ans$estimate
		pval <- ans$p.value
		ci <- round( ans$conf.int, digits=digits)
		meanDiff <- round( est[1] - est[2], digits=digits)

		# make the labels and write them
		if (annotateText) text( leftLabelLoc, yRow, label, cex=text.cex*cex)
		textTable[ yRow, 1] <<- label
		myText1 <- paste( length(x1), ", ", m1, " (", sd1, ")", sep="")
		if (annotateText) text( text1Loc, yRow, myText1, cex=text.cex*cex)
		textTable[ yRow, 2] <<- myText1
		myText2 <- paste( length(x2), ", ", m2, " (", sd2, ")", sep="")
		if (annotateText) text( text2Loc, yRow, myText2, cex=text.cex*cex)
		textTable[ yRow, 3] <<- myText2
		myText3 <- paste( meanDiff, " (", ci[1], ",", ci[2], ")", sep="")
		x.offset <- if ( nchar(myText3) > 16) leftLabelLoc*((nchar(myText3)-16)/200) else 0
		if (annotateText) text( text3Loc+x.offset, yRow, myText3, cex=text.cex*cex)
		textTable[ yRow, 4] <<- myText3
		pvalText <- if ( pval >= 0.0001) formatC( pval, format="f", digits=4) else "< 0.0001"
		if (annotateText) text( pvalLoc, yRow, pvalText, cex=text.cex)
		textTable[ yRow, 5] <<- pvalText

		# now draw the line and marks
		xx1 <- ci[1]
		xx2 <- ci[2]
		clipLeft <- clipRight <- FALSE
		if ( xx1 < xDataLim[1]) {
			clipLeft <- TRUE
			xx1 <- xDataLim[1]
		}
		if ( xx2 > xDataLim[2]) {
			clipRight <- TRUE
			xx2 <- xDataLim[2]
		}
		lines( c(xx1,xx2), c(yRow,yRow), lwd=lwd, col=col)
		# show arror head to denote more data
		if ( any( c( clipLeft,clipRight))) {
			smallX <- diff( range( xDataLim)) * 0.02
			smallY <- 0.1
			if ( clipLeft) {
				lines( c(xx1,xx1+smallX), c(yRow,yRow+smallY), lwd=lwd, col=col)
				lines( c(xx1,xx1+smallX), c(yRow,yRow-smallY), lwd=lwd, col=col)
			}
			if ( clipRight) {
				lines( c(xx2,xx2-smallX), c(yRow,yRow+smallY), lwd=lwd, col=col)
				lines( c(xx2,xx2-smallX), c(yRow,yRow-smallY), lwd=lwd, col=col)
			}
		}
		if ( meanDiff >= xx1 && meanDiff <= xx2) points( meanDiff, yRow, pch=pch, col=col, bg=col, cex=pt.cex)

		# done
		return()
	}

	result.text <- function() {

		# pass back the text we collected, invert the rows so it matches the image
		nr <-nrow( textTable)
		out <- textTable[ nr:1, ]
		rownames(out) <- NULL
		out
	}

	# all functions defined.  Return this environment
	return( environment())
}
