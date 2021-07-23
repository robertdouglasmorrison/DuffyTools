# timeHourTools.R  -  various ways to make use of the DeRisi time hour datasets


# reload the one "TimeHour Data Sets" object
`loadTimeHourData` <- function() {

	if ( ! exists( "comboSet", envir=TimeHourEnv)) {
		cat( "\nLoading TimeHour data...")
		data( TimeHourDatasets, envir=TimeHourEnv)
		cat( "\n")
	}
	return()
}


# grab the GeneIDs known by the time hour datasets
`getTimeHourGenes` <- function() {

	loadTimeHourData()
	comboSet <- get( "comboSet", envir=TimeHourEnv)
	return( comboSet$geneID)
}


`plotTimeHourCurveOneGene` <- function( gene) {

	loadTimeHourData()
	threeD7set <- get( "threeD7set", envir=TimeHourEnv)
	HB3set <- get( "HB3set", envir=TimeHourEnv)
	comboSet <- get( "comboSet", envir=TimeHourEnv)

	threeD7smooth <- threeD7set$smoothSet
	HB3smooth <- HB3set$smoothSet
	threeD7raw <- threeD7set$rawSet
	HB3raw <- HB3set$rawSet
	combosmooth <- comboSet$smoothSet
	comboraw <- comboSet$rawSet

	where3d7 <- base::match( gene, rownames( threeD7smooth), nomatch=0)
	wherehb3 <- base::match( gene, rownames( HB3smooth), nomatch=0)
	whereCombo <- base::match( gene, rownames( combosmooth), nomatch=0)
	rank3d7 <- base::match( gene, threeD7set$geneID, nomatch=NA)
	rankhb3 <- base::match( gene, HB3set$geneID, nomatch=NA)
	rankcombo <- base::match( gene, comboSet$geneID, nomatch=NA)

	N <- ncol( threeD7smooth)

	p3D7 <- where3d7[1]
	pHB3 <- wherehb3[1]
	pCombo <- whereCombo[1]
	if ( ! any( c(p3D7, pHB3, pCombo) > 0)) {
		cat( "\nGene not in Time Hour dataset: ", gene[1])
		return()
	}
	yLim <- range( c( threeD7raw[ p3D7, ], HB3raw[ pHB3, ]))

	showOtto <- showStun <- FALSE
	pOtto <- pStun <- 0
	# allow drawing the Otto curve too
	ottoSet <- get( "ottoSet", envir=TimeHourEnv)
	if ( ! is.null(ottoSet)) {
		showOtto <- TRUE
		ottoData <- ottoSet$rawSet
		whereOtto <- base::match( gene, rownames( ottoData), nomatch=0)
		pOtto <- whereOtto[1]
		if ( pOtto > 0) {
			ottoX <- ottoSet$Hours
			ottoY <- ottoData[ pOtto, ]
			ottoY <- (ottoY * ottoSet$Slope) + ottoSet$Intercept
			yLim <- range( c( yLim, ottoY))
		}
	}
	stunSet <- get( "stunSet", envir=TimeHourEnv)
	if ( ! is.null(stunSet)) {
		showStun <- TRUE
		stunData <- stunSet$rawSet
		whereStun <- base::match( gene, rownames( stunData), nomatch=0)
		pStun <- whereStun[1]
		if ( pStun > 0) {
			stunX <- stunSet$Hours
			stunY <- stunData[ pStun, ]
			stunY <- (stunY * stunSet$Slope) + stunSet$Intercept
			yLim <- range( c( yLim, stunY))
		}
	}

	if ( yLim[2] - yLim[1] < 2000) {
		yLim[1] <- yLim[1] * 0.5
	}
	yLim[2] <- yLim[2] * 2

	if ( ! (getCurrentSpecies() == "Pf3D7")) setCurrentSpecies( "Pf3D7")

	plot( 10,10, main=paste("Gene Expression Time Hour Profile:    ", gene, "\n", gene2Product(gene)), 
		xlab="Intraerythrocytic Cell Cycle Hour", ylab="Expression Level", 
		xlim=c(1,N+4), type="n", ylim=yLim, log="y", font.axis=2, font.lab=2)
	lines( x=1:N, y=threeD7smooth[ p3D7, ], type="l", col=2, lwd=3)
	lines( x=1:N, y=HB3smooth[ pHB3, ], type="l", col=3, lwd=3)
	lines( x=1:N, y=combosmooth[ pCombo, ], type="l", col=1, lwd=3)
	lines( x=1:N, y=threeD7raw[ p3D7, ], type="l", col=2, lwd=2, lty=3)
	lines( x=1:N, y=HB3raw[ pHB3, ], type="l", col=3, lwd=2, lty=3)
	lines( x=1:N, y=comboraw[ pCombo, ], type="l", col=1, lwd=2, lty=3)

	if ( pOtto > 0) {
		lines( x=ottoX, y=ottoY, type="l", col=4, lwd=2, lty=3)
		logOttoY <- log2( ottoY+1)
		ans <- smooth.spline( x=ottoX, y=logOttoY, keep.data=TRUE, 
				all.knots=TRUE, df=round(length(ottoY)*0.67),
				control.spar=list( "low"=0.0, "high"=1.0))
		logOttoY <- ans$y
		ottoY <- 2 ^ logOttoY - 1
		lines( x=ottoX, y=ottoY, type="l", col=4, lwd=3, lty=1)
	}
	if ( pStun > 0) {
		lines( x=stunX, y=stunY, type="l", col=5, lwd=2, lty=3)
		logStunY <- log2( stunY+1)
		ans <- smooth.spline( x=stunX, y=logStunY, keep.data=TRUE, 
				all.knots=TRUE, df=round(length(stunY)*0.67),
				control.spar=list( "low"=0.0, "high"=1.0))
		logStunY <- ans$y
		stunY <- 2 ^ logStunY - 1
		lines( x=stunX, y=stunY, type="l", col=5, lwd=3, lty=1)
	}

	legendText <- paste( c("3D7", "HB3", "Combo"), "  (rank=", 
			c(rank3d7[1],rankhb3[1],rankcombo[1]),")",sep="") 
	legendText <- c( legendText, "Otto RNA-seq", "Stun RNA-seq")
	legend( x="topright", legend=legendText, lwd=3, col=c(2,3,1,4,5), bg="white")
	return()
}


`getTimeHourConsensusBestHour` <- function( df, minimumIntensity=1) {

	loadTimeHourData()
	HB3set <- get( "HB3set", envir=TimeHourEnv)
	comboSet <- get( "comboSet", envir=TimeHourEnv)

	nSet <- c( 100, 250, 500, 1000)
	scoreSet <- c( "pearson") 
	TSset <- list( comboSet, HB3set)

	SCALE="max"
	finalScore <- rep( 0, times=ncol( TSset[[1]]$smoothSet))
	nScores <- 0

	for ( i in 1:4) {
	for ( j in 1:2) {
	for ( k in 1:1) {
		N <- nSet[i]
		ts <- TSset[[j]]
		SCORE <- scoreSet[k]
		ans <- bestTimeHour( df, ts, scaleMethod=SCALE, scoreMethod=SCORE, nTSgenes=N,
				minimumIntensity=minimumIntensity)
		myscore <- ans$Scores
		# various errors could return NA or NAN
		if ( all( is.na( myscore))) next
		if ( all( is.nan( myscore))) next
		finalScore <- finalScore + myscore
		nScores <- nScores + 1
	}}}

	finalScore <- finalScore / nScores
	bestHour <- which.max( finalScore)
	return( list( "BestTimePoint"=bestHour, "Scores"=finalScore))
}


`surveyBestTimeHour` <- function( df, label="", minimumIntensity=1) {

	loadTimeHourData()
	HB3set <- get( "HB3set", envir=TimeHourEnv)
	comboSet <- get( "comboSet", envir=TimeHourEnv)

	nSet <- c( 100, 250, 500, 1000)
	scoreSet <- c( "pearson", "spearman")
	TSset <- list( "Combo"=comboSet, "HB3"=HB3set)
	LOGset <- c( TRUE, FALSE)
	MAX_HRS <- ncol( comboSet$smoothSet)

	plot( 1,1, type="n", xlim=c(0,MAX_HRS+4), ylim=c( -0.3,0.95), main=paste( "Timing Evaluation:   ", label), 
		xlab="Time Hour", ylab="Score")

	finalScore <- rep( 0, times=ncol( TSset[[1]]$smoothSet))
	nScores <- 0

	for ( i in 1:4) {
	for ( j in 1:2) {
	for ( k in 1:2) {
	for ( l in 1:2) {
		N <- nSet[i]
		SCALE <- "max"
		SCORE <- scoreSet[k]
		LOG <- LOGset[l]
		TS <- names(TSset)[j]
		ans <- bestTimeHour( df, TSset[[j]], scaleMethod=SCALE, scoreMethod=SCORE, useLog=LOG, 
					nTSgenes=N, minimumIntensity=minimumIntensity)
		myscore <- ans$Scores
		lines( myscore, col=i+1, lty=(k), lwd=(j), type=c("l","p")[l], cex=0.85)
		myx <- which.max(myscore)
		myy <- max(myscore)
		text( myx,myy, base::paste( TS, SCORE, LOG, N, sep=":"), pos=3, cex=0.75, col=i+1)

		finalScore <- finalScore + myscore
		nScores <- nScores + 1
	}}}}

	finalScore <- finalScore / nScores
	lines( finalScore, col=1, lty=1, lwd=3)

	legend( "topright", c("Colors show N_genes", "Dashes  show ScoreMetric", "Widths show TimeStrain", "Points  show LogMode off"))

	bestHour <- which.max( finalScore)
	return( list( "BestTimePoint"=bestHour, "Scores"=finalScore))
}


# find the best hour for a given set of intensities against one training set.
`bestTimeHour` <- function( df, tSet, scaleMethod="max", scoreMethod="pearson", useLog=TRUE, nTSgenes=NULL,
				minimumIntensity=1) {

	# find the columns we want
	gCol <- intenCol <- 0
	for ( icol in 1:ncol(df)) {
		if ( base::match( colnames(df)[icol], c("Gene", "GeneID", "GENE_ID", "EXON_ID"),
						nomatch=0) > 0) gCol <- icol
		if ( base::match( colnames(df)[icol], c("Intensity", "Inten", "INTENSITY", "RPKM_M", "AAPKM"),
						nomatch=0) > 0) intenCol <- icol
	}

	# make a small set of just the wanted genes
	gWanted <- tSet$geneIDs
	gPointers <- tSet$genePtrs
	if ( ! is.null( nTSgenes)) {
		gWanted <- tSet$geneID[ 1:nTSgenes]
		gPointers <- tSet$genePtrs[ 1:nTSgenes]
	}

	# trap and exclude any non-genes from the dataset
	drops <- grep( "(ng)", df[ , gCol], fixed=T)
	if ( length( drops) > 0) {
		df <- df[ -drops, ]
	}

	where <- base::match( gWanted, df[ , gCol], nomatch=0)
	whoFound <- ( where > 0)
	if ( ! all( whoFound)) {
		#warning( paste( "bestTimeHour:  missing some TimePoint marker genes: ", sum( !whoFound)))
		where <- where[ whoFound]
		gPointers <- gPointers[ whoFound]
	}
	myGenes <- df[ where, gCol]
	myIntens <- df[ where, intenCol]
	if ( any( myIntens < minimumIntensity)) myIntens <- myIntens + minimumIntensity
	myDF <- data.frame( myGenes, myIntens, gPointers, stringsAsFactors=FALSE)
	colnames( myDF) <- c("Gene", "Intensity", "RowPtr")

	# do the measurement at each time
	scores <- rep( 0, times=ncol(tSet$smoothSet))
	for ( itime in 1:ncol( tSet$smoothSet)) {
		ans <- calcTimeScore( myDF, tSet$smoothSet[ ,itime], itime, scaleMethod=scaleMethod, 
					scoreMethod=scoreMethod, useLog=useLog)
		scores[itime] <- ans$score
	}

	# pearsons R is already -1 to 1, lets turn RMS deviation into the same idea
	if ( scoreMethod == "rms") {
		bigScore <- max( scores)
		scores <- (bigScore - scores) / bigScore
	}
	ord <- base::order( scores, decreasing=TRUE)

	return( list( "BestTimePoints"=ord[1:5], "Scores"=scores))
}


# score one set of marker genes against one training set time point intensity vector
`calcTimeScore` <- function( df, tVec, curtime, scaleMethod="max", scoreMethod="pearson", useLog=TRUE) {

	# get those marker gene intensities from this traning set hour vector
	nG <- nrow(df)
	tPtr <- df[ ,3]
	tSetInten <- tVec[ tPtr]

	# get the matching intensities from our sample
	mySetInten <- df[ , 2]

	# set a scale factor to put the two sets on "equal" footing
	if ( scaleMethod == "sum" ) {
		scaleFactor <- sum( tSetInten, na.rm=TRUE) / sum( mySetInten, na.rm=TRUE)
		rmsScale <- mean.default( tSetInten)
	} else if ( scaleMethod == "mean" ) {
		scaleFactor <- mean.default( tSetInten) / mean.default( mySetInten, na.rm=TRUE)
		rmsScale <- mean.default( tSetInten)
	} else if ( scaleMethod == "median" ) {
		scaleFactor <- median( tSetInten) / median( mySetInten, na.rm=TRUE)
		rmsScale <- median( tSetInten)
	} else if ( scaleMethod == "max" ) {
		scaleFactor <- max( tSetInten) / max( mySetInten, na.rm=TRUE)
		rmsScale <- max( tSetInten)
	} else if ( scaleMethod == "sum*max" ) {
		scaleFactor1 <- max( tSetInten) / max( mySetInten, na.rm=TRUE)
		scaleFactor2 <- sum( tSetInten) / sum( mySetInten, na.rm=TRUE)
		scaleFactor <- sqrt( scaleFactor1 * scaleFactor2)
		rmsScale <- mean.default( c( max( tSetInten), mean.default( tSetInten)))
	} else if ( scaleMethod == "" ) {
		scaleFactor <- 1.0
		rmsScale <- 1.0
	} else {
		stop( paste("calcTimeScore:  unknown scaling method: ", scaleMethod))
	}
	mySetInten <- mySetInten * scaleFactor

	# gracefully catch no/little intensity
	if ( diff( range( mySetInten, na.rm=T)) < 1) return( list( score=0))

	# make a score that measures similarity
	if ( scoreMethod == "rms") {
		diff <- mySetInten - tSetInten
		rms <- sqrt( mean.default( diff * diff, na.rm=TRUE))
		# convert the rms deviation to a relative measure
		out <- rms / rmsScale
	} else if ( scoreMethod == "pearson") {
		if ( useLog) {
			ans <- cor( log2(mySetInten+1), log2(tSetInten+1), use="complete", method="pearson")
		} else {
			ans <- cor( mySetInten, tSetInten, use="complete", method="pearson")
		}
		# R is already -1 to 1 range...
		out <- ans
	} else if ( scoreMethod == "spearman") {
		ans <- cor( mySetInten, tSetInten, use="complete", method="spearman")
		# R is already -1 to 1 range...
		out <- ans
	} else {
		stop( paste( "calcTimeScore:  unknown score method: ", scoreMethod))
	}

	return( list( "score"=out))
}



`plotTimeHourCurvesFromFileSet` <- function( fnames, fids, fcolors=1, fgrps="group", fLwd=2, 
				fLty=1, xMax=48, label="label goes here", legend.cex=1.0,
				minimumIntensity=1, sep="\t", geneUniverse=NULL) {

	mainText <- paste( "Time Hour Estimates: \n", label)
	plot( 1,1, type="n", xlim=c(0,xMax), ylim=c(-0.2, 1.0), main=mainText, xlab="DeRisi Time Hour", 
			ylab="Score  (combined R value)", font.lab=2, font.axis=2)
	nF <- length( fnames)
	if ( length(fcolors) < nF) fcolors <- rep( fcolors, length.out=nF)
	if ( length(fgrps) < nF) fgrps <- rep( fgrps, length.out=nF)
	if ( length(fLty) < nF) fLty <- rep( fLty, length.out=nF)
	if ( length(fLwd) < nF) fLwd <- rep( fLwd, length.out=nF)

	timeOut <- scoreOut <- vector()
	labX <- labY <-vector()
	for( i in 1:nF) {
		f <- fnames[i]; 
		if ( !file.exists(f)) {
			cat( "\nSkipping for 'file not found': ", f)
			next;
		}
		mydf <- read.delim( f, as.is=T, sep=sep); 

		# allow using a subset of genes
		if ( ! is.null( geneUniverse)) {
			geneCol <- which( colnames(mydf) %in% c("GENE_ID","GeneID","Gene"))[1]
			if ( is.na( geneCol)) stop( "plotTimeHourCurve: Unable to find gene ID column in transcriptone")
			dfGenes <- mydf[[ geneCol]]
			keep <- which( dfGenes %in% geneUniverse)
			mydf <- mydf[ keep, ]
			if ( nrow(mydf) < 100) {
				cat( "\nplotTimeHourCurve: not enough genes remain after 'geneUniverse' filtering")
				next
			}
		}

		ans <- getTimeHourConsensusBestHour( mydf, minimumIntensity=minimumIntensity); 
		if ( is.null(ans)) next
		lines( ans$Scores, col=fcolors[i], lwd=fLwd[i], lty=fLty[i]); 
		text( which.max( ans$Scores), max( ans$Scores), fids[i], col=fcolors[i], font=2, pos=3)
		timeOut[i] <- ans$BestTimePoint
		labX[i] <- which.max( ans$Scores)
		labY[i] <- scoreOut[i] <- max( ans$Scores)
		cat( "\n", fids[i], "\tTime=", timeOut[i], "\tScore=", scoreOut[i])
	}
	text( labX, labY, fids, col=fcolors, font=2, pos=3)

	myGrps <- unique.default( fgrps)
	where <- base::match( myGrps, fgrps)
	myCols <- fcolors[where]
	if ( any( myGrps != "")) legend( "bottomright", myGrps, lwd=4, col=myCols, bg="white", cex=legend.cex)

	return( data.frame( "SampleID"=fids, "BestTimeHour"=timeOut, "BestScore"=scoreOut))
}


`plotTimeHourCurvesFromMatrix` <- function( m, fcolors=1, fgrps="group", fLwd=2, 
				fLty=1, xMax=48, label="label goes here", legend.cex=1.0,
				minimumIntensity=1, geneUniverse=NULL) {

	mainText <- paste( "Time Hour Estimates: \n", label)
	plot( 1,1, type="n", xlim=c(0,xMax), ylim=c(-0.2, 1.0), main=mainText, xlab="DeRisi Time Hour", 
			ylab="Score  (combined R value)", font.lab=2, font.axis=2)
	nF <- ncol(m)
	if ( length(fcolors) < nF) fcolors <- rep( fcolors, length.out=nF)
	if ( length(fgrps) < nF) fgrps <- rep( fgrps, length.out=nF)
	if ( length(fLty) < nF) fLty <- rep( fLty, length.out=nF)
	if ( length(fLwd) < nF) fLwd <- rep( fLwd, length.out=nF)

	timeOut <- scoreOut <- vector()
	labX <- labY <-vector()

	# allow using a subset of genes
	if ( ! is.null( geneUniverse)) {
		mGenes <- rownames(m)
		keep <- which( mGenes %in% geneUniverse)
		m <- m[ keep, ]
		if ( nrow(m) < 100) {
			cat( "\nplotTimeHourCurve: not enough genes remain after 'geneUniverse' filtering")
			return(NULL)
		}
	}

	genes <- rownames(m)
	fids <- colnames(m)

	for( i in 1:nF) {
		mydf <- data.frame( "GENE_ID"=genes, "INTENSITY"=m[ , i], stringsAsFactors=F)
		ans <- getTimeHourConsensusBestHour( mydf, minimumIntensity=minimumIntensity); 
		lines( ans$Scores, col=fcolors[i], lwd=fLwd[i], lty=fLty[i]); 
		text( which.max( ans$Scores), max( ans$Scores), fids[i], col=fcolors[i], font=2, pos=3)
		timeOut[i] <- ans$BestTimePoint
		labX[i] <- which.max( ans$Scores)
		labY[i] <- scoreOut[i] <- max( ans$Scores)
		cat( "\n", fids[i], "\tTime=", timeOut[i], "\tScore=", scoreOut[i])
	}
	text( labX, labY, fids, col=fcolors, font=2, pos=3)

	myGrps <- unique.default( fgrps)
	where <- base::match( myGrps, fgrps)
	myCols <- fcolors[where]
	if ( any( myGrps != "")) legend( "bottomright", myGrps, lwd=4, col=myCols, bg="white", cex=legend.cex)

	return( data.frame( "SampleID"=fids, "BestTimeHour"=timeOut, "BestScore"=scoreOut))
}


`getTimeHourGeneFCDTT` <- function( gset, fromHour=10, toHour=16) {

	# get the change in expected Gene expression dut to time hour alone...
	loadTimeHourData()
	comboSet <- get( "comboSet", envir=TimeHourEnv)

	smoothset <- comboSet$smoothSet
	TSgenes <- rownames( smoothset)

	fcOut <- rep( 1, times=length( gset))
	hourRateOut <- rep( 1, times=length( gset))

	where <- base::match( gset, TSgenes, nomatch=0)
	nHours <- abs( toHour - fromHour)
	if ( nHours == 0) nHours <- 1

	for ( i in 1:length( gset)) {
		if ( where[i] == 0) next
		int1 <- smoothset[ where[i], fromHour]
		int2 <- smoothset[ where[i], toHour]
		fc <- int2 / int1
		rate <-  2 ^ ( log(fc,2) / nHours)
		fcOut[i] <- fc
		hourRateOut[i] <- rate
	}

	names( fcOut) <- names( hourRateOut) <- gset
	return( list( "FCDTT"=fcOut, "FCperHr"=hourRateOut))
}


`getTimeHourGeneMFCDTT` <- function( gset, fromHour=1, toHour=30) {

	# get the Maximum Fold change Due To Time hour alone, over a given time window...
	loadTimeHourData()
	comboSet <- get( "comboSet", envir=TimeHourEnv)

	smoothset <- comboSet$smoothSet
	TSgenes <- rownames( smoothset)

	mfcOut <- rep( NA, times=length( gset))
	hourRateOut <- rep( NA, times=length( gset))

	where <- base::match( gset, TSgenes, nomatch=0)
	hourRange <- (max(1,fromHour) : min( toHour, ncol(smoothset)))

	# calc the actual rate for those in the time hour dataset
	for ( i in 1:length( gset)) {
		if ( where[i] == 0) next
		myValues <- smoothset[ where[i], hourRange]
		hr1 <- which.min( myValues) + fromHour - 1
		hr2 <- which.max( myValues) + fromHour - 1
		int1 <- smoothset[ where[i], hr1]
		int2 <- smoothset[ where[i], hr2]
		fc <- int2 / int1
		nHours <- abs( hr2 - hr1)
		if (nHours == 0) nHours <- 1
		rate <-  2 ^ ( log(fc,2) / nHours)
		mfcOut[i] <- fc
		hourRateOut[i] <- rate
	}

	# take an average to give to all the 'not found' genes
	avgFC <- mean.default( mfcOut, na.rm=TRUE)
	if ( ! is.finite( avgFC)) avgFC <- 2
	avgRate <- mean.default( hourRateOut, na.rm=TRUE)
	if ( ! is.finite( avgRate)) avgRate <- avgFC / max( c( 1, (toHour-fromHour)))
	mfcOut[ is.na( mfcOut)] <- avgFC
	hourRateOut[ is.na( hourRateOut)] <- avgRate

	names( mfcOut) <- names( hourRateOut) <- gset
	return( list( "MFCDTT"=mfcOut, "MFCperHr"=hourRateOut))
}


`calcExcessFoldOverMFCDTT` <- function( gset, log2fold, fromHour=1, toHour=30) {

	# get the Excess Fold change,  i.e. how much the reported Fold exceeds the 
	#Maximum Fold change Due To Time hour alone, over a given time window...

	ans <- getTimeHourGeneMFCDTT( gset, fromHour, toHour)

	# express all folds as positives for the calculation of 'excess'
	logExcess <- abs( log2fold) - log( ans$MFCDTT, 2)
	excessFold <- 2 ^ logExcess
	pvalues <- rep( 1, times=length( gset))

	for( i in 1:length( gset)) {

		# lets assume that 95% of the time, the observed FC is due to Time Hour alone,
		# so we expect the observed FC to be less than 2 SD from mean
		mysd <- ans$MFCDTT[i]
		if ( mysd < 1) mysd <- 1
		# express the observed FC as a positive relative to "no fold change = 0"
		thisFold <- (2 ^ abs(log2fold[i]))
		# trap any 'bad values'...
		if ( is.na(thisFold) || !is.finite(thisFold)) thisFold <- 1

		pval <- pnorm( thisFold, mean=0, sd=mysd, lower.tail=F)
		# since it could be on either side of the expected range, double what we report...
		pvalues[i] <- pval * 2
	}

	# lastly, label all the excess folds from the "downn' genes as being 'down'
	isNeg <- which( log2fold < 0)
	excessFold[ isNeg] <- -( excessFold[ isNeg])

	names( excessFold) <- names( pvalues) <- gset

	return( list( "MFCDTT"=ans$MFCDTT, "ExcessFold"=excessFold, "PvalueExcessFold"=pvalues))
}


`plotMFCDTTdistribution` <- function( gset="", fromHour=1, toHour=30, foldCutoff=2.0, asLog=TRUE, 
		label="your label...", sort=FALSE, col=2, cex=1.0, geneUniverse=NULL) {

	setCurrentSpecies( "Pf3D7")

	allG <- subset( getCurrentGeneMap(), REAL_G == TRUE)$GENE_ID
	if ( ! is.null( geneUniverse)) {
		allG <- intersect( allG, geneUniverse)
		cat( "\nusing a Gene Universe of: ", length( allG))
	} 
	mfc <- getTimeHourGeneMFCDTT( allG, fromHour=fromHour, toHour=toHour)
	mfc <- as.data.frame( mfc)
	Nmfc <- length( allG)
	Ngenes <- length( gset)
	yRange <- range( mfc$MFCDTT)
	
	if ( sort) {
		ord <- base::order( mfc$MFCDTT)
		mfc <- mfc[ ord, ]
	} else {
		if ( asLog) {
			yRange[2] <- yRange[2] * 1.6
		} else {
			yRange[2] <- yRange[2] * 1.3
		}
	}

	mainText <- paste( "Maximum Fold Change Due To Time:   (MFCDTT)\n", label)
	subText <- paste( "Time Window:    Hour ", fromHour, "  to Hour ", toHour)
	logTerm <- if (asLog) "y" else ""

	myX <- 1:Nmfc
	myY <- mfc$MFCDTT 
	plot( x=myX, y=myY, log=logTerm, ylim=yRange, main=mainText, sub=subText, xlab="Genes", 
		ylab="Max Fold Change in Time Window")

	lines( c( -Nmfc, Nmfc*2), c( foldCutoff, foldCutoff), lty=1, lwd=3, col="gray48")

	who <- base::match( gset, rownames( mfc)); 
	who <- who[ !is.na(who)]
	if ( length( who) > 0) {
		points( myX[who], myY[who], col=col, pch=19, cex=cex); 
		whoOut <- intersect( which( mfc$MFCDTT > foldCutoff), who)
		genesOut <- rownames(mfc)[ whoOut]
		text( myX[whoOut], myY[whoOut], genesOut, pos=if(sort) 2 else 3, col=col, cex=0.85, font=2 )
	} else {
		genesOut <- vector()
	}

	nBig <- sum( mfc$MFCDTT > foldCutoff); 
	nSml <- sum( mfc$MFCDTT[who] > foldCutoff, na.rm=T); 
	M <- matrix( c( nSml, Ngenes-nSml, nBig, Nmfc-nBig), nrow=2, ncol=2); 
	ans <- chisq.test(M); 
	smlPct <- as.percent( nSml, big.value=Ngenes);
	bigPct <- as.percent( nBig, big.value=Nmfc);

	if ( length( who) > 0) {
	   legend( "topleft", c( 
		paste( "Given set ( N above cutoff =", nSml, ")   ", smlPct), 
		paste( "All genes ( N above cutoff =", nBig, ")   ", bigPct), 
		paste( "Chi-squared test P value =", formatC( ans$p.value, format="E", digits=2))),
		pch=c( 19,21,NA), col=c(col,1,NA), bg="white")
	} else {
	   legend( "topleft", c( 
		paste( "All genes ( N above cutoff =", nBig, ")   ", bigPct)), 
		pch=c( 21), col=c(1), bg="white")
	}


	out <- list( "Genes"=genesOut, "Pvalue"=ans$p.value)
	return( out)
}
