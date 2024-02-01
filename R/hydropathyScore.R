# hydropathyScore.R -- calculate hydrophobic / hydrophilic score of proteins

hydropathyScore <- function( pepSet, window=9, phobic.cutoff=10, philic.cutoff=-10, 
				makePlot=F, col=1, lwd=3) {

	AA_Scores <- c( "A"=1.8, "C"=2.5, "D"=-3.5, "E"=-3.5, "F"=2.8, "G"=-0.4, "H"=-3.2, "I"=4.5,
			"K"=-3.9, "L"=3.8, "M"=1.9, "N"=-3.5, "P"=-1.6, "Q"=-3.5, "R"=-4.5, "S"=-0.8,
			"T"=-0.7, "V"=4.2, "W"=-0.9, "Y"=-1.3)

	aaSet <- strsplit( pepSet, split="")
	N <- length(pepSet)
	col <- rep( col, length.out=N)
	lwd <- rep( lwd, length.out=N)
	out <- rep.int( NA, N)
	tmp1 <- tmp2 <- outDF <- vector( mode="list", length=N)
	bigNAA <- 1
	for ( k in 1:N) {
		aa <- aaSet[[k]]
		NAA <- length(aa)
		if (NAA > bigNAA) bigNAA <- NAA
		wh <- match( aa, names(AA_Scores))
		v1 <- as.numeric( AA_Scores[wh])
		v1[ is.na(v1)] <- 0
		tmp1[[k]] <- v1
		out[k] <- sum(v1)
		v2 <- movingAverage( v1, window=window, FUN=sum)
		tmp2[[k]] <- v2
		# find all regions that stay above the cutoffs, as regions of interest
		# turn the data into T/F that we can sum
		hitType <- hitStart <- hitStop <- hitSeq <- vector()
		isUP <- ifelse( v2 >= phobic.cutoff, 1, -NAA)
		for ( i in 1:NAA) {
			if ( isUP[i] < 0) next
			cs <- cumsum( isUP[i:NAA])
			pk <- which.max(cs)[1]
			if (is.na(pk)) next
			j <- (i+pk-1)
			myRange <- i:j
			if ( length(myRange) >= window) {
				# got a domain
				hitType <- c( hitType, "hydrophobic")
				hitStart <- c( hitStart, i)
				hitStop <- c( hitStop, j)
				hitSeq <- c( hitSeq, paste( aa[i:j], collapse=""))
				# skip ahead to this end
				isUP[ i:j] <- -1
			}
		}
		isDOWN <- ifelse( v2 <= philic.cutoff, 1, -NAA)
		for ( i in 1:NAA) {
			if ( isDOWN[i] < 0) next
			cs <- cumsum( isDOWN[i:NAA])
			pk <- which.max(cs)[1]
			if (is.na(pk)) next
			j <- (i+pk-1)
			myRange <- i:j
			if ( length(myRange) >= window) {
				# got a domain
				hitType <- c( hitType, "hydrophilic")
				hitStart <- c( hitStart, i)
				hitStop <- c( hitStop, j)
				hitSeq <- c( hitSeq, paste( aa[i:j], collapse=""))
				# skip ahead to this end
				isDOWN[ i:j] <- -1
			}
		}
		# save if we got any
		if ( length(hitType)) {
			smlDF <- data.frame( "Domain.Type"=hitType, "Start"=hitStart, "Stop"=hitStop, 
					"Sequence"=hitSeq, stringsAsFactors=F)
			ord <- order( smlDF$Start)
			smlDF <- smlDF[ ord, ]
			rownames(smlDF) <- 1:nrow(smlDF)
			outDF[[k]] <- smlDF
		}
	}
	if( ! is.null(names(pepSet))) {
		names(outDF) <- names(pepSet)
	}
	if ( makePlot) {
		plot( 1,1, type="n", main="Hydropathy Score", xlab="Amino Acid position", 
			ylab="(hydrophilic)  Kyte-Doolittle Score   (hydrophobic)",
			xlim=c(1,bigNAA), ylim=c(-40,50))
		lines( c(0,bigNAA), c(0,0), lwd=2, lty=2, col='grey40')
		text( bigNAA/2, -40, paste( "Moving Average window = ",window,"aa",sep=""), cex=0.75) 
		for ( k in 1:N) {
			v1 <- tmp1[[k]]
			v2 <- tmp2[[k]]
			lines( 1:length(v2), v2, col=col[k], lty=1, lwd=lwd[k])
			domains <- outDF[[k]]
			if (is.null(domains)) next
			for (i in 1:nrow(domains)) {
				smlStart <- domains$Start[i]
				smlStop <- domains$Stop[i]
				if (domains$Domain.Type[i] == "hydrophobic") {
					smlValue <- max(v2[smlStart:smlStop])
					dy <- 4
				} else {
					smlValue <- min(v2[smlStart:smlStop])
					dy <- -4
				}
				lines( c(smlStart,smlStop), rep.int( smlValue+dy,2), lwd=lwd[k]*2, col=col[k], lty=1)
			}
		}
		if( ! is.null(names(pepSet)) && N < 10) legend( 'topright', names(pepSet), col=col, lwd=lwd, bg='white')
		dev.flush()
	}
	return(outDF)
}

